"""
pubmed_parser.py
────────────────
PubMed PMID 목록을 받아 abstract + PDF 링크를 수집하고
RAG에 바로 사용 가능한 JSON으로 저장합니다.

수집 항목:
  • pmid, doi, pmc_id
  • title, authors, journal, year, pub_date
  • abstract            ← structured abstract도 섹션별로 병합
  • abstract_sections   ← {"BACKGROUND": "...", "METHODS": ...} (structured만)
  • pdf_links           ← [{source, url, type}] 우선순위 순
  • pubmed_url
  • full_text_available ← PMC free full text 여부

입력 방식:
  1) --pmids 32521135 33583502 ...   직접 지정
  2) --input medgen_results.json      MedGen 파서 출력 JSON에서 PMID 자동 수집
  3) --pmid-file pmids.txt            줄바꿈 구분 텍스트 파일

출력:
  JSON 배열: [{...article...}, ...]
  각 article이 RAG chunk로 사용 가능한 "rag_text" 필드 포함

Usage:
  python pubmed_parser.py --pmids 32521135 33583502 --output abstracts.json
  python pubmed_parser.py --input medgen_results.json --output abstracts.json
  python pubmed_parser.py --pmid-file pmids.txt --email you@example.com
"""

import re
import sys
import json
import time
import argparse
import xml.etree.ElementTree as ET
from typing import Any

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import traceback
from datetime import datetime

from utils import get_default_log_path, setup_tee_logging

# ── 상수 ──────────────────────────────────────────────────────
EUTILS_BASE   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
PUBMED_BASE   = "https://pubmed.ncbi.nlm.nih.gov"
PMC_BASE      = "https://www.ncbi.nlm.nih.gov/pmc/articles"
UNPAYWALL_BASE = "https://api.unpaywall.org/v2"

EFETCH_BATCH  = 200   # EFetch 한 번에 처리할 PMID 수
RATE_DELAY    = 0.34  # NCBI 권장 3req/s (API key 없을 때)

# ── HTTP 세션 ─────────────────────────────────────────────────
def _make_session(email: str) -> requests.Session:
    session = requests.Session()
    retry = Retry(
        total=5, backoff_factor=1.5,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    session.headers.update({
        "User-Agent": f"pubmed_parser/1.0 (mailto:{email})",
        "Accept": "application/json, text/xml, */*",
    })
    return session

SESSION: requests.Session | None = None

def _get(url: str, params: dict | None = None, delay: float = RATE_DELAY,
         **kwargs) -> requests.Response:
    time.sleep(delay)
    resp = SESSION.get(url, params=params, timeout=20, **kwargs)
    resp.raise_for_status()
    return resp


# ── PMID 수집 유틸 ────────────────────────────────────────────

def collect_pmids_from_medgen(path: str) -> list[str]:
    """MedGen 파서 출력 JSON에서 PMID 전체 수집 (중복 제거, 순서 유지)"""
    with open(path, encoding="utf-8") as f:
        records = json.load(f)

    seen: set[str] = set()
    pmids: list[str] = []

    def _add(pmid: str):
        pmid = str(pmid).strip()
        if pmid and pmid not in seen:
            seen.add(pmid)
            pmids.append(pmid)

    for rec in records:
        lit = rec.get("literature", {})

        # professional_guidelines
        for paper in lit.get("professional_guidelines", {}).get("papers", []):
            _add(paper.get("pmid", ""))

        # recent_clinical_studies (카테고리별)
        for cat_data in lit.get("recent_clinical_studies", {}).values():
            for paper in cat_data.get("papers", []):
                _add(paper.get("pmid", ""))

        # recent_systematic_reviews
        for paper in lit.get("recent_systematic_reviews", {}).get("papers", []):
            _add(paper.get("pmid", ""))

        # elink_related fallback
        for paper in lit.get("elink_related", []):
            _add(paper.get("pmid", ""))

    print(f"[MedGen JSON] {len(records)}개 레코드에서 PMID {len(pmids)}개 수집")
    return pmids


def collect_pmids_from_file(path: str) -> list[str]:
    with open(path, encoding="utf-8") as f:
        lines = f.read().splitlines()
    pmids = [l.strip() for l in lines if l.strip() and re.match(r"^\d+$", l.strip())]
    print(f"[파일] {path} → PMID {len(pmids)}개")
    return pmids


# ── PubMed EFetch XML 파서 ────────────────────────────────────

def efetch_pubmed_xml(pmids: list[str], email: str) -> str:
    """EFetch로 PubMed XML 가져오기 (배치)"""
    r = _get(
        f"{EUTILS_BASE}/efetch.fcgi",
        params={
            "db": "pubmed",
            "id": ",".join(pmids),
            "rettype": "abstract",
            "retmode": "xml",
            "email": email,
        },
        delay=RATE_DELAY,
    )
    return r.text


def _parse_abstract(article_el: ET.Element) -> tuple[str, dict[str, str]]:
    """
    Abstract 파싱.
    - 일반 abstract: <AbstractText>plain text</AbstractText>
    - Structured:    <AbstractText Label="BACKGROUND" NlmCategory="BACKGROUND">...</AbstractText>

    반환: (전체_텍스트, 섹션별_dict)
    """
    abstract_el = article_el.find(".//Abstract")
    if abstract_el is None:
        return "", {}

    sections: dict[str, str] = {}
    plain_parts: list[str] = []

    for at in abstract_el.findall("AbstractText"):
        label = at.get("Label") or at.get("NlmCategory") or ""
        # ET.tostring으로 내부 태그(sup, i, b 등) 텍스트 포함 추출
        inner = "".join(at.itertext()).strip()
        inner = re.sub(r"\s+", " ", inner)

        if label:
            sections[label] = inner
            plain_parts.append(f"{label}: {inner}")
        else:
            plain_parts.append(inner)

    full_text = " ".join(plain_parts).strip()
    return full_text, sections


def _parse_authors(article_el: ET.Element) -> list[str]:
    authors = []
    for author in article_el.findall(".//AuthorList/Author"):
        last  = (author.findtext("LastName") or "").strip()
        fore  = (author.findtext("ForeName") or author.findtext("Initials") or "").strip()
        cname = (author.findtext("CollectiveName") or "").strip()
        if cname:
            authors.append(cname)
        elif last:
            authors.append(f"{last} {fore}".strip())
    return authors


def _parse_pub_date(article_el: ET.Element) -> tuple[str, str]:
    """(year, pub_date_str) 반환"""
    # PubDate 우선 탐색
    for pd in article_el.findall(".//PubDate"):
        year  = pd.findtext("Year") or ""
        month = pd.findtext("Month") or ""
        day   = pd.findtext("Day") or ""
        med   = pd.findtext("MedlineDate") or ""
        if year:
            return year, f"{year} {month} {day}".strip()
        if med:
            m = re.match(r"(\d{4})", med)
            return (m.group(1) if m else ""), med
    return "", ""


def _parse_ids(pubmed_data_el: ET.Element) -> dict[str, str]:
    """ArticleIdList에서 pubmed/pmc/doi/mid 추출"""
    ids: dict[str, str] = {}
    for aid in pubmed_data_el.findall(".//ArticleIdList/ArticleId"):
        id_type = aid.get("IdType", "")
        val = (aid.text or "").strip()
        if id_type and val:
            ids[id_type] = val
    return ids


def _build_pdf_links(ids: dict[str, str], pmc_free: bool) -> list[dict[str, str]]:
    """
    PDF 링크 우선순위:
    1. PMC free full text PDF (open access, 가장 확실)
    2. PMC full text HTML → PDF 변환 URL
    3. DOI → Unpaywall open access PDF (별도 조회 필요 → url만 제공)
    4. PubMed abstract 페이지 (PDF 없어도 참조용)
    """
    links: list[dict[str, str]] = []
    pmid = ids.get("pubmed", "")
    pmc  = ids.get("pmc", "").replace("PMC", "").strip()
    doi  = ids.get("doi", "")

    # 1. PMC PDF (free full text)
    if pmc:
        links.append({
            "source": "PMC",
            "url": f"{PMC_BASE}/PMC{pmc}/pdf/",
            "type": "pdf",
            "note": "PMC free full text PDF",
        })
        links.append({
            "source": "PMC",
            "url": f"{PMC_BASE}/PMC{pmc}/",
            "type": "html",
            "note": "PMC full text HTML",
        })

    # 2. DOI 링크 (출판사 페이지, 유료일 수 있음)
    if doi:
        links.append({
            "source": "DOI",
            "url": f"https://doi.org/{doi}",
            "type": "publisher",
            "note": "Publisher page (may require subscription)",
        })
        # Unpaywall open access 조회 URL (런타임에 GET 필요)
        links.append({
            "source": "Unpaywall",
            "url": f"{UNPAYWALL_BASE}/{doi}?email=user@example.com",
            "type": "oa_check",
            "note": "Check Unpaywall for open access PDF",
        })

    # 3. PubMed 페이지
    if pmid:
        links.append({
            "source": "PubMed",
            "url": f"{PUBMED_BASE}/{pmid}/",
            "type": "abstract",
            "note": "PubMed abstract page",
        })

    return links


def _check_unpaywall(doi: str, email: str) -> dict | None:
    """
    Unpaywall API로 open access PDF URL 조회.
    실패 시 None 반환 (네트워크 오류 무시).
    """
    if not doi:
        return None
    try:
        r = _get(
            f"{UNPAYWALL_BASE}/{doi}",
            params={"email": email},
            delay=0.5,
        )
        data = r.json()
        oa_url = data.get("best_oa_location", {})
        pdf_url = oa_url.get("url_for_pdf") or oa_url.get("url")
        if pdf_url:
            return {
                "source": "Unpaywall",
                "url": pdf_url,
                "type": "pdf",
                "note": f"Open access PDF via Unpaywall (license: {data.get('license', 'unknown')})",
            }
    except Exception:
        pass
    return None


def parse_pubmed_xml(xml_text: str, email: str,
                     use_unpaywall: bool = True) -> list[dict[str, Any]]:
    """
    EFetch XML 전체를 파싱하여 article 목록 반환.
    """
    try:
        root = ET.fromstring(xml_text.encode("utf-8"))
    except ET.ParseError as e:
        print(f"[Error] XML 파싱 실패: {e}")
        return []

    articles: list[dict[str, Any]] = []

    for pubmed_article in root.findall(".//PubmedArticle"):
        medline = pubmed_article.find("MedlineCitation")
        pubmed_data = pubmed_article.find("PubmedData")

        if medline is None:
            continue

        # PMID
        pmid_el = medline.find("PMID")
        pmid = pmid_el.text.strip() if pmid_el is not None and pmid_el.text else ""

        article_el = medline.find("Article")
        if article_el is None:
            articles.append({"pmid": pmid, "error": "No Article element"})
            continue

        # 제목
        title_el = article_el.find("ArticleTitle")
        title = "".join(title_el.itertext()).strip() if title_el is not None else ""
        title = re.sub(r"\s+", " ", title)

        # Abstract
        abstract, abstract_sections = _parse_abstract(article_el)

        # 저자
        authors = _parse_authors(article_el)

        # 저널
        journal = ""
        journal_el = article_el.find(".//Journal/Title")
        if journal_el is not None:
            journal = (journal_el.text or "").strip()
        # ISOAbbreviation fallback
        if not journal:
            iso = article_el.find(".//Journal/ISOAbbreviation")
            if iso is not None:
                journal = (iso.text or "").strip()

        # 발행일
        year, pub_date = _parse_pub_date(article_el)

        # DOI (Article 레벨)
        doi = ""
        for loc in article_el.findall("ELocationID"):
            if loc.get("EIdType") == "doi":
                doi = (loc.text or "").strip()
                break

        # ArticleId (PMC, doi 보완)
        ids: dict[str, str] = {"pubmed": pmid}
        if pubmed_data is not None:
            ids.update(_parse_ids(pubmed_data))
        if doi:
            ids["doi"] = doi  # Article 레벨이 우선

        pmc_id = ids.get("pmc", "").replace("PMC", "").strip()

        # MeSH 키워드
        mesh_terms: list[str] = []
        for mh in medline.findall(".//MeshHeadingList/MeshHeading"):
            desc = mh.find("DescriptorName")
            if desc is not None and desc.text:
                mesh_terms.append(desc.text.strip())

        # KeywordList
        keywords: list[str] = []
        for kw in medline.findall(".//KeywordList/Keyword"):
            if kw.text:
                keywords.append(kw.text.strip())

        # Publication types
        pub_types: list[str] = []
        for pt in article_el.findall(".//PublicationTypeList/PublicationType"):
            if pt.text:
                pub_types.append(pt.text.strip())

        # PDF 링크 빌드
        pdf_links = _build_pdf_links(ids, bool(pmc_id))

        # Unpaywall open access PDF 조회 (doi 있고 PMC 없을 때)
        if use_unpaywall and ids.get("doi") and not pmc_id:
            oa = _check_unpaywall(ids["doi"], email)
            if oa:
                pdf_links.insert(0, oa)  # 가장 앞에 삽입

        # RAG용 텍스트 (title + abstract 전문)
        rag_text = _build_rag_text(
            pmid=pmid, title=title, authors=authors, journal=journal,
            year=year, abstract=abstract, mesh_terms=mesh_terms,
            keywords=keywords,
        )

        article: dict[str, Any] = {
            "pmid":               pmid,
            "doi":                ids.get("doi", ""),
            "pmc_id":             f"PMC{pmc_id}" if pmc_id else "",
            "title":              title,
            "authors":            authors,
            "journal":            journal,
            "year":               year,
            "pub_date":           pub_date,
            "publication_types":  pub_types,
            "abstract":           abstract,
            "abstract_sections":  abstract_sections,   # structured abstract 섹션별
            "mesh_terms":         mesh_terms,
            "keywords":           keywords,
            "full_text_available": bool(pmc_id),
            "pdf_links":          pdf_links,
            "pubmed_url":         f"{PUBMED_BASE}/{pmid}/",
            "rag_text":           rag_text,
        }
        articles.append(article)

    return articles


def _build_rag_text(
    pmid: str, title: str, authors: list[str], journal: str,
    year: str, abstract: str, mesh_terms: list[str], keywords: list[str],
) -> str:
    """
    RAG 청크 텍스트.
    벡터 DB 임베딩에 적합하도록 메타데이터 + 본문을 단일 문자열로 구성.
    """
    author_str = ", ".join(authors[:6])
    if len(authors) > 6:
        author_str += f" et al. ({len(authors)} authors)"

    parts = [
        f"Title: {title}",
        f"Authors: {author_str}",
        f"Journal: {journal} ({year})",
        f"PMID: {pmid}",
    ]
    if abstract:
        parts.append(f"Abstract: {abstract}")
    if mesh_terms:
        parts.append(f"MeSH Terms: {'; '.join(mesh_terms[:15])}")
    if keywords:
        parts.append(f"Keywords: {'; '.join(keywords[:10])}")

    return "\n".join(parts)


# ── 배치 처리 ─────────────────────────────────────────────────

def fetch_articles(
    pmids: list[str],
    email: str,
    use_unpaywall: bool = True,
    batch_size: int = EFETCH_BATCH,
) -> list[dict[str, Any]]:
    """PMID 목록을 배치로 나눠 EFetch 후 파싱."""
    all_articles: list[dict[str, Any]] = []
    total = len(pmids)

    for i in range(0, total, batch_size):
        batch = pmids[i : i + batch_size]
        print(f"  [EFetch] {i+1}~{min(i+batch_size, total)}/{total} ({len(batch)}개)")
        try:
            xml_text = efetch_pubmed_xml(batch, email)
        except Exception as e:
            print(f"  [Error] EFetch 실패: {e}")
            # 오류 시 개별 PMID 기록
            for pmid in batch:
                all_articles.append({"pmid": pmid, "error": str(e)})
            continue

        parsed = parse_pubmed_xml(xml_text, email, use_unpaywall=use_unpaywall)
        all_articles.extend(parsed)

        # 파싱된 PMID 확인
        parsed_pmids = {a.get("pmid") for a in parsed}
        missing = set(batch) - parsed_pmids
        for pmid in missing:
            all_articles.append({"pmid": pmid, "error": "Not found in EFetch response"})

    return all_articles


# ── 통계 출력 ─────────────────────────────────────────────────

def print_summary(articles: list[dict[str, Any]]) -> None:
    total      = len(articles)
    with_abs   = sum(1 for a in articles if a.get("abstract"))
    with_pmc   = sum(1 for a in articles if a.get("pmc_id"))
    errors     = sum(1 for a in articles if a.get("error"))
    with_oa    = sum(
        1 for a in articles
        if any(l.get("type") in ("pdf", "html") and l.get("source") != "PubMed"
               for l in a.get("pdf_links", []))
    )

    print(f"\n{'='*50}")
    print(f"  총 PMID        : {total}")
    print(f"  Abstract 수집  : {with_abs} ({with_abs/total*100:.1f}%)" if total else "")
    print(f"  PMC full text  : {with_pmc} ({with_pmc/total*100:.1f}%)" if total else "")
    print(f"  OA/PDF 링크    : {with_oa}")
    print(f"  오류           : {errors}")
    print(f"{'='*50}")


# ── CLI ───────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="PubMed PMID → Abstract + PDF 링크 수집기 (RAG용)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
예시:
  # PMID 직접 지정
  python pubmed_parser.py --pmids 32521135 33583502 --output abstracts.json

  # MedGen 출력 JSON에서 자동 수집
  python pubmed_parser.py --input medgen_results.json --output abstracts.json

  # 텍스트 파일 (줄당 PMID 1개)
  python pubmed_parser.py --pmid-file pmids.txt --output abstracts.json

  # Unpaywall 조회 비활성화 (속도 우선)
  python pubmed_parser.py --pmids 32521135 --no-unpaywall --output abstracts.json
        """,
    )
    parser.add_argument("--pmids",       nargs="+",          help="PMID 목록")
    parser.add_argument("--input",       type=str,           help="MedGen 파서 출력 JSON 경로")
    parser.add_argument("--pmid-file",   type=str,           help="PMID 목록 텍스트 파일")
    parser.add_argument("--output",      type=str, default="pubmed_abstracts.json")
    parser.add_argument("--email",       type=str, default="your.email@example.com",
                        help="NCBI API 이메일")
    parser.add_argument("--no-unpaywall", action="store_true",
                        help="Unpaywall open access PDF 조회 비활성화")
    parser.add_argument("--batch-size",  type=int, default=EFETCH_BATCH)
    parser.add_argument(
        "--log",
        type=str,
        default="",
        help="로그 파일 경로. 미지정 시 output 파일명을 기준으로 생성",
    )
    args = parser.parse_args()

    # PMID 수집
    pmids: list[str] = []
    if args.pmids:
        pmids = [p.strip() for p in args.pmids if p.strip().isdigit()]
        print(f"[CLI] PMID {len(pmids)}개")
    elif args.input:
        pmids = collect_pmids_from_medgen(args.input)
    elif getattr(args, "pmid_file", None):
        pmids = collect_pmids_from_file(args.pmid_file)
    else:
        parser.error("--pmids, --input, --pmid-file 중 하나는 필수입니다.")

    if not pmids:
        print("[Warning] 수집된 PMID 없음. 종료.")
        sys.exit(0)

    # 세션 초기화
    global SESSION
    SESSION = _make_session(args.email)

    log_path = args.log or get_default_log_path(args.output)
    tee_logger = setup_tee_logging(log_path)
    
    try:
        print(f"\n▶ PubMed 수집 시작: {len(pmids)}개 PMID")
        articles = fetch_articles(
            pmids=pmids,
            email=args.email,
            use_unpaywall=not args.no_unpaywall,
            batch_size=args.batch_size,
        )

        print_summary(articles)

        with open(args.output, "w", encoding="utf-8") as f:
            json.dump(articles, f, ensure_ascii=False, indent=2)
        print(f"\n▶ 저장 완료: {args.output}")

        # 샘플 미리보기
        for a in articles[:2]:
            if a.get("error"):
                print(f"\n  [오류] PMID={a['pmid']}: {a['error']}")
                continue
            print(f"\n  PMID={a['pmid']} | {a['title'][:60]}")
            print(f"  Abstract: {a['abstract'][:100]}...")
            print(f"  PDF links ({len(a['pdf_links'])}개):")
            for l in a['pdf_links'][:3]:
                print(f"    [{l['source']}][{l['type']}] {l['url'][:80]}")

    except Exception:
        print("\n[ERROR] 실행 중 오류가 발생했습니다.")
        traceback.print_exc()
        raise

    finally:
        tee_logger.close()

if __name__ == "__main__":
    main()