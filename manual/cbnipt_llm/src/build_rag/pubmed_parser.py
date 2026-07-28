"""
pubmed_parser.py
────────────────
PMID 목록 → abstract + PDF 링크 수집 + PMC ID 검증 + YAML 저장

주요 기능:
  1. EFetch XML 파싱 (abstract, 저자, 저널, MeSH, DOI, PMC ID)
  2. PMC ID 역검증 (PMC → PMID 역조회로 mismatch 탐지)
  3. JSON 저장 (pmid keyed dict)
  4. YAML 저장 (pmid keyed, review 공란 포함)

PMC ID 검증 이유:
  NCBI EFetch의 ArticleIdList에서 pmc 값이
  Author Manuscript 버전이거나 메타데이터 오류인 경우
  실제 PMID와 다른 논문의 PMC 페이지를 가리킬 수 있음.
  역조회로 PMID 일치 여부 확인 후 불일치 시 pdf_links에서 PMC 제외.

Usage:
  python pubmed_parser.py --pmids 38745141 33053283 --output out/
  python pubmed_parser.py --pmid-file pmids.txt --output out/ --no-verify-pmc
"""

import re
import sys
import json
import time
import argparse
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

import requests
import yaml
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# ── 상수 ──────────────────────────────────────────────────────
EUTILS  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
PUBMED  = "https://pubmed.ncbi.nlm.nih.gov"
PMC_BASE = "https://www.ncbi.nlm.nih.gov/pmc/articles"

SESSION: requests.Session | None = None


def _make_session(email: str) -> requests.Session:
    s = requests.Session()
    retry = Retry(total=5, backoff_factor=1.5,
                  status_forcelist=[429, 500, 502, 503, 504],
                  allowed_methods=["GET"])
    s.mount("https://", HTTPAdapter(max_retries=retry))
    s.mount("http://",  HTTPAdapter(max_retries=retry))
    s.headers.update({"User-Agent": f"pubmed-parser/1.0 (mailto:{email})"})
    return s


def _get(url: str, params: dict | None = None,
         delay: float = 0.35) -> requests.Response:
    time.sleep(delay)
    r = SESSION.get(url, params=params, timeout=20)
    r.raise_for_status()
    return r


# ════════════════════════════════════════════════════════════════
# PMC ID 검증
# ════════════════════════════════════════════════════════════════

# 검증 캐시 (같은 PMC ID 중복 조회 방지)
_pmc_verify_cache: dict[str, str] = {}  # pmc_id → verified_pmid or ""


def verify_pmc_id(pmc_id: str, expected_pmid: str, email: str) -> dict:
    """
    PMC ID를 역조회하여 해당 PMID와 일치하는지 검증.

    NCBI EFetch의 ArticleIdList에서 pmc 값이
    Author Manuscript(AM) 또는 메타데이터 오류인 경우
    실제 논문의 PMC ID와 다를 수 있음.

    반환:
      {
        "pmc_id":       원본 PMC ID,
        "verified":     True/False,
        "matched_pmid": 역조회로 찾은 PMID (없으면 ""),
        "pmc_type":     "publisher" | "author_manuscript" | "mismatch" | "not_found",
        "note":         설명 문자열
      }
    """
    pmc_num = pmc_id.replace("PMC", "").strip()
    cache_key = pmc_num

    if cache_key in _pmc_verify_cache:
        matched = _pmc_verify_cache[cache_key]
    else:
        try:
            r = _get(f"{EUTILS}/esearch.fcgi",
                     params={"db": "pubmed",
                             "term": f"PMC{pmc_num}[pmc]",
                             "retmode": "json",
                             "email": email},
                     delay=0.35)
            id_list = r.json().get("esearchresult", {}).get("idlist", [])
            matched = id_list[0] if id_list else ""
            _pmc_verify_cache[cache_key] = matched
        except Exception as e:
            print(f"    [PMC verify] PMC{pmc_num} 역조회 실패: {e}")
            return {
                "pmc_id": pmc_id, "verified": False,
                "matched_pmid": "", "pmc_type": "not_found",
                "note": f"역조회 실패: {e}",
            }

    if not matched:
        return {
            "pmc_id": pmc_id, "verified": False,
            "matched_pmid": "", "pmc_type": "not_found",
            "note": f"PMC{pmc_num}에 해당하는 PMID 없음",
        }

    if matched == str(expected_pmid):
        return {
            "pmc_id": pmc_id, "verified": True,
            "matched_pmid": matched, "pmc_type": "publisher",
            "note": "PMC ID 검증 통과 (PMID 일치)",
        }
    else:
        # PMID 불일치 → AM 버전이거나 메타데이터 오류
        return {
            "pmc_id": pmc_id, "verified": False,
            "matched_pmid": matched,
            "pmc_type": "author_manuscript_or_mismatch",
            "note": (f"PMC{pmc_num}의 역조회 PMID={matched}가 "
                     f"현재 PMID={expected_pmid}와 불일치. "
                     f"Author Manuscript 또는 메타데이터 오류 가능성."),
        }


def _build_pdf_links(pmid: str, doi: str, pmc_id: str,
                     pmc_verify: dict | None) -> tuple[list[dict], bool]:
    """
    PDF 링크 목록 생성.
    pmc_verify 결과에 따라 PMC 링크 포함/제외 결정.

    반환: (pdf_links, full_text_available)
    """
    links: list[dict] = []
    full_text = False

    if pmc_id:
        if pmc_verify and not pmc_verify["verified"]:
            # PMC ID 불일치 → 링크는 남기되 경고 note 추가
            links.append({
                "source": "PMC",
                "type":   "pdf",
                "url":    f"{PMC_BASE}/{pmc_id}/pdf/",
                "note":   f"⚠ {pmc_verify['note']}",
                "verified": False,
            })
            links.append({
                "source": "PMC",
                "type":   "html",
                "url":    f"{PMC_BASE}/{pmc_id}/",
                "note":   f"⚠ {pmc_verify['note']}",
                "verified": False,
            })
            # full_text_available은 False (내용 불일치 가능성)
        else:
            # 검증 통과 또는 검증 미실시
            links.append({
                "source":   "PMC",
                "type":     "pdf",
                "url":      f"{PMC_BASE}/{pmc_id}/pdf/",
                "note":     "PMC free full text PDF",
                "verified": True if pmc_verify else None,
            })
            links.append({
                "source":   "PMC",
                "type":     "html",
                "url":      f"{PMC_BASE}/{pmc_id}/",
                "note":     "PMC full text HTML",
                "verified": True if pmc_verify else None,
            })
            full_text = True

    if doi:
        links.append({
            "source": "DOI",
            "type":   "publisher",
            "url":    f"https://doi.org/{doi}",
            "note":   "Publisher page (may require subscription)",
            "verified": None,
        })

    links.append({
        "source":   "PubMed",
        "type":     "abstract",
        "url":      f"{PUBMED}/{pmid}/",
        "note":     "PubMed abstract page",
        "verified": None,
    })

    return links, full_text


# ════════════════════════════════════════════════════════════════
# EFetch XML 파싱
# ════════════════════════════════════════════════════════════════

def _parse_abstract(article_el: ET.Element) -> tuple[str, dict]:
    """Abstract → (full_text, sections_dict)"""
    abs_el = article_el.find(".//Abstract")
    if abs_el is None:
        return "", {}
    sections: dict[str, str] = {}
    parts: list[str] = []
    for at in abs_el.findall("AbstractText"):
        label = at.get("Label") or at.get("NlmCategory") or ""
        text  = re.sub(r"\s+", " ", "".join(at.itertext()).strip())
        if label:
            sections[label] = text
            parts.append(f"{label}: {text}")
        else:
            parts.append(text)
    return " ".join(parts), sections


def _parse_authors(article_el: ET.Element) -> list[str]:
    authors = []
    for au in article_el.findall(".//AuthorList/Author"):
        last  = (au.findtext("LastName") or "").strip()
        fore  = (au.findtext("ForeName") or
                 au.findtext("Initials") or "").strip()
        cname = (au.findtext("CollectiveName") or "").strip()
        if cname:
            authors.append(cname)
        elif last:
            authors.append(f"{last} {fore}".strip())
    return authors


def fetch_abstracts(pmids: list[str], email: str,
                    verify_pmc: bool = True,
                    batch_size: int = 100) -> dict[str, dict]:
    """
    PMID 목록 → {pmid: article_dict}

    article_dict 필드:
      pmid, doi, pmc_id, pmc_verify
      title, authors, journal, year, pub_date
      publication_types
      abstract, abstract_sections
      mesh_terms, keywords
      full_text_available
      pdf_links (PMC 검증 결과 반영)
      pubmed_url
    """
    if not pmids:
        return {}

    lookup: dict[str, dict] = {}
    total = len(pmids)

    for i in range(0, total, batch_size):
        batch = pmids[i: i + batch_size]
        print(f"  [EFetch] {i+1}~{min(i+batch_size, total)}/{total}")

        try:
            r = _get(f"{EUTILS}/efetch.fcgi",
                     params={"db": "pubmed", "id": ",".join(batch),
                             "rettype": "abstract", "retmode": "xml",
                             "email": email})
            root = ET.fromstring(r.text.encode("utf-8"))
        except Exception as e:
            print(f"  [EFetch] 오류: {e}")
            for pmid in batch:
                lookup[pmid] = {"pmid": pmid, "error": str(e)}
            continue

        for pa in root.findall(".//PubmedArticle"):
            mc  = pa.find("MedlineCitation")
            art = mc.find("Article") if mc is not None else None
            if mc is None or art is None:
                continue

            pmid_el = mc.find("PMID")
            pmid = pmid_el.text.strip() if pmid_el is not None else ""

            # 제목
            title_el = art.find("ArticleTitle")
            title = re.sub(r"\s+", " ",
                           "".join(title_el.itertext()).strip()
                           if title_el is not None else "")

            # Abstract
            abstract, abs_sections = _parse_abstract(art)

            # 저자
            authors = _parse_authors(art)

            # 저널
            journal = (art.findtext(".//Journal/Title") or
                       art.findtext(".//Journal/ISOAbbreviation") or "")

            # 발행일
            pub_date_el = art.find(".//PubDate")
            year, pub_date = "", ""
            if pub_date_el is not None:
                year = pub_date_el.findtext("Year") or ""
                if not year:
                    med = pub_date_el.findtext("MedlineDate") or ""
                    m = re.match(r"(\d{4})", med)
                    year = m.group(1) if m else ""
                mo  = pub_date_el.findtext("Month") or ""
                day = pub_date_el.findtext("Day") or ""
                pub_date = f"{year} {mo} {day}".strip()

            # 출판 타입
            pub_types = [
                pt.text.strip()
                for pt in art.findall(".//PublicationTypeList/PublicationType")
                if pt.text
            ]

            # DOI + PMC ID
            doi, pmc_id = "", ""
            pd_el = pa.find("PubmedData")
            if pd_el is not None:
                for aid in pd_el.findall(".//ArticleIdList/ArticleId"):
                    t = aid.get("IdType", "")
                    v = (aid.text or "").strip()
                    if t == "doi" and not doi:
                        doi = v
                    elif t == "pmc" and not pmc_id:
                        pmc_id = v
            # Article 레벨 DOI 우선
            for loc in art.findall("ELocationID"):
                if loc.get("EIdType") == "doi":
                    doi = (loc.text or "").strip()
                    break

            # MeSH
            mesh_terms = [
                mh.findtext("DescriptorName") or ""
                for mh in mc.findall(".//MeshHeadingList/MeshHeading")
            ]
            mesh_terms = [m for m in mesh_terms if m]

            # Keyword
            keywords = [
                kw.text.strip()
                for kw in mc.findall(".//KeywordList/Keyword")
                if kw.text
            ]

            # PMC ID 검증
            pmc_verify: dict | None = None
            if pmc_id and verify_pmc:
                print(f"    [PMC verify] PMID={pmid} PMC={pmc_id}")
                pmc_verify = verify_pmc_id(pmc_id, pmid, email)
                status = "✓" if pmc_verify["verified"] else "✗"
                print(f"    {status} {pmc_verify['note'][:70]}")

            # PDF 링크 빌드
            pdf_links, full_text = _build_pdf_links(
                pmid, doi, pmc_id, pmc_verify
            )

            lookup[pmid] = {
                "pmid":              pmid,
                "doi":               doi,
                "pmc_id":            pmc_id,
                "pmc_verify":        pmc_verify,
                "title":             title,
                "authors":           authors,
                "journal":           journal,
                "year":              year,
                "pub_date":          pub_date,
                "publication_types": pub_types,
                "abstract":          abstract,
                "abstract_sections": abs_sections,
                "mesh_terms":        mesh_terms,
                "keywords":          keywords,
                "full_text_available": full_text,
                "pdf_links":         pdf_links,
                "pubmed_url":        f"{PUBMED}/{pmid}/",
            }

    return lookup


# ════════════════════════════════════════════════════════════════
# YAML 저장
# ════════════════════════════════════════════════════════════════

def _best_pdf_url(pdf_links: list[dict]) -> str:
    """verified=True인 PDF URL 우선, 없으면 publisher, 없으면 PubMed"""
    order = {"pdf": 0, "html": 1, "publisher": 2, "abstract": 3}
    verified = [l for l in pdf_links if l.get("verified") is True]
    if verified:
        verified.sort(key=lambda x: order.get(x.get("type", ""), 9))
        return verified[0]["url"]
    for t in ("publisher", "abstract"):
        for l in pdf_links:
            if l.get("type") == t:
                return l["url"]
    return ""


def to_yaml_dict(article: dict) -> dict:
    """
    article dict → YAML 저장 형식.
    메타데이터 + review 공란 포함.
    """
    pmc_v    = article.get("pmc_verify") or {}
    pmc_note = pmc_v.get("note", "") if pmc_v else ""
    pmc_ok   = pmc_v.get("verified", None) if pmc_v else None

    # PDF 링크 정리 (verified 여부 명시)
    pdf_entries = []
    for l in article.get("pdf_links", []):
        entry = {"source": l["source"], "type": l["type"], "url": l["url"]}
        if l.get("verified") is not None:
            entry["verified"] = l["verified"]
        if "⚠" in l.get("note", ""):
            entry["warning"] = l["note"].replace("⚠ ", "")
        pdf_entries.append(entry)

    return {
        # ── 서지 정보 ──────────────────────────────────────────
        "pmid":    article.get("pmid", ""),
        "doi":     article.get("doi", ""),
        "pmc_id":  article.get("pmc_id", ""),
        "title":   article.get("title", ""),
        "authors": article.get("authors", []),
        "journal": article.get("journal", ""),
        "year":    article.get("year", ""),
        "pub_date":article.get("pub_date", ""),
        "publication_types": article.get("publication_types", []),
        "pubmed_url": article.get("pubmed_url", ""),

        # ── Full text / PDF ────────────────────────────────────
        "full_text_available": article.get("full_text_available", False),
        "pmc_verified": pmc_ok,
        "pmc_verify_note": pmc_note or None,
        "pdf_links": pdf_entries,
        "best_pdf_url": _best_pdf_url(article.get("pdf_links", [])),

        # ── 내용 ──────────────────────────────────────────────
        "abstract": article.get("abstract", ""),
        "abstract_sections": article.get("abstract_sections", {}) or None,
        "mesh_terms": article.get("mesh_terms", []),
        "keywords":   article.get("keywords", []),

        # ── Review 공란 (직접 채워넣는 항목) ───────────────────
        "review": {
            "relevance": None,       # high / medium / low / exclude
            "study_design": None,    # RCT / cohort / case-control / review / meta-analysis / case-report / other
            "sample_size": None,     # 숫자 또는 설명
            "key_findings": None,    # 이 논문의 핵심 결과 (1-3문장)
            "limitations": [],       # 이 논문 방법의 한계
            "apply_to": None,        # |
                                     # 어떤 분석/스크립트에 어떻게 적용할지
                                     # 예) BAF scoring의 bin size 기준 500kb 사용
            "data_availability": {
                "geo_id":  None,     # GEO accession (예: GSE12345)
                "sra_id":  None,     # SRA accession (예: SRP123456)
                "github":  None,     # 코드 저장소 URL
                "other":   None,     # 기타 데이터 접근 경로
            },
            "note": None,            # |
                                     # 자유 메모
                                     # 예) 한국인 코호트 데이터 포함
        },
    }


def save_yaml(articles: dict[str, dict], output_path: str):
    """
    {pmid: article_dict} → pmid keyed YAML 파일.
    각 PMID가 최상위 키.
    """
    yaml_data: dict[str, dict] = {}

    for pmid, article in articles.items():
        if article.get("error"):
            yaml_data[pmid] = {
                "pmid":  pmid,
                "error": article["error"],
                "review": {"relevance": None, "note": None},
            }
        else:
            yaml_data[pmid] = to_yaml_dict(article)

    with open(output_path, "w", encoding="utf-8") as f:
        # 헤더 주석
        f.write("# PubMed Article Records\n")
        f.write("# Generated by pubmed_parser.py\n")
        f.write("# \n")
        f.write("# review 항목을 직접 채워 프로젝트 노트로 활용하세요.\n")
        f.write("# relevance: high / medium / low / exclude\n")
        f.write("# study_design: RCT / cohort / case-control / review / meta-analysis / case-report\n")
        f.write("#\n\n")

        yaml.dump(
            yaml_data,
            f,
            allow_unicode=True,
            default_flow_style=False,
            sort_keys=False,
            indent=2,
            width=120,
        )

    print(f"  [YAML] 저장: {output_path} ({len(yaml_data)}개 PMID)")


# ════════════════════════════════════════════════════════════════
# CLI
# ════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="PubMed PMID → abstract + PDF 링크 + PMC 검증 + YAML",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
예시:
  python pubmed_parser.py --pmids 38745141 33053283 --output out/

  # PMC 검증 없이 빠르게
  python pubmed_parser.py --pmids 38745141 --output out/ --no-verify-pmc

  # 파일에서 PMID 목록
  python pubmed_parser.py --pmid-file pmids.txt --output out/
        """
    )
    parser.add_argument("--pmids",         nargs="+", help="PMID 목록")
    parser.add_argument("--pmid-file",     type=str,  help="PMID 파일 (줄당 1개)")
    parser.add_argument("--output",        type=str,  default="output",
                        help="출력 디렉토리")
    parser.add_argument("--prefix",        type=str,  default="pubmed",
                        help="출력 파일 이름 prefix")
    parser.add_argument("--email",         type=str,  default="your@email.com")
    parser.add_argument("--no-verify-pmc", action="store_true",
                        help="PMC ID 역검증 생략 (속도 우선)")
    parser.add_argument("--batch-size",    type=int,  default=100)
    args = parser.parse_args()

    # PMID 수집
    pmids: list[str] = []
    if args.pmids:
        pmids = [p.strip() for p in args.pmids if p.strip().isdigit()]
    elif args.pmid_file:
        with open(args.pmid_file, encoding="utf-8") as f:
            pmids = [l.strip() for l in f if l.strip().isdigit()]
    else:
        parser.error("--pmids 또는 --pmid-file 필수")

    print(f"▶ PMID {len(pmids)}개 처리 시작")

    global SESSION
    SESSION = _make_session(args.email)

    Path(args.output).mkdir(parents=True, exist_ok=True)

    articles = fetch_abstracts(
        pmids=pmids,
        email=args.email,
        verify_pmc=not args.no_verify_pmc,
        batch_size=args.batch_size,
    )

    # JSON 저장 (pmid keyed)
    json_path = Path(args.output) / f"{args.prefix}_abstracts.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(articles, f, ensure_ascii=False, indent=2)
    print(f"  [JSON] 저장: {json_path}")

    # YAML 저장 (pmid keyed + review 공란)
    yaml_path = Path(args.output) / f"{args.prefix}_abstracts.yaml"
    save_yaml(articles, str(yaml_path))

    # 요약
    verified_ok  = sum(1 for a in articles.values()
                       if (a.get("pmc_verify") or {}).get("verified") is True)
    verified_fail= sum(1 for a in articles.values()
                       if (a.get("pmc_verify") or {}).get("verified") is False)
    with_abs     = sum(1 for a in articles.values() if a.get("abstract"))
    full_text    = sum(1 for a in articles.values() if a.get("full_text_available"))

    print(f"\n▶ 완료")
    print(f"  총 PMID       : {len(articles)}")
    print(f"  abstract 보유 : {with_abs}")
    print(f"  full text     : {full_text}")
    print(f"  PMC 검증 통과 : {verified_ok}")
    print(f"  PMC 검증 실패 : {verified_fail}  ← pdf_links에 ⚠ 표시")
    print(f"  JSON : {json_path}")
    print(f"  YAML : {yaml_path}")


if __name__ == "__main__":
    main()