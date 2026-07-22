"""
medgen_parser.py
─────────────────────
NCBI MedGen 질병 페이지 파서 (실제 HTML DOM 구조 기반)

수집 항목:
  • 기본 정보  : MedGen UID, 질병명, Concept ID, 타입, definition, synonyms
  • Cross-refs : MONDO, OMIM, Orphanet, SNOMED-CT (medgenTable에서 직접 추출)
  • 논문 섹션 :
      professional_guidelines (ID_105)
      etiology / diagnosis / therapy / prognosis /
      clinical_prediction_guides (ID_103 내 h3.subhead 기준)
      recent_systematic_reviews (ID_104)
    각 카테고리마다:
      papers: [{pmid, title, authors, journal, citation, pmc_free, url}]
      see_all: {count, url}    ← "See all (N)" 링크 처리

Usage:
    python medgen_parser.py --disease "Down syndrome" --max-results 3 --output out.json
    python medgen_parser.py --uids 4385 --output out.json
"""

import re
import sys
import json
import time
import argparse
from typing import Any

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from bs4 import BeautifulSoup, Tag
import traceback
from datetime import datetime

from utils import get_default_log_path, setup_tee_logging

# ── 상수 ──────────────────────────────────────────────────────
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
MEDGEN_BASE = "https://www.ncbi.nlm.nih.gov/medgen"

# h3.subhead 텍스트 → JSON key 매핑 (소문자 부분매칭)
SUBHEAD_MAP = {
    "etiology":                 "etiology",
    "diagnosis":                "diagnosis",
    "therapy":                  "therapy",
    "prognosis":                "prognosis",
    "clinical prediction":      "clinical_prediction_guides",
}

# ── HTTP 세션 ─────────────────────────────────────────────────
def _make_session() -> requests.Session:
    session = requests.Session()
    retry = Retry(
        total=5,
        backoff_factor=1.5,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    session.headers.update({
        "User-Agent": (
            "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
            "AppleWebKit/537.36 (KHTML, like Gecko) "
            "Chrome/124.0.0.0 Safari/537.36"
        ),
        "Accept": "text/html,application/xhtml+xml;q=0.9,*/*;q=0.8",
        "Accept-Language": "en-US,en;q=0.9",
    })
    return session

SESSION = _make_session()

def _get(url: str, params: dict | None = None, delay: float = 0.4) -> requests.Response:
    time.sleep(delay)
    resp = SESSION.get(url, params=params, timeout=20)
    resp.raise_for_status()
    return resp

def _txt(tag: Tag | None, sep: str = " ") -> str:
    if tag is None:
        return ""
    return re.sub(r"\s+", " ", tag.get_text(separator=sep, strip=True))


# ── E-utilities ───────────────────────────────────────────────
def esearch_medgen(disease_name: str, max_results: int, email: str) -> list[str]:
    r = _get(
        f"{EUTILS_BASE}/esearch.fcgi",
        params={
            "db": "medgen",
            "term": f'"{disease_name}"',
            "retmax": max_results,
            "sort": "relevance",
            "retmode": "json",
            "email": email,
        },
    )
    id_list: list[str] = r.json().get("esearchresult", {}).get("idlist", [])
    print(f"[ESearch] '{disease_name}' → {len(id_list)} UIDs: {id_list}")
    return id_list

def elink_medgen_pubmed(uid: str, email: str, max_papers: int = 30) -> list[str]:
    r = _get(
        f"{EUTILS_BASE}/elink.fcgi",
        params={
            "dbfrom": "medgen", "db": "pubmed",
            "id": uid, "retmode": "json", "email": email,
        },
    )
    pmids: list[str] = []
    for ls in r.json().get("linksets", []):
        for lsdb in ls.get("linksetdbs", []):
            if lsdb.get("dbto") == "pubmed":
                pmids.extend(str(l["id"]) for l in lsdb.get("links", []))
    return pmids[:max_papers]

def esummary_pubmed(pmids: list[str], email: str) -> list[dict]:
    if not pmids:
        return []
    r = _get(
        f"{EUTILS_BASE}/esummary.fcgi",
        params={"db": "pubmed", "id": ",".join(pmids),
                "retmode": "json", "email": email},
    )
    result_map = r.json().get("result", {})
    articles = []
    for pmid in pmids:
        e = result_map.get(pmid, {})
        if not e or pmid == "uids":
            continue
        pub_date = e.get("pubdate", "")
        articles.append({
            "pmid": pmid,
            "title": e.get("title", ""),
            "authors": [a.get("name", "") for a in e.get("authors", [])[:6]],
            "journal": e.get("source", ""),
            "year": pub_date[:4] if pub_date else "",
            "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
        })
    return articles


# ── HTML 파서: 기본 정보 ───────────────────────────────────────

def parse_title(soup: BeautifulSoup) -> str:
    """div.MedGenTitleText에서 질병명 추출 (highlight span 포함 텍스트 정제)"""
    div = soup.find("div", class_="MedGenTitleText")
    if not div:
        h1 = soup.find("h1", class_="medgenTitle")
        return _txt(h1) if h1 else ""
    return re.sub(r"\s+", " ", div.get_text(separator=" ", strip=True))


def parse_concept_info(soup: BeautifulSoup) -> dict:
    """dl.rprtid에서 MedGen UID, Concept ID, 타입 추출"""
    result = {"uid": "", "concept_id": "", "type": ""}
    dl = soup.find("dl", class_="rprtid")
    if not dl:
        return result
    dts = dl.find_all("dt")
    dds = dl.find_all("dd")
    for dt, dd in zip(dts, dds):
        key = _txt(dt).lower()
        val = _txt(dd)
        if "medgen uid" in key:
            result["uid"] = val.strip()
        elif "concept id" in key:
            result["concept_id"] = val.strip()
        elif val and not result["type"]:
            # 세 번째 dd는 타입 (Disease or Syndrome 등)
            if result["concept_id"]:
                result["type"] = val.strip()
    return result


def parse_definition(soup: BeautifulSoup) -> str:
    """
    실제 DOM: div.portlet#ID_100 > div.portlet_content.ln
    → h1#Definition 헤딩의 부모 portlet 안의 portlet_content
    """
    # 방법 1: id="ID_100" portlet 직접 탐색
    portlet = soup.find("div", id="ID_100")
    if portlet:
        content = portlet.find("div", class_="portlet_content")
        if content:
            return re.sub(r"\s+", " ", content.get_text(separator=" ", strip=True))

    # 방법 2: h1#Definition → 부모 portlet → portlet_content
    h1 = soup.find("h1", id="Definition")
    if h1:
        portlet_head = h1.find_parent("div", class_="portlet_head")
        if portlet_head:
            portlet = portlet_head.find_parent("div", class_="portlet")
            if portlet:
                content = portlet.find("div", class_="portlet_content")
                if content:
                    return re.sub(r"\s+", " ", content.get_text(separator=" ", strip=True))

    # 방법 3: meta description fallback
    meta = soup.find("meta", attrs={"name": "description"})
    if meta and meta.get("content"):
        return meta["content"].strip()

    return ""


def parse_synonyms(soup: BeautifulSoup) -> list[str]:
    """medgenTable의 'Synonyms:' 행에서 추출"""
    table = soup.find("table", class_="medgenTable")
    if not table:
        return []
    for tr in table.find_all("tr"):
        tds = tr.find_all("td")
        if len(tds) >= 2 and "synonym" in tds[0].get_text().lower():
            raw = tds[1].get_text(separator=";", strip=True)
            return [s.strip() for s in raw.split(";") if s.strip()]
    return []


def parse_cross_references(soup: BeautifulSoup) -> dict[str, str]:
    """
    medgenTable에서 MONDO, OMIM, Orphanet, SNOMED-CT 추출.
    각 행의 td[0] 텍스트가 DB명, td[1]의 a 태그 텍스트/href가 ID.
    """
    xrefs: dict[str, str] = {}
    table = soup.find("table", class_="medgenTable")
    if not table:
        return xrefs

    for tr in table.find_all("tr"):
        tds = tr.find_all("td")
        if len(tds) < 2:
            continue
        label = _txt(tds[0]).strip().rstrip(":").strip()
        a = tds[1].find("a")
        val = _txt(a) if a else _txt(tds[1])
        val = val.strip()
        if not val:
            continue

        label_lower = label.lower()
        if "monarch" in label_lower or "mondo" in label_lower:
            # "MONDO:0008608" 형태
            xrefs["MONDO"] = val
        elif "omim" in label_lower:
            xrefs["OMIM"] = val
        elif "orphanet" in label_lower:
            # "ORPHA870" → "870"
            xrefs["Orphanet"] = re.sub(r"^ORPHA", "", val)
        elif "snomed" in label_lower:
            # "Trisomy 21 (737542000)" → "737542000"
            m = re.search(r"\((\d+)\)", val)
            xrefs["SNOMED-CT"] = m.group(1) if m else val

    return xrefs


# ── HTML 파서: 논문 ────────────────────────────────────────────

def _parse_see_all(container: Tag) -> dict:
    """
    'See all (N)' 링크 파싱.
    실제 DOM: <div><a data-ga-label="See all (731)" href="...">See all (731)</a></div>
    → {"count": 731, "url": "https://pubmed..."}
    반환이 없으면 {} 반환.
    """
    for a in container.find_all("a", attrs={"data-ga-label": re.compile(r"See all", re.I)}):
        label = a.get("data-ga-label", "")
        m = re.search(r"See all \((\d+)\)", label)
        if m:
            href = a.get("href", "")
            if href.startswith("/"):
                href = "https://www.ncbi.nlm.nih.gov" + href
            return {"count": int(m.group(1)), "url": href}
    return {}


def _parse_paper_block(nl_div: Tag, detail_div: Tag | None) -> dict | None:
    """
    논문 한 편 파싱.
    nl_div  : <div class="nl"><a>제목</a></div>
    detail_div: <div class="portlet_content ln">저자<br/><span>저널</span> ... PMID: <a>숫자</a></div>
    """
    # 제목
    a_title = nl_div.find("a")
    if not a_title:
        return None
    title = re.sub(r"\s+", " ", a_title.get_text(strip=True))

    # "See all" 링크는 제목이 아님 → 건너뜀
    ga_label = a_title.get("data-ga-label", "")
    if re.match(r"See all", ga_label, re.I):
        return None

    # PMID
    pmid = ""
    href = a_title.get("href", "")
    m = re.search(r"/pubmed/(\d+)", href)
    if m:
        pmid = m.group(1)

    if not detail_div:
        return {
            "pmid": pmid, "title": title,
            "authors": [], "journal": "", "citation": "",
            "pmc_free": False,
            "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else "",
        }

    # 저자
    author_span = detail_div.find("span", class_="medgenPMauthor")
    authors = []
    if author_span:
        raw_authors = _txt(author_span)
        authors = [a.strip() for a in raw_authors.split(",") if a.strip()]

    # 저널 + 인용정보
    journal_span = detail_div.find("span", class_="medgenPMjournal")
    journal = _txt(journal_span) if journal_span else ""

    # 인용 (저널 이후 텍스트: 연도/권호/doi 등)
    citation_parts = []
    if journal_span:
        for sib in journal_span.next_siblings:
            if isinstance(sib, Tag):
                cls = sib.get("class", [])
                if "medgenPMauthor" in cls or "medgenPMjournal" in cls:
                    break
                sib_txt = sib.get_text(strip=True)
                if sib_txt and sib_txt not in ("PMID:", "PMID: "):
                    citation_parts.append(sib_txt)
            else:
                raw = str(sib).strip()
                if raw:
                    citation_parts.append(raw)
    citation = re.sub(r"\s+", " ", " ".join(citation_parts)).strip()

    # PMID (detail_div 내 재확인)
    if not pmid:
        pmid_a = detail_div.find("a", attrs={"data-ga-action": "PMID"})
        if pmid_a:
            pmid = _txt(pmid_a).strip()

    # PMC Free 여부
    pmc_free = bool(detail_div.find("a", class_="PubMedFree"))

    return {
        "pmid": pmid,
        "title": title,
        "authors": authors,
        "journal": journal,
        "citation": citation,
        "pmc_free": pmc_free,
        "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else "",
    }


def _parse_portlet_papers(portlet_content: Tag) -> dict:
    """
    portlet_content 안에서 div.nl → 제목, 바로 다음 div.portlet_content.ln → 상세
    h3.subhead 기준으로 카테고리 분리.

    반환: {
        "category_key": {
            "papers": [...],
            "see_all": {"count": N, "url": "..."}
        },
        ...
    }
    단일 카테고리이면 key="papers" 없이 {"papers":[], "see_all":{}} 반환.
    """
    result: dict = {}
    current_cat = "_default"

    children = list(portlet_content.children)
    i = 0
    while i < len(children):
        node = children[i]

        if not isinstance(node, Tag):
            i += 1
            continue

        # h3.subhead → 카테고리 전환
        if node.name == "h3" and "subhead" in node.get("class", []):
            h3_txt = _txt(node).lower()
            current_cat = "_default"
            for key_pat, cat_key in SUBHEAD_MAP.items():
                if key_pat in h3_txt:
                    current_cat = cat_key
                    break
            # 매핑 안 된 subhead (PubMed, Curated 등)는 그대로 소문자로
            if current_cat == "_default":
                current_cat = re.sub(r"\s+", "_", h3_txt.strip())
            if current_cat not in result:
                result[current_cat] = {"papers": [], "see_all": {}}
            i += 1
            continue

        # div.nl → 논문 제목 블록
        if node.name == "div" and "nl" in node.get("class", []):
            # See all 링크인지 확인
            see_all = _parse_see_all(node)
            if see_all:
                if current_cat not in result:
                    result[current_cat] = {"papers": [], "see_all": {}}
                result[current_cat]["see_all"] = see_all
                i += 1
                continue

            # 다음 형제에서 detail div 찾기
            detail_div = None
            j = i + 1
            while j < len(children):
                next_node = children[j]
                if not isinstance(next_node, Tag):
                    j += 1
                    continue
                classes = next_node.get("class", [])
                if "portlet_content" in classes and "ln" in classes:
                    detail_div = next_node
                    i = j  # detail_div까지 소비
                break

            paper = _parse_paper_block(node, detail_div)
            if paper:
                if current_cat not in result:
                    result[current_cat] = {"papers": [], "see_all": {}}
                result[current_cat]["papers"].append(paper)

        # 독립 div (See all 감싸는 경우)
        elif node.name == "div" and node.find(
            "a", attrs={"data-ga-label": re.compile(r"See all", re.I)}
        ):
            see_all = _parse_see_all(node)
            if see_all:
                if current_cat not in result:
                    result[current_cat] = {"papers": [], "see_all": {}}
                result[current_cat]["see_all"] = see_all

        i += 1

    return result


def parse_professional_guidelines(soup: BeautifulSoup) -> dict:
    """
    div#ID_105 portlet → Professional guidelines 파싱.
    PubMed 논문 목록 + Curated 가이드라인(외부 URL) 포함.
    """
    portlet = soup.find("div", id="ID_105")
    if not portlet:
        return {"papers": [], "curated": [], "see_all": {}}

    content = portlet.find("div", class_="portlet_content")
    if not content:
        return {"papers": [], "curated": [], "see_all": {}}

    parsed = _parse_portlet_papers(content)

    # PubMed 논문 (subhead "PubMed" 하위)
    pubmed_section = parsed.get("pubmed", parsed.get("_default", {"papers": [], "see_all": {}}))
    papers = pubmed_section.get("papers", [])
    see_all = pubmed_section.get("see_all", {})

    # See all 이 다른 섹션에 있을 수 있음
    for v in parsed.values():
        if v.get("see_all"):
            see_all = v["see_all"]
            break

    # Curated 가이드라인 (h3.nl.vspace > a 형태)
    curated = []
    for h3 in content.find_all("h3", class_="nl"):
        a = h3.find("a")
        if a and a.get("href"):
            curated.append({
                "title": re.sub(r"\s+", " ", a.get_text(strip=True)),
                "url": a["href"],
            })

    return {"papers": papers, "curated": curated, "see_all": see_all}


def parse_recent_clinical_studies(soup: BeautifulSoup) -> dict:
    """
    div#ID_103 → Recent clinical studies
    내부에 h3.subhead 기준으로 Etiology/Diagnosis/Therapy/Prognosis/Clinical prediction guides 분리.
    """
    portlet = soup.find("div", id="ID_103")
    if not portlet:
        return {}
    content = portlet.find("div", class_="portlet_content")
    if not content:
        return {}
    return _parse_portlet_papers(content)


def parse_systematic_reviews(soup: BeautifulSoup) -> dict:
    """div#ID_104 → Recent systematic reviews"""
    portlet = soup.find("div", id="ID_104")
    if not portlet:
        return {"papers": [], "see_all": {}}
    content = portlet.find("div", class_="portlet_content")
    if not content:
        return {"papers": [], "see_all": {}}
    parsed = _parse_portlet_papers(content)
    # 단일 카테고리
    all_papers = []
    see_all = {}
    for v in parsed.values():
        all_papers.extend(v.get("papers", []))
        if v.get("see_all"):
            see_all = v["see_all"]
    return {"papers": all_papers, "see_all": see_all}


# ── 메인 레코드 수집 ───────────────────────────────────────────

def fetch_medgen_record(uid: str, email: str) -> dict[str, Any]:
    url = f"{MEDGEN_BASE}/{uid}"
    print(f"  [Fetch] UID={uid}")

    try:
        resp = _get(url, delay=1.5)
    except requests.HTTPError as e:
        print(f"  [Error] HTTP {e.response.status_code}")
        return {"medgen_uid": uid, "error": str(e), "url": url}
    except Exception as e:
        print(f"  [Error] {e}")
        return {"medgen_uid": uid, "error": str(e), "url": url}

    soup = BeautifulSoup(resp.text, "lxml")

    # 기본 정보
    title       = parse_title(soup)
    concept     = parse_concept_info(soup)
    definition  = parse_definition(soup)
    synonyms    = parse_synonyms(soup)
    xrefs       = parse_cross_references(soup)

    # 논문 섹션
    guidelines   = parse_professional_guidelines(soup)
    clinical     = parse_recent_clinical_studies(soup)
    sys_reviews  = parse_systematic_reviews(soup)

    # clinical 내 카테고리 정리 및 papers 없는 경우 ELink fallback
    total_papers = (
        len(guidelines.get("papers", []))
        + sum(len(v.get("papers", [])) for v in clinical.values())
        + len(sys_reviews.get("papers", []))
    )
    elink_fallback = []
    if total_papers == 0:
        print(f"  [Info] 논문 파싱 결과 없음 → ELink fallback")
        pmids = elink_medgen_pubmed(uid, email)
        if pmids:
            time.sleep(0.35)
            elink_fallback = esummary_pubmed(pmids, email)

    record = {
        "medgen_uid":        uid,
        "disease_name":      title,
        "concept_id":        concept.get("concept_id", ""),
        "disease_type":      concept.get("type", ""),
        "definition":        definition,
        "synonyms":          synonyms,
        "cross_references":  xrefs,
        "literature": {
            "professional_guidelines": guidelines,
            "recent_clinical_studies": clinical,
            "recent_systematic_reviews": sys_reviews,
            **({"elink_related": elink_fallback} if elink_fallback else {}),
        },
        "url": url,
    }
    return record


def fetch_medgen_details(
    disease_name: str,
    max_results: int,
    email: str,
    explicit_uids: list[str] | None = None,
) -> list[dict[str, Any]]:
    uid_list = explicit_uids or esearch_medgen(disease_name, max_results, email)
    if not uid_list:
        print("[Warning] UID 없음.")
        return []

    results = []
    for uid in uid_list:
        record = fetch_medgen_record(uid, email)
        results.append(record)
        time.sleep(1.0)
    return results


# ── CLI ───────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="NCBI MedGen 질병 정보 + 카테고리별 논문 파서",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
예시:
  python medgen_parser.py --disease "Down syndrome" --max-results 3 --output out.json
  python medgen_parser.py --uids 4385 --output out.json
  python medgen_parser.py --uids 4385 836 --output multi.json
        """,
    )
    parser.add_argument("--disease",     type=str, default="")
    parser.add_argument("--uids",        type=str, nargs="+")
    parser.add_argument("--max-results", type=int, default=1)
    parser.add_argument("--output",      type=str, default="medgen_results.json")
    parser.add_argument("--email",       type=str, default="your.email@example.com")
    parser.add_argument(
        "--log",
        type=str,
        default="",
        help="로그 파일 경로. 미지정 시 output 파일명을 기준으로 생성",
    )
    args = parser.parse_args()

    if not args.disease and not args.uids:
        parser.error("--disease 또는 --uids 중 하나는 필수입니다.")
    
    log_path = args.log or get_default_log_path(args.output)
    tee_logger = setup_tee_logging(log_path)

    try:
        print("=" * 80)
        print(f"실행 시작: {datetime.now():%Y-%m-%d %H:%M:%S}")
        print(f"결과 파일: {args.output}")
        print(f"로그 파일: {log_path}")
        print("=" * 80)

        print(
            f"▶ MedGen 수집 시작 "
            f"| disease='{args.disease}' "
            f"uids={args.uids}"
        )

        results = fetch_medgen_details(
            disease_name=args.disease,
            max_results=args.max_results,
            email=args.email,
            explicit_uids=args.uids,
        )

        print(f"\n▶ 수집 완료: {len(results)}건")

        with open(args.output, "w", encoding="utf-8") as f:
            json.dump(
                results,
                f,
                ensure_ascii=False,
                indent=2,
            )

        print(f"▶ 저장: {args.output}")

        for rec in results:
            lit = rec.get("literature", {})
            gl = lit.get("professional_guidelines", {})
            cs = lit.get("recent_clinical_studies", {})
            sr = lit.get("recent_systematic_reviews", {})

            n_gl = len(gl.get("papers", []))
            n_cs = sum(
                len(category.get("papers", []))
                for category in cs.values()
            )
            n_sr = len(sr.get("papers", []))
            cats = list(cs.keys())

            print(
                f"  UID={rec.get('medgen_uid')} "
                f"| {rec.get('disease_name')} "
                f"| def={'O' if rec.get('definition') else 'X'} "
                f"| xrefs={list(rec.get('cross_references', {}).keys())} "
                f"| guidelines={n_gl} "
                f"| clinical_cats={cats}({n_cs}편) "
                f"| sysrev={n_sr}"
            )

        print("=" * 80)
        print(f"실행 완료: {datetime.now():%Y-%m-%d %H:%M:%S}")
        print("=" * 80)

    except Exception:
        print("\n[ERROR] 실행 중 오류가 발생했습니다.")
        traceback.print_exc()
        raise

    finally:
        tee_logger.close()

if __name__ == "__main__":
    main()