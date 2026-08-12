"""
parsers/pubmed_parser.py
─────────────────────────
두 가지 독립적인 기능을 제공:

  [fetch]    PMID 목록 → articles dict
  [discover] 질환명    → PubMed 검색 → PMID 발굴

INPUT  (fetch)   : pmids: list[str]
INPUT  (discover): disease_name: str, synonyms: list[str]
OUTPUT           : articles: {pmid: article_dict}

article_dict:
  pmid, doi, pmc_id, pmc_verify,
  title, authors, journal, year, publication_types,
  abstract, abstract_sections,
  mesh_terms, keywords,
  full_text_available, pdf_links, pubmed_url

단독 실행:
  # PMID로 fetch
  python -m nipt_pipeline.parsers.pubmed_parser \\
      --pmids 38745141 33053283 --output articles.json

  # 질환명으로 discover
  python -m nipt_pipeline.parsers.pubmed_parser \\
      --discover "Down syndrome" --output articles.json

  # 둘 다
  python -m nipt_pipeline.parsers.pubmed_parser \\
      --pmids 38745141 --discover "Down syndrome" --output articles.json
"""

import re
import sys
import json
import time
import argparse
import pathlib
import xml.etree.ElementTree as ET
from typing import Any, Iterable

from nipt_pipeline.utils import EUTILS, PUBMED, PMC_BASE, http_get

# ── PMC 역검증 캐시 ───────────────────────────────────────────
_pmc_cache: dict[str, str] = {}


def verify_pmc_id(pmc_id: str, expected_pmid: str, email: str) -> dict:
    """PMC ID 역조회 → PMID 일치 여부 검증."""
    pmc_num = pmc_id.replace("PMC", "").strip()
    if pmc_num not in _pmc_cache:
        try:
            r = http_get(f"{EUTILS}/esearch.fcgi",
                         params={"db": "pubmed",
                                 "term": f"PMC{pmc_num}[pmc]",
                                 "retmode": "json", "email": email})
            ids = r.json().get("esearchresult", {}).get("idlist", [])
            _pmc_cache[pmc_num] = ids[0] if ids else ""
        except Exception as e:
            return {"pmc_id": pmc_id, "verified": False,
                    "matched_pmid": None, "pmc_type": "error",
                    "note": f"역조회 실패: {e}"}

    matched = _pmc_cache[pmc_num]
    if not matched:
        return {"pmc_id": pmc_id, "verified": False,
                "matched_pmid": None, "pmc_type": "not_found",
                "note": f"PMC{pmc_num}에 해당 PMID 없음"}
    if matched == str(expected_pmid):
        return {"pmc_id": pmc_id, "verified": True,
                "matched_pmid": matched, "pmc_type": "publisher",
                "note": "검증 통과"}
    return {"pmc_id": pmc_id, "verified": False,
            "matched_pmid": matched,
            "pmc_type": "author_manuscript_or_mismatch",
            "note": (f"역조회={matched} ≠ 현재={expected_pmid}. "
                     "AM 버전 또는 메타데이터 오류 가능성.")}


def _build_pdf_links(pmid: str, doi: str, pmc_id: str,
                     pmc_verify: dict | None) -> tuple[list[dict], bool]:
    links: list[dict] = []
    full_text = False

    if pmc_id:
        verified = pmc_verify.get("verified") if pmc_verify else None
        warn     = (pmc_verify.get("note")
                    if pmc_verify and not pmc_verify.get("verified") else None)
        links += [
            {"source": "PMC", "type": "pdf",
             "url": f"{PMC_BASE}/{pmc_id}/pdf/",
             "verified": verified,
             **({"warning": warn} if warn else {})},
            {"source": "PMC", "type": "html",
             "url": f"{PMC_BASE}/{pmc_id}/",
             "verified": verified,
             **({"warning": warn} if warn else {})},
        ]
        if verified is not False:
            full_text = True

    if doi:
        links.append({"source": "DOI", "type": "publisher",
                      "url": f"https://doi.org/{doi}"})
    links.append({"source": "PubMed", "type": "abstract",
                  "url": f"{PUBMED}/{pmid}/"})
    return links, full_text


def _parse_abstract(art: ET.Element) -> tuple[str, dict]:
    abs_el = art.find(".//Abstract")
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


def _parse_authors(art: ET.Element) -> list[str]:
    out = []
    for au in art.findall(".//AuthorList/Author"):
        last  = (au.findtext("LastName") or "").strip()
        fore  = (au.findtext("ForeName") or au.findtext("Initials") or "").strip()
        cname = (au.findtext("CollectiveName") or "").strip()
        if cname:
            out.append(cname)
        elif last:
            out.append(f"{last} {fore}".strip())
    return out


def _parse_article_element(pa: ET.Element,
                            verify_pmc: bool = False,
                            email: str = "") -> dict | None:
    mc  = pa.find("MedlineCitation")
    art = mc.find("Article") if mc is not None else None
    if mc is None or art is None:
        return None

    pmid_el = mc.find("PMID")
    pmid = pmid_el.text.strip() if pmid_el is not None else ""

    title_el = art.find("ArticleTitle")
    title = re.sub(r"\s+", " ",
                   "".join(title_el.itertext()).strip()
                   if title_el is not None else "")

    abstract, abs_sections = _parse_abstract(art)
    authors = _parse_authors(art)
    journal = (art.findtext(".//Journal/Title") or
               art.findtext(".//Journal/ISOAbbreviation") or "")

    pd_el = art.find(".//PubDate")
    year  = ""
    if pd_el is not None:
        year = pd_el.findtext("Year") or ""
        if not year:
            med = pd_el.findtext("MedlineDate") or ""
            m   = re.match(r"(\d{4})", med)
            year = m.group(1) if m else ""

    pub_types = [pt.text.strip()
                 for pt in art.findall(".//PublicationTypeList/PublicationType")
                 if pt.text]

    doi, pmc_id = "", ""
    pd_xml = pa.find("PubmedData")
    if pd_xml is not None:
        for aid in pd_xml.findall(".//ArticleIdList/ArticleId"):
            t, v = aid.get("IdType", ""), (aid.text or "").strip()
            if t == "doi" and not doi:
                doi = v
            elif t == "pmc" and not pmc_id:
                pmc_id = v
    for loc in art.findall("ELocationID"):
        if loc.get("EIdType") == "doi":
            doi = (loc.text or "").strip()
            break

    mesh_terms = [mh.findtext("DescriptorName") or ""
                  for mh in mc.findall(".//MeshHeadingList/MeshHeading")]
    keywords   = [kw.text.strip()
                  for kw in mc.findall(".//KeywordList/Keyword") if kw.text]

    pmc_verify = None
    if pmc_id and verify_pmc and email:
        pmc_verify = verify_pmc_id(pmc_id, pmid, email)

    pdf_links, full_text = _build_pdf_links(pmid, doi, pmc_id, pmc_verify)

    return {
        "pmid":               pmid,
        "doi":                doi,
        "pmc_id":             pmc_id,
        "pmc_verify":         pmc_verify,
        "title":              title,
        "authors":            authors,
        "journal":            journal,
        "year":               year,
        "publication_types":  pub_types,
        "abstract":           abstract,
        "abstract_sections":  abs_sections or None,
        "mesh_terms":         [m for m in mesh_terms if m],
        "keywords":           keywords,
        "full_text_available":full_text,
        "pdf_links":          pdf_links,
        "pubmed_url":         f"{PUBMED}/{pmid}/",
    }


# ════════════════════════════════════════════════════════════════
# FETCH: PMID 목록 → articles
# ════════════════════════════════════════════════════════════════

def fetch_articles(pmids: list[str], email: str,
                   verify_pmc: bool = True,
                   batch_size: int = 100) -> dict[str, dict]:
    """
    INPUT : pmids (중복 허용, 자동 제거)
    OUTPUT: {pmid: article_dict}
    """
    # 중복 제거 (순서 유지)
    seen: set[str] = set()
    unique = [p for p in pmids if p and p not in seen and not seen.add(p)]
    if not unique:
        return {}

    lookup: dict[str, dict] = {}
    total = len(unique)
    for i in range(0, total, batch_size):
        batch = unique[i: i + batch_size]
        print(f"  [PubMed fetch] {i+1}~{min(i+batch_size, total)}/{total}")
        try:
            r = http_get(f"{EUTILS}/efetch.fcgi",
                         params={"db": "pubmed", "id": ",".join(batch),
                                 "rettype": "abstract", "retmode": "xml",
                                 "email": email})
            root = ET.fromstring(r.text.encode("utf-8"))
        except Exception as e:
            print(f"  [PubMed fetch] 오류: {e}")
            for pmid in batch:
                lookup[pmid] = {"pmid": pmid, "error": str(e)}
            continue

        for pa in root.findall(".//PubmedArticle"):
            article = _parse_article_element(pa, verify_pmc=verify_pmc,
                                              email=email)
            if article and article.get("pmid"):
                lookup[article["pmid"]] = article

    print(f"  [PubMed fetch] 완료: {len(lookup)}/{total}")
    return lookup


# ════════════════════════════════════════════════════════════════
# DISCOVER: 질환명 → PMID 발굴
# ════════════════════════════════════════════════════════════════

def _dedup(vals: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    return [v for v in vals if v and v not in seen and not seen.add(v)]


def _build_disease_query(disease_name: str,
                          synonyms: list[str] | None,
                          field: str = "Title/Abstract",
                          max_synonyms: int = 6) -> str:
    names = _dedup([disease_name, *(synonyms or [])])[:max_synonyms + 1]
    cleaned = [re.sub(r'"', "", n).strip() for n in names if len(n.strip()) >= 4]
    if not cleaned:
        return ""
    clauses = [f'"{n}"[{field}]' for n in cleaned]
    return "(" + " OR ".join(clauses) + ")"


def discover_pmids(disease_name: str,
                   synonyms: list[str] | None,
                   email: str,
                   limits: dict[str, int] | None = None) -> dict[str, dict]:
    """
    질환명 기반 PubMed bucket 검색.

    INPUT : disease_name, synonyms
    OUTPUT: {bucket: {query, pmids}}
      bucket 종류: review, genetics, genotype_phenotype, cohort, general
    """
    default_limits = {
        "review":              5,
        "genetics":           10,
        "genotype_phenotype": 10,
        "cohort":              5,
        "general":            10,
    }
    if limits:
        default_limits.update(limits)

    dq  = _build_disease_query(disease_name, synonyms)
    if not dq:
        return {}
    flt = "English[lang] NOT comment[pt] NOT letter[pt]"

    buckets: dict[str, str] = {
        "review": (
            f"{dq} AND (review[pt] OR systematic review[pt]) AND {flt}"
        ),
        "genetics": (
            f"{dq} AND ("
            "gene[Title/Abstract] OR genes[Title/Abstract] OR "
            "genetic[Title/Abstract] OR haploinsufficiency[Title/Abstract] OR "
            "deletion[Title/Abstract] OR duplication[Title/Abstract] OR "
            '"copy number"[Title/Abstract]'
            f") AND {flt}"
        ),
        "genotype_phenotype": (
            f"{dq} AND ("
            '"genotype-phenotype"[Title/Abstract] OR '
            '"critical region"[Title/Abstract] OR '
            "mapping[Title/Abstract] OR phenotype[Title/Abstract]"
            f") AND {flt}"
        ),
        "cohort": (
            f"{dq} AND ("
            "cohort[Title/Abstract] OR patients[Title/Abstract] OR "
            '"case series"[Title/Abstract]'
            f") AND {flt}"
        ),
        "general": f"{dq} AND {flt}",
    }

    result: dict[str, dict] = {}
    for bucket, query in buckets.items():
        limit = default_limits.get(bucket, 0)
        if not limit:
            continue
        try:
            r = http_get(f"{EUTILS}/esearch.fcgi",
                         params={"db": "pubmed", "term": query,
                                 "retmax": limit, "sort": "relevance",
                                 "retmode": "json", "email": email})
            pmids = r.json().get("esearchresult", {}).get("idlist", [])
            result[bucket] = {"query": query, "pmids": pmids}
            print(f"  [PubMed discover] {bucket}: {len(pmids)}편")
        except Exception as e:
            result[bucket] = {"query": query, "pmids": [], "error": str(e)}
        time.sleep(0.35)

    return result


def collect_discovery_pmids(discovery: dict[str, dict]) -> list[str]:
    """discover_pmids 결과에서 PMID 병합/중복 제거."""
    all_pmids: list[str] = []
    for bd in discovery.values():
        all_pmids.extend(bd.get("pmids", []))
    return _dedup(all_pmids)


def collect_literature_pmids(literature: dict) -> list[str]:
    """medgen_parser literature에서 PMID 수집."""
    all_pmids: list[str] = []
    for cat_data in literature.values():
        if isinstance(cat_data, dict):
            for paper in cat_data.get("papers", []):
                p = str(paper.get("pmid", "")).strip()
                if p:
                    all_pmids.append(p)
    return _dedup(all_pmids)


def discover_and_fetch(disease_name: str,
                       synonyms: list[str] | None,
                       email: str,
                       medgen_literature: dict | None = None,
                       limits: dict[str, int] | None = None,
                       verify_pmc: bool = False) -> dict[str, Any]:
    """
    convenience: discover → merge MedGen PMIDs → fetch
    반환: {discovery, pmids, articles}
    """
    discovery = discover_pmids(disease_name, synonyms, email, limits)
    disc_pmids = collect_discovery_pmids(discovery)
    lit_pmids  = collect_literature_pmids(medgen_literature or {})
    pmids      = _dedup([*lit_pmids, *disc_pmids])

    print(f"  [PubMed] 총 {len(pmids)}개 PMID (MedGen={len(lit_pmids)}, discover={len(disc_pmids)})")
    articles = fetch_articles(pmids, email, verify_pmc=verify_pmc)
    return {"discovery": discovery, "pmids": pmids, "articles": articles}


# ── 단독 실행 ─────────────────────────────────────────────────
def _main():
    ap = argparse.ArgumentParser(
        description="PubMed fetch(PMID) / discover(질환명) → articles JSON",
        epilog=(
            "# PMID fetch\n"
            "python -m nipt_pipeline.parsers.pubmed_parser "
            "--pmids 38745141 33053283 --output articles.json\n\n"
            "# 질환명 discover\n"
            "python -m nipt_pipeline.parsers.pubmed_parser "
            "--discover 'Down syndrome' --output articles.json"
        ))
    ap.add_argument("--pmids",    nargs="*", default=[])
    ap.add_argument("--discover", type=str,  default=None,
                    help="질환명 (PubMed discovery 모드)")
    ap.add_argument("--synonyms", nargs="*", default=[])
    ap.add_argument("--email",    default="your@email.com")
    ap.add_argument("--no-verify-pmc", action="store_true")
    ap.add_argument("--output",   default="-")
    args = ap.parse_args()

    from nipt_pipeline.utils import init_session
    init_session(args.email)

    pmids: list[str] = [p.strip() for p in args.pmids if p.strip().isdigit()]
    discovery: dict = {}

    if args.discover:
        result = discover_and_fetch(
            disease_name=args.discover,
            synonyms=args.synonyms or None,
            email=args.email,
            verify_pmc=not args.no_verify_pmc,
        )
        # 직접 지정 PMID 추가
        if pmids:
            extra = fetch_articles(pmids, args.email,
                                   verify_pmc=not args.no_verify_pmc)
            result["articles"].update(extra)
        out_data = result
    else:
        if not pmids:
            ap.error("--pmids 또는 --discover 중 하나 필수")
        articles = fetch_articles(pmids, args.email,
                                  verify_pmc=not args.no_verify_pmc)
        out_data = {"pmids": pmids, "articles": articles}

    out = json.dumps(out_data, ensure_ascii=False, indent=2)
    if args.output == "-":
        print(out)
    else:
        pathlib.Path(args.output).write_text(out, encoding="utf-8")
        print(f"저장: {args.output}", file=sys.stderr)


import json    # noqa
import pathlib # noqa

if __name__ == "__main__":
    _main()