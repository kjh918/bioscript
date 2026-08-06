"""
parsers/pubmed_parser.py
─────────────────────────
PubMed EFetch 파싱 + PMC ID 역검증.
"""

import re
import time
import xml.etree.ElementTree as ET
from typing import Any

from resource_parse_pipeline.utils import EUTILS, PUBMED, PMC_BASE, http_get

_pmc_verify_cache: dict[str, str] = {}


def verify_pmc_id(pmc_id: str, expected_pmid: str, email: str) -> dict:
    """PMC ID 역조회 → PMID 일치 여부 검증."""
    pmc_num = pmc_id.replace("PMC", "").strip()
    if pmc_num in _pmc_verify_cache:
        matched = _pmc_verify_cache[pmc_num]
    else:
        try:
            r = http_get(f"{EUTILS}/esearch.fcgi",
                         params={"db": "pubmed",
                                 "term": f"PMC{pmc_num}[pmc]",
                                 "retmode": "json", "email": email})
            ids = r.json().get("esearchresult", {}).get("idlist", [])
            matched = ids[0] if ids else ""
            _pmc_verify_cache[pmc_num] = matched
        except Exception as e:
            return {"pmc_id": pmc_id, "verified": False,
                    "matched_pmid": "", "pmc_type": "error",
                    "note": f"역조회 실패: {e}"}

    if not matched:
        return {"pmc_id": pmc_id, "verified": False,
                "matched_pmid": "", "pmc_type": "not_found",
                "note": f"PMC{pmc_num}에 해당하는 PMID 없음"}
    if matched == str(expected_pmid):
        return {"pmc_id": pmc_id, "verified": True,
                "matched_pmid": matched, "pmc_type": "publisher",
                "note": "PMC ID 검증 통과"}
    return {"pmc_id": pmc_id, "verified": False,
            "matched_pmid": matched,
            "pmc_type": "author_manuscript_or_mismatch",
            "note": (f"PMC{pmc_num} 역조회={matched} ≠ {expected_pmid}. "
                     "Author Manuscript 또는 메타데이터 오류 가능성.")}


def _build_pdf_links(pmid: str, doi: str, pmc_id: str,
                     pmc_verify: dict | None) -> tuple[list[dict], bool]:
    links: list[dict] = []
    full_text = False
    if pmc_id:
        ok = pmc_verify["verified"] if pmc_verify else None
        warn = pmc_verify.get("note") if pmc_verify and not (pmc_verify.get("verified")) else None
        links += [
            {"source": "PMC", "type": "pdf",
             "url": f"{PMC_BASE}/{pmc_id}/pdf/",
             "verified": ok, **({"warning": warn} if warn else {})},
            {"source": "PMC", "type": "html",
             "url": f"{PMC_BASE}/{pmc_id}/",
             "verified": ok, **({"warning": warn} if warn else {})},
        ]
        if ok is not False:
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


def _parse_single_article(pa: ET.Element,
                           verify_pmc: bool = False,
                           email: str = "") -> dict | None:
    """PubmedArticle XML element → dict"""
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
    year = ""
    if pd_el is not None:
        year = pd_el.findtext("Year") or ""
        if not year:
            med = pd_el.findtext("MedlineDate") or ""
            m = re.match(r"(\d{4})", med)
            year = m.group(1) if m else ""

    pub_types = [
        pt.text.strip()
        for pt in art.findall(".//PublicationTypeList/PublicationType")
        if pt.text
    ]

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

    mesh_terms = [
        mh.findtext("DescriptorName") or ""
        for mh in mc.findall(".//MeshHeadingList/MeshHeading")
    ]
    keywords = [
        kw.text.strip()
        for kw in mc.findall(".//KeywordList/Keyword") if kw.text
    ]

    pmc_verify = None
    if pmc_id and verify_pmc and email:
        pmc_verify = verify_pmc_id(pmc_id, pmid, email)

    pdf_links, full_text = _build_pdf_links(pmid, doi, pmc_id, pmc_verify)

    return {
        "pmid":              pmid,
        "doi":               doi,
        "pmc_id":            pmc_id,
        "pmc_verify":        pmc_verify,
        "title":             title,
        "authors":           authors,
        "journal":           journal,
        "year":              year,
        "publication_types": pub_types,
        "abstract":          abstract,
        "abstract_sections": abs_sections,
        "mesh_terms":        [m for m in mesh_terms if m],
        "keywords":          keywords,
        "full_text_available": full_text,
        "pdf_links":         pdf_links,
        "pubmed_url":        f"{PUBMED}/{pmid}/",
    }


def fetch_pubmed_abstracts(
    pmids: list[str],
    email: str,
    verify_pmc: bool = True,
    batch_size: int = 100,
) -> dict[str, dict]:
    """PMID 목록 → {pmid: article_dict}"""
    if not pmids:
        return {}
    lookup: dict[str, dict] = {}
    for i in range(0, len(pmids), batch_size):
        batch = pmids[i: i + batch_size]
        print(f"  [PubMed] EFetch {i+1}~{min(i+batch_size, len(pmids))}/{len(pmids)}")
        try:
            r = http_get(f"{EUTILS}/efetch.fcgi",
                         params={"db": "pubmed", "id": ",".join(batch),
                                 "rettype": "abstract", "retmode": "xml",
                                 "email": email})
            root = ET.fromstring(r.text.encode("utf-8"))
        except Exception as e:
            print(f"  [PubMed] 오류: {e}")
            for pmid in batch:
                lookup[pmid] = {"pmid": pmid, "error": str(e)}
            continue

        for pa in root.findall(".//PubmedArticle"):
            article = _parse_single_article(pa, verify_pmc=verify_pmc,
                                             email=email)
            if article and article.get("pmid"):
                lookup[article["pmid"]] = article
                if verify_pmc and article.get("pmc_verify"):
                    v = article["pmc_verify"]
                    mark = "✓" if v["verified"] else "✗"
                    print(f"    PMC {mark} PMID={article['pmid']} "
                          f"PMC={article['pmc_id']} → {v['note'][:60]}")
    return lookup


def collect_pmids_from_literature(literature: dict) -> list[str]:
    """literature dict에서 PMID 전체 수집 (중복 제거, 순서 유지)"""
    seen: set[str] = set()
    pmids: list[str] = []
    for cat_data in literature.values():
        if isinstance(cat_data, dict):
            for paper in cat_data.get("papers", []):
                p = str(paper.get("pmid", "")).strip()
                if p and p not in seen:
                    seen.add(p)
                    pmids.append(p)
    return pmids