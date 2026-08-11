"""
parsers/pubmed_parser.py
─────────────────────────
PubMed ESearch/EFetch 파싱 + PMC ID 역검증 + 질환 기반 문헌 discovery.

주요 흐름:
  1) disease name + synonyms -> PubMed ESearch (bucket별 discovery)
  2) PMID merge/dedup
  3) PubMed EFetch -> title/abstract/MeSH/keyword/PMCID
  4) 확보 문헌에서 candidate gene symbol 추출
  5) candidate gene x disease 후속 검색에 재사용 가능
"""

import re
import time
import xml.etree.ElementTree as ET
from typing import Any, Iterable

from resource_parse_pipeline.utils import EUTILS, PUBMED, PMC_BASE, http_get

_pmc_verify_cache: dict[str, str] = {}


def verify_pmc_id(pmc_id: str, expected_pmid: str, email: str) -> dict:
    """PMC ID 역조회 -> PMID 일치 여부 검증."""
    pmc_num = pmc_id.replace("PMC", "").strip()
    if pmc_num in _pmc_verify_cache:
        matched = _pmc_verify_cache[pmc_num]
    else:
        try:
            r = http_get(
                f"{EUTILS}/esearch.fcgi",
                params={
                    "db": "pubmed",
                    "term": f"PMC{pmc_num}[pmc]",
                    "retmode": "json",
                    "email": email,
                },
            )
            ids = r.json().get("esearchresult", {}).get("idlist", [])
            matched = ids[0] if ids else ""
            _pmc_verify_cache[pmc_num] = matched
        except Exception as e:
            return {
                "pmc_id": pmc_id,
                "verified": False,
                "matched_pmid": "",
                "pmc_type": "error",
                "note": f"역조회 실패: {e}",
            }

    if not matched:
        return {
            "pmc_id": pmc_id,
            "verified": False,
            "matched_pmid": "",
            "pmc_type": "not_found",
            "note": f"PMC{pmc_num}에 해당하는 PMID 없음",
        }
    if matched == str(expected_pmid):
        return {
            "pmc_id": pmc_id,
            "verified": True,
            "matched_pmid": matched,
            "pmc_type": "publisher",
            "note": "PMC ID 검증 통과",
        }
    return {
        "pmc_id": pmc_id,
        "verified": False,
        "matched_pmid": matched,
        "pmc_type": "author_manuscript_or_mismatch",
        "note": (
            f"PMC{pmc_num} 역조회={matched} != {expected_pmid}. "
            "Author Manuscript 또는 메타데이터 오류 가능성."
        ),
    }


def _build_pdf_links(
    pmid: str,
    doi: str,
    pmc_id: str,
    pmc_verify: dict | None,
) -> tuple[list[dict], bool]:
    links: list[dict] = []
    full_text = False
    if pmc_id:
        ok = pmc_verify["verified"] if pmc_verify else None
        warn = (
            pmc_verify.get("note")
            if pmc_verify and not pmc_verify.get("verified")
            else None
        )
        links += [
            {
                "source": "PMC",
                "type": "pdf",
                "url": f"{PMC_BASE}/{pmc_id}/pdf/",
                "verified": ok,
                **({"warning": warn} if warn else {}),
            },
            {
                "source": "PMC",
                "type": "html",
                "url": f"{PMC_BASE}/{pmc_id}/",
                "verified": ok,
                **({"warning": warn} if warn else {}),
            },
        ]
        if ok is not False:
            full_text = True
    if doi:
        links.append(
            {
                "source": "DOI",
                "type": "publisher",
                "url": f"https://doi.org/{doi}",
            }
        )
    links.append(
        {
            "source": "PubMed",
            "type": "abstract",
            "url": f"{PUBMED}/{pmid}/",
        }
    )
    return links, full_text


def _parse_abstract(art: ET.Element) -> tuple[str, dict]:
    abs_el = art.find(".//Abstract")
    if abs_el is None:
        return "", {}
    sections: dict[str, str] = {}
    parts: list[str] = []
    for at in abs_el.findall("AbstractText"):
        label = at.get("Label") or at.get("NlmCategory") or ""
        text = re.sub(r"\s+", " ", "".join(at.itertext()).strip())
        if label:
            sections[label] = text
            parts.append(f"{label}: {text}")
        else:
            parts.append(text)
    return " ".join(parts), sections


def _parse_authors(art: ET.Element) -> list[str]:
    out = []
    for au in art.findall(".//AuthorList/Author"):
        last = (au.findtext("LastName") or "").strip()
        fore = (au.findtext("ForeName") or au.findtext("Initials") or "").strip()
        cname = (au.findtext("CollectiveName") or "").strip()
        if cname:
            out.append(cname)
        elif last:
            out.append(f"{last} {fore}".strip())
    return out


def _parse_single_article(
    pa: ET.Element,
    verify_pmc: bool = False,
    email: str = "",
) -> dict | None:
    """PubmedArticle XML element -> dict"""
    mc = pa.find("MedlineCitation")
    art = mc.find("Article") if mc is not None else None
    if mc is None or art is None:
        return None

    pmid_el = mc.find("PMID")
    pmid = pmid_el.text.strip() if pmid_el is not None else ""

    title_el = art.find("ArticleTitle")
    title = re.sub(
        r"\s+",
        " ",
        "".join(title_el.itertext()).strip() if title_el is not None else "",
    )

    abstract, abs_sections = _parse_abstract(art)
    authors = _parse_authors(art)
    journal = (
        art.findtext(".//Journal/Title")
        or art.findtext(".//Journal/ISOAbbreviation")
        or ""
    )

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
        for kw in mc.findall(".//KeywordList/Keyword")
        if kw.text
    ]

    pmc_verify = None
    if pmc_id and verify_pmc and email:
        pmc_verify = verify_pmc_id(pmc_id, pmid, email)

    pdf_links, full_text = _build_pdf_links(pmid, doi, pmc_id, pmc_verify)

    return {
        "pmid": pmid,
        "doi": doi,
        "pmc_id": pmc_id,
        "pmc_verify": pmc_verify,
        "title": title,
        "authors": authors,
        "journal": journal,
        "year": year,
        "publication_types": pub_types,
        "abstract": abstract,
        "abstract_sections": abs_sections,
        "mesh_terms": [m for m in mesh_terms if m],
        "keywords": keywords,
        "full_text_available": full_text,
        "pdf_links": pdf_links,
        "pubmed_url": f"{PUBMED}/{pmid}/",
    }


def fetch_pubmed_abstracts(
    pmids: list[str],
    email: str,
    verify_pmc: bool = True,
    batch_size: int = 100,
) -> dict[str, dict]:
    """PMID 목록 -> {pmid: article_dict}"""
    if not pmids:
        return {}
    lookup: dict[str, dict] = {}
    for i in range(0, len(pmids), batch_size):
        batch = pmids[i : i + batch_size]
        print(
            f"  [PubMed] EFetch {i+1}~"
            f"{min(i+batch_size, len(pmids))}/{len(pmids)}"
        )
        try:
            r = http_get(
                f"{EUTILS}/efetch.fcgi",
                params={
                    "db": "pubmed",
                    "id": ",".join(batch),
                    "rettype": "abstract",
                    "retmode": "xml",
                    "email": email,
                },
            )
            root = ET.fromstring(r.text.encode("utf-8"))
        except Exception as e:
            print(f"  [PubMed] 오류: {e}")
            for pmid in batch:
                lookup[pmid] = {"pmid": pmid, "error": str(e)}
            continue

        for pa in root.findall(".//PubmedArticle"):
            article = _parse_single_article(
                pa,
                verify_pmc=verify_pmc,
                email=email,
            )
            if article and article.get("pmid"):
                lookup[article["pmid"]] = article
                if verify_pmc and article.get("pmc_verify"):
                    v = article["pmc_verify"]
                    mark = "OK" if v["verified"] else "FAIL"
                    print(
                        f"    PMC {mark} PMID={article['pmid']} "
                        f"PMC={article['pmc_id']} -> {v['note'][:60]}"
                    )
    return lookup


def collect_pmids_from_literature(literature: dict) -> list[str]:
    """literature dict에서 PMID 전체 수집 (중복 제거, 순서 유지)."""
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


# ── PubMed disease discovery ──────────────────────────────────


def _dedup_keep_order(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        value = str(value or "").strip()
        if value and value not in seen:
            seen.add(value)
            out.append(value)
    return out


def search_pubmed(
    query: str,
    email: str,
    max_results: int = 20,
    sort: str = "relevance",
) -> list[str]:
    """임의 PubMed query -> PMID 목록."""
    if not query.strip() or max_results <= 0:
        return []

    r = http_get(
        f"{EUTILS}/esearch.fcgi",
        params={
            "db": "pubmed",
            "term": query,
            "retmax": max_results,
            "sort": sort,
            "retmode": "json",
            "email": email,
        },
    )
    return r.json().get("esearchresult", {}).get("idlist", [])


def _quote_pubmed_phrase(text: str) -> str:
    """PubMed phrase query에 안전하게 넣기 위한 최소 escaping."""
    text = re.sub(r"\s+", " ", str(text or "").strip())
    text = text.replace('"', "")
    return text


def build_disease_query(
    disease_name: str,
    synonyms: list[str] | None = None,
    field: str = "Title/Abstract",
    max_synonyms: int = 8,
) -> str:
    """
    canonical disease name + synonym을 OR로 묶은 PubMed query 생성.

    너무 일반적인 짧은 synonym은 noise를 만들 수 있으므로 4자 미만은 제외한다.
    """
    names = [disease_name, *(synonyms or [])]
    cleaned: list[str] = []
    for name in names:
        phrase = _quote_pubmed_phrase(name)
        if len(phrase) < 4:
            continue
        if phrase.lower() not in {x.lower() for x in cleaned}:
            cleaned.append(phrase)
        if len(cleaned) >= max_synonyms + 1:
            break

    if not cleaned:
        return ""

    clauses = [f'"{name}"[{field}]' for name in cleaned]
    return "(" + " OR ".join(clauses) + ")"


def discover_disease_pmids(
    disease_name: str,
    synonyms: list[str] | None,
    email: str,
    bucket_limits: dict[str, int] | None = None,
    include_general: bool = True,
    sleep_seconds: float = 0.34,
) -> dict[str, dict[str, Any]]:
    """
    질환명만으로 PubMed discovery 검색을 수행한다.

    반환 예:
      {
        "review": {"query": "...", "pmids": [...]},
        "genotype_phenotype": {...},
        "genetics": {...},
        "cohort": {...},
        "general": {...},
      }

    목적은 모든 논문 수집이 아니라 gene/region/phenotype discovery용 anchor
    문헌을 제한적으로 확보하는 것이다.
    """
    limits = {
        "review": 5,
        "genotype_phenotype": 10,
        "genetics": 10,
        "cohort": 5,
        "general": 10,
    }
    if bucket_limits:
        limits.update(bucket_limits)

    disease_query = build_disease_query(disease_name, synonyms)
    if not disease_query:
        return {}

    common_filter = "English[lang] NOT comment[pt] NOT letter[pt]"

    queries: dict[str, str] = {
        "review": (
            f"{disease_query} AND (review[pt] OR systematic review[pt]) "
            f"AND {common_filter}"
        ),
        "genotype_phenotype": (
            f"{disease_query} AND ("
            'genotype[Title/Abstract] OR phenotype[Title/Abstract] OR '
            '"genotype-phenotype"[Title/Abstract] OR '
            '"critical region"[Title/Abstract] OR mapping[Title/Abstract]'
            f") AND {common_filter}"
        ),
        "genetics": (
            f"{disease_query} AND ("
            'gene[Title/Abstract] OR genes[Title/Abstract] OR '
            'genetic[Title/Abstract] OR haploinsufficiency[Title/Abstract] OR '
            'dosage[Title/Abstract] OR deletion[Title/Abstract] OR '
            'duplication[Title/Abstract] OR "copy number"[Title/Abstract]'
            f") AND {common_filter}"
        ),
        "cohort": (
            f"{disease_query} AND ("
            'cohort[Title/Abstract] OR patients[Title/Abstract] OR '
            '"case series"[Title/Abstract]'
            f") AND {common_filter}"
        ),
    }
    if include_general:
        queries["general"] = f"{disease_query} AND {common_filter}"

    result: dict[str, dict[str, Any]] = {}
    for bucket, query in queries.items():
        limit = max(0, int(limits.get(bucket, 0)))
        if limit == 0:
            continue
        try:
            pmids = search_pubmed(
                query=query,
                email=email,
                max_results=limit,
                sort="relevance",
            )
            result[bucket] = {"query": query, "pmids": pmids}
            print(f"  [PubMed discovery] {bucket}: {len(pmids)} PMID")
        except Exception as e:
            result[bucket] = {"query": query, "pmids": [], "error": str(e)}
            print(f"  [PubMed discovery] {bucket} 오류: {e}")
        if sleep_seconds > 0:
            time.sleep(sleep_seconds)

    return result


def collect_pmids_from_discovery(discovery: dict[str, dict[str, Any]]) -> list[str]:
    """discover_disease_pmids 결과를 bucket 순서대로 merge/dedup."""
    merged: list[str] = []
    for bucket_data in discovery.values():
        if isinstance(bucket_data, dict):
            merged.extend(bucket_data.get("pmids", []) or [])
    return _dedup_keep_order(merged)


def merge_medgen_and_discovery_pmids(
    medgen_literature: dict | None,
    discovery: dict[str, dict[str, Any]] | None,
) -> list[str]:
    """MedGen 노출 PMID + 독립 PubMed discovery PMID merge/dedup."""
    medgen_pmids = collect_pmids_from_literature(medgen_literature or {})
    discovery_pmids = collect_pmids_from_discovery(discovery or {})
    return _dedup_keep_order([*medgen_pmids, *discovery_pmids])


# ── Candidate gene discovery ──────────────────────────────────


def extract_candidate_genes(
    articles: dict[str, dict] | Iterable[dict],
    gene_symbols: Iterable[str],
    include_fields: tuple[str, ...] = ("title", "abstract", "keywords"),
    min_mentions: int = 1,
) -> dict[str, dict[str, Any]]:
    """
    확보된 PubMed article에서 gene symbol 후보를 추출한다.

    gene_symbols에는 질환 관련 gene 목록이 아니라 HGNC/NCBI Gene 등에서 받은
    '전체 허용 gene symbol 사전'을 넘기는 것을 전제로 한다.

    짧고 일반 단어와 충돌하기 쉬운 symbol도 존재하므로, 실제 운영에서는
    HGNC alias/withdrawn symbol 및 ambiguity blacklist를 함께 적용하는 것을 권장.
    """
    if isinstance(articles, dict):
        article_list = list(articles.values())
    else:
        article_list = list(articles)

    symbols = _dedup_keep_order(str(x).upper() for x in gene_symbols)
    # 긴 symbol 먼저 검사해 디버깅/출력 안정성 확보
    symbols.sort(key=lambda x: (-len(x), x))

    evidence: dict[str, dict[str, Any]] = {}

    for article in article_list:
        if not isinstance(article, dict) or article.get("error"):
            continue

        pieces: list[str] = []
        for field in include_fields:
            value = article.get(field, "")
            if isinstance(value, list):
                pieces.extend(str(v) for v in value)
            elif value:
                pieces.append(str(value))
        text = " ".join(pieces)
        if not text:
            continue

        pmid = str(article.get("pmid", ""))
        title = str(article.get("title", ""))

        for symbol in symbols:
            # 영숫자 gene symbol 경계. '-'가 포함된 일부 symbol도 허용.
            pattern = rf"(?<![A-Za-z0-9]){re.escape(symbol)}(?![A-Za-z0-9])"
            matches = re.findall(pattern, text, flags=re.I)
            if not matches:
                continue

            item = evidence.setdefault(
                symbol,
                {
                    "gene_symbol": symbol,
                    "mention_count": 0,
                    "paper_count": 0,
                    "pmids": [],
                    "papers": [],
                },
            )
            item["mention_count"] += len(matches)
            if pmid and pmid not in item["pmids"]:
                item["pmids"].append(pmid)
                item["paper_count"] += 1
                item["papers"].append(
                    {
                        "pmid": pmid,
                        "title": title,
                        "mentions": len(matches),
                    }
                )

    return {
        gene: data
        for gene, data in sorted(
            evidence.items(),
            key=lambda kv: (-kv[1]["paper_count"], -kv[1]["mention_count"], kv[0]),
        )
        if data["mention_count"] >= min_mentions
    }


def build_gene_disease_query(
    disease_name: str,
    synonyms: list[str] | None,
    gene_symbol: str,
) -> str:
    """후속 validation용 gene x disease PubMed query 생성."""
    disease_query = build_disease_query(disease_name, synonyms)
    symbol = re.sub(r"[^A-Za-z0-9_.-]", "", str(gene_symbol or "").strip())
    if not disease_query or not symbol:
        return ""
    return (
        f"{disease_query} AND {symbol}[Title/Abstract] "
        "AND English[lang] NOT comment[pt] NOT letter[pt]"
    )


def discover_gene_disease_pmids(
    disease_name: str,
    synonyms: list[str] | None,
    candidate_genes: Iterable[str],
    email: str,
    max_results_per_gene: int = 10,
    sleep_seconds: float = 0.34,
) -> dict[str, dict[str, Any]]:
    """candidate gene 각각에 대해 gene x disease validation PubMed 검색."""
    out: dict[str, dict[str, Any]] = {}
    for gene in _dedup_keep_order(str(g).upper() for g in candidate_genes):
        query = build_gene_disease_query(disease_name, synonyms, gene)
        if not query:
            continue
        try:
            pmids = search_pubmed(
                query=query,
                email=email,
                max_results=max_results_per_gene,
                sort="relevance",
            )
            out[gene] = {"query": query, "pmids": pmids}
            print(f"  [PubMed gene validation] {gene}: {len(pmids)} PMID")
        except Exception as e:
            out[gene] = {"query": query, "pmids": [], "error": str(e)}
        if sleep_seconds > 0:
            time.sleep(sleep_seconds)
    return out


def discover_and_fetch_disease_literature(
    disease_name: str,
    synonyms: list[str] | None,
    email: str,
    medgen_literature: dict | None = None,
    bucket_limits: dict[str, int] | None = None,
    verify_pmc: bool = False,
) -> dict[str, Any]:
    """
    convenience wrapper:
      disease -> discovery -> MedGen PMID merge -> EFetch

    discovery 단계에서는 verify_pmc=False를 기본으로 두어 API 호출량을 줄인다.
    full-text로 실제 사용할 anchor paper만 후속 단계에서 PMC 검증하는 것이 효율적이다.
    """
    discovery = discover_disease_pmids(
        disease_name=disease_name,
        synonyms=synonyms,
        email=email,
        bucket_limits=bucket_limits,
    )
    pmids = merge_medgen_and_discovery_pmids(
        medgen_literature=medgen_literature,
        discovery=discovery,
    )
    articles = fetch_pubmed_abstracts(
        pmids=pmids,
        email=email,
        verify_pmc=verify_pmc,
    )
    return {
        "disease_name": disease_name,
        "synonyms": synonyms or [],
        "discovery": discovery,
        "pmids": pmids,
        "articles": articles,
    }