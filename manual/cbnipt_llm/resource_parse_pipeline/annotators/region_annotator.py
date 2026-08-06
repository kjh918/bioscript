"""
annotators/region_annotator.py
────────────────────────────────
CoreRegion / PrimaryTargetRegion annotation:
  - ClinGen region dosage overlap
  - PubMed ESearch (critical_region / genotype_phenotype 쿼리)
"""

import re
import time
from typing import Any

from resource_parse_pipeline.utils import EUTILS, http_get


# ── PubMed 쿼리 ───────────────────────────────────────────────

def _build_region_queries(syndrome: str, feature_name: str,
                           feature_type: str, chrom: str) -> dict[str, str]:
    chrom_num = chrom.replace("chr", "").upper()
    queries: dict[str, str] = {}

    if feature_type in ("CoreRegion", "PrimaryTargetRegion"):
        queries["critical_region"] = (
            f'("{feature_name}"[tiab] OR "chromosome {chrom_num}"[tiab]) '
            f'AND ("{syndrome}"[MeSH Terms] OR "{syndrome}"[tiab]) '
            f'AND ("critical region"[tiab] OR "minimal region"[tiab] OR '
            f'"smallest region"[tiab] OR "deletion"[tiab])'
        )
        queries["genotype_phenotype"] = (
            f'("{feature_name}"[tiab]) '
            f'AND ("phenotype"[tiab] OR "genotype-phenotype"[tiab] OR '
            f'"clinical features"[tiab]) '
            f'AND ("{syndrome}"[MeSH Terms] OR "{syndrome}"[tiab])'
        )
    elif feature_type == "PartialChromosome":
        queries["partial_region"] = (
            f'("{feature_name}"[tiab] OR '
            f'"partial trisomy {chrom_num}"[tiab] OR '
            f'"partial monosomy {chrom_num}"[tiab]) '
            f'AND ("{syndrome}"[MeSH Terms] OR "{syndrome}"[tiab]) '
            f'AND ("phenotype"[tiab] OR "critical region"[tiab])'
        )

    return queries


def fetch_region_literature(
    syndrome: str,
    feature_name: str,
    feature_type: str,
    chrom: str,
    email: str,
    max_results: int = 10,
) -> dict[str, list[dict]]:
    """CoreRegion/PrimaryTargetRegion → 카테고리별 PubMed 논문"""
    import xml.etree.ElementTree as ET
    from resource_parse_pipeline.parsers.pubmed_parser import _parse_single_article

    queries = _build_region_queries(syndrome, feature_name, feature_type, chrom)
    literature: dict[str, list[dict]] = {}

    for qtype, query in queries.items():
        try:
            r = http_get(f"{EUTILS}/esearch.fcgi",
                         params={"db": "pubmed", "term": query,
                                 "retmax": max_results, "sort": "relevance",
                                 "retmode": "json", "email": email})
            pmids = r.json().get("esearchresult", {}).get("idlist", [])
            if not pmids:
                literature[qtype] = []
                continue

            time.sleep(0.35)
            r2 = http_get(f"{EUTILS}/efetch.fcgi",
                          params={"db": "pubmed", "id": ",".join(pmids),
                                  "rettype": "abstract", "retmode": "xml",
                                  "email": email})
            root = ET.fromstring(r2.text.encode("utf-8"))
            articles = [
                _parse_single_article(pa)
                for pa in root.findall(".//PubmedArticle")
            ]
            literature[qtype] = [a for a in articles if a]
            print(f"    [{qtype}] {len(literature[qtype])}편")

        except Exception as e:
            print(f"    [{qtype}] 오류: {e}")
            literature[qtype] = []
        time.sleep(0.35)

    return literature


# ── 배치 처리 ─────────────────────────────────────────────────

def annotate_regions(
    region_entries: list[dict],
    syndrome: str,
    email: str,
    clingen_region_db: list[dict] | None = None,
    max_lit: int = 10,
) -> list[dict]:
    """
    CoreRegion / PrimaryTargetRegion 목록 전체 annotation.
    ClinGen overlap + PubMed 논문 추가.
    """
    from resource_parse_pipeline.parsers.clingen_parser import get_region_overlaps

    for entry in region_entries:
        fname  = entry.get("feature_name", "")
        ftype  = entry.get("feature_type", "CoreRegion")
        chrom  = entry.get("chromosome", "")

        # ClinGen region overlap
        if clingen_region_db and chrom and entry.get("start") and entry.get("end"):
            overlaps = get_region_overlaps(
                chrom, entry["start"], entry["end"], clingen_region_db)
            entry["clingen_region"] = overlaps
        else:
            entry["clingen_region"] = []

        # PubMed 논문
        print(f"  [Region] '{fname}' 논문 수집")
        entry["feature_literature"] = fetch_region_literature(
            syndrome, fname, ftype, chrom, email, max_results=max_lit)

    return region_entries