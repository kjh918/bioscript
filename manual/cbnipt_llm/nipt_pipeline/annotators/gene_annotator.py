"""
annotators/gene_annotator.py
──────────────────────────────
CoreGene annotation:
  - NCBI Gene (ESearch + ESummary): gene_id, full_name, summary, map_location
  - UniProt: protein_function, domains, disease_associations
  - Human 검증 3단계 방어 (마우스 ortholog 혼입 방지)
"""

import re
import time
from typing import Any

from resource_parse_pipeline.utils import EUTILS, http_get


# ── NCBI Gene ─────────────────────────────────────────────────

def fetch_ncbi_gene(gene_symbol: str, email: str,
                    hgnc_id: str = "",
                    cytoband: str = "") -> dict[str, Any]:
    """
    NCBI Gene에서 Human(taxid=9606) 유전자만 조회.

    방어 전략:
      1순위: HGNC ID → Gene Accession 직접 검색
      2순위: Gene Name + Organism + 염색체 번호 필터
      3순위: Gene Name + Organism → taxid=9606 검증
    """
    def _esummary(gene_id: str) -> dict:
        time.sleep(0.35)
        r = http_get(f"{EUTILS}/esummary.fcgi",
                     params={"db": "gene", "id": gene_id,
                             "retmode": "json", "email": email})
        return r.json().get("result", {}).get(gene_id, {})

    def _is_human(data: dict) -> bool:
        org = data.get("organism", {})
        return (org.get("taxid", 0) == 9606 or
                "homo sapiens" in org.get("scientificname", "").lower())

    def _build(gene_id: str, data: dict) -> dict:
        ginfo = (data.get("genomicinfo", [{}])[0]
                 if data.get("genomicinfo") else {})
        return {
            "ncbi_gene_id":   gene_id,
            "gene_full_name": data.get("description", ""),
            "gene_summary":   data.get("summary", ""),
            "hgnc_id":        data.get("nomenclatureauthorityid", ""),
            "map_location":   data.get("maplocation", ""),
            "chromosome":     data.get("chromosome", ""),
            "chr_start":      ginfo.get("chrstart", ""),
            "chr_stop":       ginfo.get("chrstop", ""),
            "chr_accver":     ginfo.get("chraccver", ""),
            "ncbi_gene_url":  f"https://www.ncbi.nlm.nih.gov/gene/{gene_id}",
        }

    try:
        # 1순위: HGNC ID
        hgnc_num = re.sub(r"^HGNC:", "", hgnc_id).strip() if hgnc_id else ""
        if hgnc_num:
            r = http_get(f"{EUTILS}/esearch.fcgi",
                         params={"db": "gene",
                                 "term": f"HGNC:{hgnc_num}[Gene Accession]",
                                 "retmax": 1, "retmode": "json",
                                 "email": email})
            ids = r.json().get("esearchresult", {}).get("idlist", [])
            if ids:
                data = _esummary(ids[0])
                if _is_human(data):
                    print(f"    [NCBI Gene] {gene_symbol} → HGNC:{hgnc_num} "
                          f"id={ids[0]} loc={data.get('maplocation','')}")
                    return _build(ids[0], data)

        # 2순위: Gene Name + Organism + CHR
        chrom_num = ""
        if cytoband:
            m = re.match(r"^(\d+|X|Y)", str(cytoband).replace("chr", ""))
            if m:
                chrom_num = m.group(1)

        if chrom_num:
            query = (f"{gene_symbol}[Gene Name] AND "
                     f"Homo sapiens[Organism] AND {chrom_num}[CHR]")
            r = http_get(f"{EUTILS}/esearch.fcgi",
                         params={"db": "gene", "term": query,
                                 "retmax": 5, "sort": "relevance",
                                 "retmode": "json", "email": email})
            for gene_id in r.json().get("esearchresult", {}).get("idlist", []):
                data = _esummary(gene_id)
                if _is_human(data):
                    print(f"    [NCBI Gene] {gene_symbol} → chr{chrom_num} "
                          f"id={gene_id} loc={data.get('maplocation','')}")
                    return _build(gene_id, data)

        # 3순위: Gene Name + Organism → taxid 검증
        r = http_get(f"{EUTILS}/esearch.fcgi",
                     params={"db": "gene",
                             "term": (f"{gene_symbol}[Gene Name] AND "
                                      f"Homo sapiens[Organism]"),
                             "retmax": 5, "sort": "relevance",
                             "retmode": "json", "email": email})
        for gene_id in r.json().get("esearchresult", {}).get("idlist", []):
            data = _esummary(gene_id)
            if _is_human(data):
                print(f"    [NCBI Gene] {gene_symbol} → organism filter "
                      f"id={gene_id} loc={data.get('maplocation','')}")
                return _build(gene_id, data)

        print(f"    [NCBI Gene] {gene_symbol} → Human 결과 없음")
        return {}

    except Exception as e:
        print(f"    [NCBI Gene] {gene_symbol} 오류: {e}")
        return {}


# ── 배치 처리 ─────────────────────────────────────────────────

def annotate_core_genes(
    gene_entries: list[dict],
    email: str,
    uniprot_local_db=None,
) -> list[dict]:
    """
    CoreGene 목록 전체 annotation.

    Parameters
    ----------
    gene_entries   : parse_markerset_features()['core_genes'] 목록
    email          : NCBI API 이메일
    uniprot_local_db : UniProtLocalDB 인스턴스 (없으면 REST/ExPASy 사용)

    각 entry에 ncbi_gene, uniprot 필드를 추가해서 반환.
    """
    from resource_parse_pipeline.parsers.uniprot_parser import fetch_uniprot

    for entry in gene_entries:
        sym      = entry.get("feature_name", "")
        hgnc_id  = entry.get("hgnc_id", "")
        cytoband = entry.get("cytoband", "") or entry.get("size_mb", "")

        print(f"  → {sym}  hgnc={hgnc_id}  cytoband={cytoband}")

        # NCBI Gene
        ncbi = fetch_ncbi_gene(sym, email,
                               hgnc_id=hgnc_id, cytoband=str(cytoband))
        entry["ncbi_gene"] = ncbi

        # cytoband 보완
        if not entry.get("cytoband") and ncbi.get("map_location"):
            entry["cytoband"] = ncbi["map_location"]

        # UniProt
        entry["uniprot"] = fetch_uniprot(
            sym,
            local_db=uniprot_local_db,
            email=email,
        )

        time.sleep(0.5)

    return gene_entries