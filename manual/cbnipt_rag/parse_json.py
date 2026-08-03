"""
parse_json.py
JSON 파일(질병 1개) → List[chunk dict] 변환

chunk 유형:
  - gene_meta    : gene_symbol + tier + MOI + cytoband + GenCC 분류
  - gene_summary : NCBI gene_summary 전문
  - abstract     : PubMed abstract 1개 = 1 chunk
"""

import json
import hashlib
from pathlib import Path


def _uid(*parts: str) -> str:
    """재현 가능한 고정 ID (upsert 시 중복 방지용)"""
    raw = "|".join(parts)
    return hashlib.md5(raw.encode()).hexdigest()


def parse_disease_json(json_path: str) -> list[dict]:
    """
    Returns:
        List of {
            "id"      : str  (MD5, Qdrant point id 용),
            "text"    : str  (임베딩 대상),
            "metadata": dict (Qdrant payload),
        }
    """
    data = json.loads(Path(json_path).read_text())
    syndrome = data.get("syndrome", "unknown")
    chunks = []

    for gene in data.get("genes", []):
        symbol      = gene.get("gene_symbol", "")
        tier        = gene.get("evidence_tier", "")
        moi         = gene.get("gencc_moi", "")
        cytoband    = gene.get("cytoband", "")
        gencc_class = gene.get("gencc_classification", "")
        submitters  = gene.get("gencc_submitters", "")
        gencc_pmids = gene.get("gencc_pmids", "")
        ncbi        = gene.get("ncbi_gene", {})
        ncbi_id     = ncbi.get("ncbi_gene_id", "")

        # ── chunk 1: gene_meta ──────────────────────────────────────
        meta_text = (
            f"Gene: {symbol}\n"
            f"Syndrome: {syndrome}\n"
            f"Evidence tier: {tier}\n"
            f"GenCC classification: {gencc_class}\n"
            f"Mode of inheritance: {moi}\n"
            f"Cytoband: {cytoband}\n"
            f"GenCC submitters: {submitters}\n"
            f"GenCC PMIDs: {gencc_pmids}"
        )
        chunks.append({
            "id": _uid(syndrome, symbol, "gene_meta"),
            "text": meta_text,
            "metadata": {
                "chunk_type":  "gene_meta",
                "syndrome":    syndrome,
                "gene_symbol": symbol,
                "tier":        tier,
                "moi":         moi,
                "cytoband":    cytoband,
                "ncbi_gene_id": ncbi_id,
            },
        })

        # ── chunk 2: gene_summary ───────────────────────────────────
        summary = ncbi.get("gene_summary", "").strip()
        if summary:
            summary_text = (
                f"Gene: {symbol} ({ncbi.get('gene_full_name', '')})\n"
                f"Syndrome: {syndrome}\n"
                f"NCBI Gene Summary:\n{summary}"
            )
            chunks.append({
                "id": _uid(syndrome, symbol, "gene_summary"),
                "text": summary_text,
                "metadata": {
                    "chunk_type":  "gene_summary",
                    "syndrome":    syndrome,
                    "gene_symbol": symbol,
                    "ncbi_gene_id": ncbi_id,
                },
            })

        # ── chunk 3: abstracts ──────────────────────────────────────
        for paper in gene.get("literature", {}).get("gene_disease", []):
            pmid     = paper.get("pmid", "")
            title    = paper.get("title", "")
            abstract = paper.get("abstract", "").strip()
            pmc_id   = paper.get("pmc_id", "")

            if not abstract:
                continue

            abs_text = (
                f"Syndrome: {syndrome}\n"
                f"Gene: {symbol}\n"
                f"Title: {title}\n"
                f"Abstract:\n{abstract}"
            )
            chunks.append({
                "id": _uid(syndrome, symbol, pmid),
                "text": abs_text,
                "metadata": {
                    "chunk_type":  "abstract",
                    "syndrome":    syndrome,
                    "gene_symbol": symbol,
                    "pmid":        pmid,
                    "pmc_id":      pmc_id,
                    "title":       title,
                    "full_text":   paper.get("full_text_available", False),
                },
            })

    return chunks


# ── CLI 실행 시 검증 출력 ────────────────────────────────────────────
if __name__ == "__main__":
    import sys, textwrap

    json_path = sys.argv[1] if len(sys.argv) > 1 else "data/digeorge_syndrome.json"
    chunks = parse_disease_json(json_path)

    print(f"\n✅ 총 {len(chunks)}개 chunk 생성\n")
    for c in chunks:
        print(f"[{c['metadata']['chunk_type']:12s}] id={c['id'][:8]}… "
              f"| {c['metadata'].get('gene_symbol','')} "
              f"| pmid={c['metadata'].get('pmid','')}")
        preview = textwrap.shorten(c["text"], width=120, placeholder="…")
        print(f"  {preview}\n")