"""parsers/gencc_parser.py - GenCC TSV 파서"""
import csv
from pathlib import Path
from collections import defaultdict

def load_gencc(tsv_path: str) -> dict[str, list[dict]]:
    """GenCC TSV → {gene_symbol: [entries]}"""
    lookup: dict[str, list[dict]] = defaultdict(list)
    if not tsv_path or not Path(tsv_path).exists():
        return lookup
    with open(tsv_path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sym = row.get("gene_symbol", "").strip()
            if not sym:
                continue
            lookup[sym].append({
                "disease_curie":  row.get("disease_curie","").strip(),
                "disease_title":  row.get("disease_title","").strip(),
                "classification": row.get("classification_title","").strip(),
                "moi":            row.get("moi_title","").strip(),
                "submitter":      row.get("submitter_title","").strip(),
                "pmids":          row.get("submitted_as_pmids","").strip(),
                "report_url":     row.get("submitted_as_public_report_url","").strip(),
            })
    print(f"  [GenCC] {len(lookup)}개 유전자 로드")
    return dict(lookup)