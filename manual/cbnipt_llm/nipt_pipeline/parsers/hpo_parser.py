"""parsers/hpo_parser.py"""
from pathlib import Path
from collections import defaultdict

def load_hpo(path: str) -> dict[str, list[dict]]:
    lookup: dict[str, list[dict]] = defaultdict(list)
    if not path or not Path(path).exists():
        print(f"  [HPO] 파일 없음: {path}")
        return dict(lookup)
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            if parts[0].startswith("HP:"):
                hpo_id, hpo_name = parts[0], parts[1]
                gene_sym = parts[3] if len(parts) > 3 else ""
            else:
                gene_sym = parts[1] if len(parts) > 1 else ""
                hpo_id   = ""
                hpo_name = parts[4] if len(parts) > 4 else ""
            if gene_sym:
                lookup[gene_sym].append({"hpo_id": hpo_id, "hpo_name": hpo_name})
    print(f"  [HPO] {len(lookup)}개 유전자 로드")
    return dict(lookup)