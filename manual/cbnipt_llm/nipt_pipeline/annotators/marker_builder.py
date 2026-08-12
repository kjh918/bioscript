"""
annotators/marker_builder.py
──────────────────────────────
disease_info → markers[]

MedGen 파싱 결과를 기반으로 마커 후보를 구성한다.
Genomic position은 references/hg38.json의 chromosome/cytoband 정보를 사용한다.
NCBI Gene API는 사용하지 않는다.

INPUT : disease_info (medgen_parser 출력)
        optional: clingen_gene_db, gencc_db, hpo_db (dict)
OUTPUT: markers list
"""

import sys
import json
import argparse
import pathlib
from typing import Any

from nipt_pipeline.utils import cytoband_pos, full_chrom_pos


# ── annotation ──────────────────────────────────────────────────

def _annotation_from_dbs(
    gene_symbol: str,
    clingen_gene_db: dict,
    gencc_db: dict,
    hpo_db: dict,
) -> dict:
    """외부 DB에서 gene annotation 구성."""
    ann: dict[str, Any] = {}

    cd = clingen_gene_db.get(gene_symbol, {})
    if cd:
        ann["hi_score"] = cd.get("hi_score") or None
        ann["hi_score_label"] = cd.get("hi_score_label") or None
        ann["hi_disease"] = cd.get("hi_disease_name") or None
        ann["ts_score"] = cd.get("ts_score") or None
        ann["hi_pmids"] = cd.get("hi_pmids") or []
        ann["hgnc_id"] = cd.get("hgnc_id") or None

    gc_list = gencc_db.get(gene_symbol, [])
    if gc_list:
        best = max(
            gc_list,
            key=lambda g: {
                "Definitive": 4,
                "Strong": 3,
                "Moderate": 2,
                "Limited": 1,
            }.get(g.get("classification", ""), 0),
        )
        ann["gencc_classification"] = best.get("classification") or None
        ann["gencc_moi"] = best.get("moi") or None
        ann["gencc_submitters"] = list(
            {
                g.get("submitter", "")
                for g in gc_list
                if g.get("submitter")
            }
        )

    hpo_terms = hpo_db.get(gene_symbol, [])[:15]
    if hpo_terms:
        ann["hpo_terms"] = hpo_terms

    return ann


def _gene_confidence(
    gene_symbol: str,
    clingen_gene_db: dict,
    gencc_db: dict,
) -> str:
    """ClinGen / GenCC 근거를 이용한 CoreGene confidence."""
    hi = str(clingen_gene_db.get(gene_symbol, {}).get("hi_score", ""))

    gc_list = gencc_db.get(gene_symbol, [])
    gc_classes = {
        g.get("classification", "")
        for g in gc_list
        if g.get("classification")
    }

    if hi == "3" or gc_classes & {"Definitive", "Strong"}:
        return "high"
    if hi == "2" or "Moderate" in gc_classes:
        return "medium"

    # MedGen Gene(location) 직접 명시 유전자는 기본 high
    return "high"


def _normalize_gene_cytoband(cytoband: str | None) -> tuple[str | None, str | None]:
    """
    cytoband 문자열에서 chromosome과 band를 분리한다.

    Examples
    --------
    22q11.21 -> ("chr22", "22q11.21")
    chr22q11.21 -> ("chr22", "22q11.21")
    Xp22.31 -> ("chrX", "Xp22.31")
    """
    if not cytoband:
        return None, None

    value = str(cytoband).strip()
    if not value:
        return None, None

    import re

    m = re.match(r"^(?:chr)?(\d+|X|Y)([pq].+)$", value, re.I)
    if not m:
        return None, value

    chrom_no = m.group(1).upper()
    band_part = m.group(2)
    return f"chr{chrom_no}", f"{chrom_no}{band_part}"


# ── marker builder ──────────────────────────────────────────────

def build_markers(
    disease_info: dict,
    clingen_gene_db: dict | None = None,
    gencc_db: dict | None = None,
    hpo_db: dict | None = None,
) -> list[dict]:
    """
    disease_info → markers list.

    CoreGene:
        MedGen Gene(location) 직접 명시 유전자.
        genomic_pos는 gene 자체 좌표가 아니라 해당 cytoband 범위를 사용한다.

    CoreRegion:
        genomic_targets 중 cytogenetic band.
        references/hg38.json의 cytobands에서 좌표를 찾는다.
        정확한 band가 없으면 utils.cytoband_pos()가 하위 band를 병합한다.

    TargetChromosome:
        genomic_targets 중 전체 염색체.
        references/hg38.json의 chromosome 전체 좌표를 사용한다.
    """
    cg_db = clingen_gene_db or {}
    gc_db = gencc_db or {}
    hp_db = hpo_db or {}

    markers: list[dict] = []

    # ── CoreGene ────────────────────────────────────────────────
    for gene_loc in disease_info.get("gene_locations", []):
        sym = gene_loc.get("gene_symbol")
        if not sym:
            continue

        cytoband = gene_loc.get("cytoband")
        chrom, normalized_cytoband = _normalize_gene_cytoband(cytoband)

        print(f"  [marker] CoreGene: {sym}")

        geo_pos = None
        if chrom and normalized_cytoband:
            geo_pos = cytoband_pos(chrom, normalized_cytoband)
            if geo_pos is None:
                print(
                    f"    [marker] cytoband 좌표 없음: "
                    f"{normalized_cytoband}"
                )

        ann = _annotation_from_dbs(sym, cg_db, gc_db, hp_db)

        # MedGen parser가 이미 제공한 ID가 있으면 보존
        if gene_loc.get("ncbi_gene_id"):
            ann["ncbi_gene_id"] = gene_loc.get("ncbi_gene_id")

        confidence = _gene_confidence(sym, cg_db, gc_db)

        source = ["medgen_direct"]
        if cg_db.get(sym):
            source.append("ClinGen")
        if gc_db.get(sym):
            source.append("GenCC")
        if hp_db.get(sym):
            source.append("HPO")
        if geo_pos:
            source.append("hg38_cytoband")

        markers.append(
            {
                "feature_type": "CoreGene",
                "feature_name": sym,
                "chromosome": chrom,
                "cytoband": normalized_cytoband or cytoband,
                "genomic_pos": geo_pos,
                "event": None,
                "mechanism": None,
                "confidence": confidence,
                "source": source,
                "annotation": ann if ann else None,
            }
        )

    # ── CoreRegion / TargetChromosome ──────────────────────────
    for tgt in disease_info.get("genomic_targets", []):
        ftype = tgt.get("type")
        chrom = tgt.get("chromosome")
        band = tgt.get("cytoband")
        ev = tgt.get("event")
        mech = tgt.get("mechanism")

        if not ftype or not chrom:
            continue

        source = [tgt.get("source", "name_or_synonym")]

        print(f"  [marker] {ftype}: {band or chrom}")

        if ftype == "TargetChromosome":
            geo_pos = full_chrom_pos(chrom)
            if geo_pos:
                source.append("hg38_chromosome")

        elif ftype == "CoreRegion":
            geo_pos = cytoband_pos(chrom, band) if band else None
            if geo_pos:
                source.append("hg38_cytoband")
            elif band:
                print(f"    [marker] cytoband 좌표 없음: {band}")

        else:
            geo_pos = None

        markers.append(
            {
                "feature_type": ftype,
                "feature_name": chrom if ftype == "TargetChromosome" else band,
                "chromosome": chrom,
                "cytoband": band,
                "genomic_pos": geo_pos,
                "event": ev,
                "mechanism": mech,
                "confidence": "medium",
                "source": source,
                "annotation": None,
            }
        )

    # confidence 순 정렬: high → medium → low
    order = {"high": 0, "medium": 1, "low": 2}
    markers.sort(
        key=lambda x: (
            order.get(x.get("confidence"), 3),
            x.get("feature_type") != "CoreGene",
        )
    )

    return markers


# ── CLI ─────────────────────────────────────────────────────────

def _main():
    ap = argparse.ArgumentParser(
        description="disease_info → markers JSON (local hg38 cytoband 기반)",
        epilog=(
            "python -m nipt_pipeline.annotators.marker_builder "
            "--disease-info disease_info.json --output markers.json"
        ),
    )

    ap.add_argument(
        "--disease-info",
        required=True,
        help="medgen_parser 출력 JSON",
    )
    ap.add_argument(
        "--clingen-gene",
        default="",
        help="ClinGen_gene_curation_list_GRCh38.tsv",
    )
    ap.add_argument(
        "--gencc",
        default="",
        help="gencc_submissions.tsv",
    )
    ap.add_argument(
        "--hpo",
        default="",
        help="phenotype_to_genes.txt",
    )
    ap.add_argument(
        "--output",
        default="-",
    )

    args = ap.parse_args()

    disease_info = json.loads(
        pathlib.Path(args.disease_info).read_text(encoding="utf-8")
    )

    cg_db: dict = {}
    gc_db: dict = {}
    hp_db: dict = {}

    if args.clingen_gene:
        from nipt_pipeline.parsers.clingen_parser import parse_clingen_gene
        cg_db = parse_clingen_gene(args.clingen_gene)

    if args.gencc:
        from nipt_pipeline.parsers.gencc_parser import load_gencc
        gc_db = load_gencc(args.gencc)

    if args.hpo:
        from nipt_pipeline.parsers.hpo_parser import load_hpo
        hp_db = load_hpo(args.hpo)

    markers = build_markers(
        disease_info=disease_info,
        clingen_gene_db=cg_db,
        gencc_db=gc_db,
        hpo_db=hp_db,
    )

    out = json.dumps(markers, ensure_ascii=False, indent=2)

    if args.output == "-":
        print(out)
    else:
        pathlib.Path(args.output).write_text(out, encoding="utf-8")
        print(
            f"저장: {args.output} ({len(markers)}개 마커)",
            file=sys.stderr,
        )


if __name__ == "__main__":
    _main()