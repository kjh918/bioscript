"""
annotators/chromosome_annotator.py
────────────────────────────────────
Autosome Abnormality / Sex Chromosome Abnormality 처리.

TargetChromosome feature_type → GRCh38 전체 좌표 마커 + 질환 특이 annotation.

기존 파이프라인은 CoreGene 기반이라 TargetChromosome에 아무것도 붙지 않았음.
이 모듈이 그 공백을 채운다:
  1. GRCh38 전체 염색체 좌표 (하드코딩)
  2. 질환 타입별 annotation (trisomy/monosomy/sex_chr)
  3. ClinGen region overlap (해당 염색체 전체 검색)
  4. 임상 마커 정보 (PAR 영역, 염색체 특이 gene density 등)
"""

import re
from typing import Any

from resource_parse_pipeline.utils import get_chrom_marker, GRCH38_CHROM_COORDS

# ── 질환별 타입 분류 ──────────────────────────────────────────
# nipt_id → {mechanism, karyotype, affected_chromosomes, clinical_notes}
CHROMOSOMAL_DISEASE_META: dict[str, dict] = {
    # Autosome Abnormality
    "NIPT_TRISOMY9": {
        "mechanism":            "trisomy",
        "karyotype":            "47,+9",
        "affected_chromosomes": ["chr9"],
        "copy_number":          3,
        "clinical_notes": (
            "Trisomy 9 is rare, usually mosaic. Features include "
            "distinctive facial features, cardiac defects, and joint abnormalities. "
            "Full trisomy 9 is typically lethal."
        ),
    },
    "NIPT_TRISOMY13": {
        "mechanism":            "trisomy",
        "karyotype":            "47,+13",
        "affected_chromosomes": ["chr13"],
        "copy_number":          3,
        "clinical_notes": (
            "Patau syndrome. Severe intellectual disability, heart defects, "
            "cleft lip/palate, polydactyly. Most die within first year."
        ),
    },
    "NIPT_TRISOMY18": {
        "mechanism":            "trisomy",
        "karyotype":            "47,+18",
        "affected_chromosomes": ["chr18"],
        "copy_number":          3,
        "clinical_notes": (
            "Edwards syndrome. Severe intellectual disability, heart defects, "
            "clenched fists, rocker-bottom feet. 90% die within first year."
        ),
    },
    "NIPT_TRISOMY21": {
        "mechanism":            "trisomy",
        "karyotype":            "47,+21",
        "affected_chromosomes": ["chr21"],
        "copy_number":          3,
        "clinical_notes": (
            "Down syndrome. Most common chromosomal aneuploidy. "
            "Intellectual disability, characteristic facial features, "
            "hypotonia. HR-DSCR at 21q22.13 is the critical region."
        ),
    },
    "NIPT_TRISOMY22": {
        "mechanism":            "trisomy",
        "karyotype":            "47,+22",
        "affected_chromosomes": ["chr22"],
        "copy_number":          3,
        "clinical_notes": (
            "Trisomy 22 mosaicism syndrome. Variable phenotype. "
            "Features include growth retardation, preauricular skin tags, "
            "cardiac defects, and genital abnormalities."
        ),
    },
    # Sex Chromosome Abnormality
    "NIPT_45X": {
        "mechanism":            "monosomy",
        "karyotype":            "45,X",
        "affected_chromosomes": ["chrX"],
        "copy_number":          1,
        "sex_chromosome_type":  "X_monosomy",
        "clinical_notes": (
            "Turner syndrome. Short stature, gonadal dysgenesis, "
            "webbed neck, cardiovascular defects (bicuspid aortic valve, "
            "coarctation of aorta). 45,X often mosaic."
        ),
        "par_regions": {
            "PAR1": {"chrom": "chrX", "start": 60001,    "end": 2699520},
            "PAR2": {"chrom": "chrX", "start": 154931044,"end": 155260560},
        },
        "key_regions": {
            "SHOX":  {"chrom": "chrX", "start": 585079, "end": 619522,
                      "note": "Short stature homeobox gene, PAR1 region"},
            "XIST":  {"chrom": "chrX", "start": 73820650, "end": 73852753,
                      "note": "X-inactivation specific transcript"},
        },
    },
    "NIPT_47XXY": {
        "mechanism":            "sex_chromosome_aneuploidy",
        "karyotype":            "47,XXY",
        "affected_chromosomes": ["chrX", "chrY"],
        "copy_number_X":        2,
        "copy_number_Y":        1,
        "sex_chromosome_type":  "XXY",
        "clinical_notes": (
            "Klinefelter syndrome. Hypogonadism, infertility, tall stature, "
            "gynecomastia, learning difficulties. Extra X chromosome. "
            "1 in 500-1000 male births."
        ),
        "par_regions": {
            "PAR1": {"chrom": "chrX", "start": 60001, "end": 2699520},
            "PAR2": {"chrom": "chrX", "start": 154931044, "end": 155260560},
        },
        "key_regions": {
            "SHOX":  {"chrom": "chrX", "start": 585079, "end": 619522,
                      "note": "Dosage-sensitive, PAR1"},
            "KAL1":  {"chrom": "chrX", "start": 865253, "end": 1069966,
                      "note": "Xp22.31, anosmin-1"},
        },
    },
    "NIPT_47XYY": {
        "mechanism":            "sex_chromosome_aneuploidy",
        "karyotype":            "47,XYY",
        "affected_chromosomes": ["chrY"],
        "copy_number_Y":        2,
        "sex_chromosome_type":  "XYY",
        "clinical_notes": (
            "Double Y syndrome. Tall stature, learning difficulties, "
            "behavioral issues. Usually not identified until adulthood. "
            "Extra Y chromosome."
        ),
        "par_regions": {
            "PAR1": {"chrom": "chrY", "start": 10001,    "end": 2781479},
            "PAR2": {"chrom": "chrY", "start": 56887903, "end": 57217415},
        },
        "key_regions": {
            "SRY":   {"chrom": "chrY", "start": 2786855, "end": 2787699,
                      "note": "Sex-determining region Y"},
            "AMELY": {"chrom": "chrY", "start": 6865936, "end": 6914948,
                      "note": "Amelogenin Y, used in sex determination"},
        },
    },
}


def annotate_target_chromosome(
    nipt_id: str,
    chrom: str,
    clingen_regions: list[dict] | None = None,
) -> dict[str, Any]:
    """
    TargetChromosome feature → 전체 annotation 반환.

    Parameters
    ----------
    nipt_id          : NIPT_TRISOMY21 등
    chrom            : chr21, chrX 등
    clingen_regions  : parse_clingen_region() 결과 (optional)
    """
    result: dict[str, Any] = {
        "feature_name": chrom,
        "feature_type": "TargetChromosome",
        "chromosome":   chrom,
    }

    # GRCh38 좌표
    marker = get_chrom_marker(chrom)
    result["grch38_coords"] = marker

    # 질환 특이 meta
    meta = CHROMOSOMAL_DISEASE_META.get(nipt_id, {})
    if meta:
        result["mechanism"]             = meta.get("mechanism", "")
        result["karyotype"]             = meta.get("karyotype", "")
        result["copy_number"]           = meta.get("copy_number", "")
        result["clinical_notes"]        = meta.get("clinical_notes", "")
        result["sex_chromosome_type"]   = meta.get("sex_chromosome_type", "")
        result["par_regions"]           = meta.get("par_regions", {})
        result["key_regions"]           = meta.get("key_regions", {})
        result["affected_chromosomes"]  = meta.get("affected_chromosomes", [chrom])
    else:
        # meta 없으면 기본값
        result["mechanism"]             = _infer_mechanism(nipt_id)
        result["affected_chromosomes"]  = [chrom]

    # ClinGen region overlap (해당 염색체 전체)
    if clingen_regions and marker:
        overlaps = [
            r for r in clingen_regions
            if r.get("chrom") == chrom and
            r.get("start", 0) <= marker["end"] and
            r.get("end", 0) >= marker["start"]
        ]
        # dosage sensitive region만 필터 (hi_score >= 2)
        dosage_sensitive = [
            r for r in overlaps
            if r.get("hi_score", "0") in ("2", "3")
        ]
        result["clingen_overlaps"]       = len(overlaps)
        result["clingen_dosage_regions"] = dosage_sensitive[:10]
    else:
        result["clingen_overlaps"]       = 0
        result["clingen_dosage_regions"] = []

    return result


def annotate_all_chromosomal(
    features: dict,
    nipt_id: str,
    clingen_regions: list[dict] | None = None,
) -> list[dict]:
    """
    parse_markerset_features() 결과에서
    target_chromosomes + partial_chromosomes 모두 처리.
    """
    annotated = []

    for entry in features.get("target_chromosomes", []):
        chrom = entry.get("chromosome", entry.get("feature_name", ""))
        ann = annotate_target_chromosome(nipt_id, chrom, clingen_regions)
        ann.update({k: v for k, v in entry.items()
                    if k not in ann})  # 마커셋 원본 필드 보존
        annotated.append(ann)

    for entry in features.get("partial_chromosomes", []):
        # PartialChromosome: 18p, 18q12-q21 등 → 좌표 추정
        fname = entry.get("feature_name", "")
        chrom = entry.get("chromosome", "")
        ann = {
            "feature_name":  fname,
            "feature_type":  "PartialChromosome",
            "chromosome":    chrom,
            "mechanism":     "partial_" + _infer_mechanism(nipt_id),
            "karyotype":     f"partial {nipt_id}",
            "grch38_coords": _partial_coords(fname, chrom),
        }
        ann.update({k: v for k, v in entry.items()
                    if k not in ann})
        annotated.append(ann)

    return annotated


def _infer_mechanism(nipt_id: str) -> str:
    nid = nipt_id.upper()
    if "TRISOMY" in nid:
        return "trisomy"
    if "45X" in nid or "MONOSOMY" in nid:
        return "monosomy"
    if "47XX" in nid or "47XY" in nid:
        return "sex_chromosome_aneuploidy"
    return "chromosomal_abnormality"


def _partial_coords(feature_name: str, chrom: str) -> dict:
    """
    18p, 18q12-q21 같은 partial chromosome 표기 → 좌표 추정.
    정확한 좌표는 UCSC나 Ensembl에서 가져와야 하지만
    여기서는 arm 기반 대략 좌표 반환.
    """
    from nipt_pipeline.utils import GRCH38_CHROM_COORDS, GRCH38_CENTROMERE
    key = chrom if chrom.startswith("chr") else f"chr{chrom}"
    coords = GRCH38_CHROM_COORDS.get(key, {})
    cen    = GRCH38_CENTROMERE.get(key, 0)
    if not coords:
        return {}

    fn = feature_name.lower()
    # p arm
    if re.search(r"^chr?\d+p$", fn, re.I) or fn.endswith("p"):
        return {"start": 1, "end": cen,
                "arm": "p", "note": "p arm (approximate)"}
    # q arm
    if re.search(r"^chr?\d+q$", fn, re.I) or fn.endswith("q"):
        return {"start": cen + 1, "end": coords["end"],
                "arm": "q", "note": "q arm (approximate)"}
    # specific band (18q12-q21 등)
    return {"start": 1, "end": coords["end"],
            "note": f"partial region {feature_name} (coordinate not resolved)"}
