"""
sex_calling.py
──────────────
태아 성별 판별과 성염색체 이수성(Turner/Klinefelter/Jacob(XYY)/TripleX 등)
판정을 하나의 모듈로 통합합니다.

[통합 배경]
기존에는 동일한 판정이 서로 다른 스케일의 기준으로 두 군데에 따로 있었습니다:
  - chrom_calling.analyze_all_chromosomes(): raw_count 비율 기반 하드컷오프
  - clinical_diagnosis.diagnose_clinical_markers(): copy_number_signal(CN) 기반 range

이 모듈은 후자(CN scale, 상염색체 baseline=2.0)를 단일 기준으로 채택하고,
config.yaml 의 CFG["SEX_ANEUPLOIDY"] 테이블을 유일한 근거로 사용합니다.
모든 호출부는 x_cn / y_cn (copy_number_signal 기준 chrX/chrY 중앙값)을
동일한 스케일로 전달해야 합니다.

[동작 변경 안내]
과거 raw_count 비율 기준 threshold(sex_xxx_ratio 등)와 CN 기준 threshold는
정확히 같은 값이 아니었기 때문에, 통합 이후 일부 경계값 근처 샘플의 판정이
달라질 수 있습니다. 검증이 필요합니다.
"""

import numpy as np

from rules import CFG


def _normalize_sex_tag(sex_tag):
    """'Male (XY)' / 'Female (XX)' / 'XY' / 'XX' 등 다양한 표기를 'XY'/'XX'/'Unknown'으로 통일."""
    if sex_tag in ("XY", "Male (XY)", "Male"):
        return "XY"
    if sex_tag in ("XX", "Female (XX)", "Female"):
        return "XX"
    return "Unknown"


def determine_fetal_sex(x_cn, y_cn, cfg=CFG):
    """
    copy_number_signal 기준 chrX/chrY 중앙값(x_cn, y_cn)으로 태아 성별 판별.
    (normalization.normalize_by_chrom_with_sex, pipeline._determine_fetal_sex
     에 각각 하드코딩되어 있던 1.5 / 0.2 기준을 config로 이동)
    """
    x_female_min = cfg["x_female_min_cn"]
    y_male_min = cfg["y_male_min_cn"]

    if x_cn > x_female_min and y_cn < y_male_min:
        return "XX"
    if x_cn < x_female_min and y_cn > y_male_min:
        return "XY"
    return "Unknown"


def classify_sex_chromosome(chrom, sex_tag, x_cn, y_cn, cfg=CFG):
    """
    chrX / chrY 에 대한 이수성 하드콜 판정.

    Parameters
    ----------
    chrom   : "chrX" 또는 "chrY" (그 외 염색체는 항상 None 반환)
    sex_tag : determine_fetal_sex() 등이 반환한 성별 태그 (다양한 표기 허용)
    x_cn, y_cn : copy_number_signal 기준 chrX / chrY 중앙값 (상염색체 baseline=2.0)

    Returns
    -------
    dict(call, pattern_name, score) 매칭되는 규칙이 있을 때, 없으면 None
    """
    if chrom not in ("chrX", "chrY"):
        return None

    sex_norm = _normalize_sex_tag(sex_tag)
    if sex_norm == "Unknown":
        return None

    for _, rule in cfg["SEX_ANEUPLOIDY"].items():
        if rule["applies_to"] != sex_norm:
            continue
        if rule["row_chrom"] != chrom:
            continue

        if "x_min" in rule and x_cn < rule["x_min"]:
            continue
        if "x_max" in rule and x_cn > rule["x_max"]:
            continue
        if "y_min" in rule and y_cn < rule["y_min"]:
            continue
        if "y_max" in rule and y_cn > rule["y_max"]:
            continue

        return {
            "call": rule["call"],
            "pattern_name": rule["pattern_name"],
            "score": rule["score"],
        }

    return None