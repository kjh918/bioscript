"""
fetal_fraction.py
=================
두 가지 방법으로 Fetal Fraction (FF) 을 추정합니다.

방법 1 — SeqFF (short/long count ratio 기반)
  Larsen et al. 2017 선형 모델
  FF = α + β × (Σshort / Σlong)   [autosome mappability-pass bin 만 사용]

방법 2 — Y 염색체 (남아 샘플)
  chrY_mean_count / autosome_median_count 로 Y fraction 산출
  Y fraction > Y_MALE_THRESHOLD 이면 남아로 판정
  FF ≈ 2 × y_fraction   (Y 는 fetal 에만 존재 → haploid 보정)

출력 fetal_fraction.json
------------------------
{
  "seqff_ff":       0.12,
  "seqff_sl_ratio": 0.45,
  "y_ff":           0.11,
  "y_fraction":     0.055,
  "is_male_fetus":  true,
  "consensus_ff":   0.115,
  "method_used":    ["seqff", "y_chr"]
}
"""

from __future__ import annotations

import json
import logging
import os

import numpy as np
import pandas as pd

from nipt_fragmentomics.core.constants import (
    SEQFF_ALPHA, SEQFF_BETA, Y_MALE_THRESHOLD,
)

log = logging.getLogger(__name__)

SEX_CHROMS = {"chrX", "chrY"}


# ─────────────────────────────────────────────────────────────────────
# SeqFF
# ─────────────────────────────────────────────────────────────────────
def estimate_seqff(
    df: pd.DataFrame,
    alpha: float = SEQFF_ALPHA,
    beta:  float = SEQFF_BETA,
) -> tuple[float, float]:
    """
    autosome mappability-pass bin 의 short/long count 합산 비율로 FF 추정.

    Returns
    -------
    (ff_estimate, sl_ratio)  — 추정 불가 시 (nan, nan)
    """
    auto = df[~df["chrom"].isin(SEX_CHROMS) & df["mappability_pass"]]
    if auto.empty:
        return float("nan"), float("nan")

    short_sum = float(auto["short_count"].sum())
    long_sum  = float(auto["long_count"].sum())

    if long_sum == 0:
        return float("nan"), float("nan")

    sl_ratio = short_sum / long_sum
    ff       = float(np.clip(alpha + beta * sl_ratio, 0.0, 1.0))
    return ff, sl_ratio


# ─────────────────────────────────────────────────────────────────────
# Y 염색체 기반
# ─────────────────────────────────────────────────────────────────────
def estimate_y_chr(
    df: pd.DataFrame,
    threshold: float = Y_MALE_THRESHOLD,
) -> tuple[float, float, bool]:
    """
    chrY coverage 대비 autosome median 으로 Y fraction 계산.

    Returns
    -------
    (ff_estimate, y_fraction, is_male_fetus)
    """
    y_rows    = df[df["chrom"] == "chrY"]
    auto_rows = df[~df["chrom"].isin(SEX_CHROMS) & df["mappability_pass"]]

    if y_rows.empty or auto_rows.empty:
        return float("nan"), float("nan"), False

    y_mean       = float(y_rows["short_count"].mean())
    auto_median  = float(np.median(auto_rows["short_count"].values))

    if auto_median == 0:
        return float("nan"), float("nan"), False

    y_frac   = y_mean / auto_median
    is_male  = y_frac > threshold

    if not is_male:
        return float("nan"), float(y_frac), False

    # haploid Y → FF ≈ 2 × y_fraction
    ff = float(np.clip(2.0 * y_frac, 0.0, 1.0))
    return ff, float(y_frac), True


# ─────────────────────────────────────────────────────────────────────
# 공개 인터페이스
# ─────────────────────────────────────────────────────────────────────
def run(
    corrected_path: str,
    out_path:       str,
    alpha: float = SEQFF_ALPHA,
    beta:  float = SEQFF_BETA,
) -> dict:
    """
    bins_corrected.parquet → fetal_fraction.json

    Returns
    -------
    result dict (json 과 동일 내용)
    """
    df = pd.read_parquet(corrected_path)

    seqff_ff, sl_ratio      = estimate_seqff(df, alpha=alpha, beta=beta)
    y_ff, y_frac, is_male   = estimate_y_chr(df)

    log.info("SeqFF  FF=%.4f  sl_ratio=%.4f", seqff_ff, sl_ratio)
    if is_male:
        log.info("Y-chr  FF=%.4f  y_frac=%.6f  (남아)", y_ff, y_frac)
    else:
        log.info("Y-chr  is_male=False  y_frac=%.6f",
                 y_frac if np.isfinite(y_frac) else -1)

    estimates    = [v for v in [seqff_ff, y_ff] if np.isfinite(v)]
    consensus_ff = float(np.mean(estimates)) if estimates else float("nan")
    methods      = (["seqff"] if np.isfinite(seqff_ff) else []) + \
                   (["y_chr"]  if is_male else [])

    def _safe(v, digits=4):
        return round(float(v), digits) if np.isfinite(v) else None

    result = {
        "seqff_ff":       _safe(seqff_ff),
        "seqff_sl_ratio": _safe(sl_ratio),
        "y_ff":           _safe(y_ff),
        "y_fraction":     _safe(y_frac, 6),
        "is_male_fetus":  is_male,
        "consensus_ff":   _safe(consensus_ff),
        "method_used":    methods,
    }

    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    log.info("fetal_fraction 저장: %s  consensus FF=%s",
             out_path, result["consensus_ff"])
    return result
