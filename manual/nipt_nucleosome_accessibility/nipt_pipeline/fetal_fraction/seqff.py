"""
fetal_fraction/seqff.py
=======================
SeqFF: autosome bin의 GC-corrected short/total ratio → 선형모델 FF 추정

공식 (Kim et al. 2015 / Papageorgiou 2011 기반)
------------------------------------------------
short_fraction_i = short_gc_i / total_gc_i   (per bin)
FF (%) = intercept + slope × mean(short_fraction of reference bins)

기본 계수 (내장)
----------------
intercept : -72.6570
slope     :  97.5714
출처: SeqFF 원저자 공개 계수 (cfDNA short fragment 비율 기반)

- 입력: run_normalize 출력 normalized.tsv
  컬럼: chrom, bin_start, bin_end, gc, mappability,
         total_raw, short_raw, long_raw,
         total_gc, short_gc, long_gc,
         zscore_total, zscore_short, zscore_long
- reference bin: chr15, chr18, chr21, chrX, chrY 제외 autosome
- total_gc == 0 bin 제외
"""
from __future__ import annotations

import logging
from typing import Optional

import numpy as np

from cnv.utils import read_tsv
from .utils import FetalFractionResult

log = logging.getLogger(__name__)

# 기본 내장 계수 (SeqFF, Kim et al. 2015)
_DEFAULT_INTERCEPT: float = -72.6570
_DEFAULT_SLOPE:     float =  97.5714

# reference bin 제외 염색체
_EXCLUDED = {"chr15", "chr18", "chr21", "chrX", "chrY"}


def estimate(
    normalized_tsv: str,
    sample_id:      str   = "",
    intercept:      float = _DEFAULT_INTERCEPT,
    slope:          float = _DEFAULT_SLOPE,
) -> FetalFractionResult:
    """
    Parameters
    ----------
    normalized_tsv : run_normalize 출력 TSV
    sample_id      : 결과 메타 표기용
    intercept      : 선형모델 절편 (기본값: SeqFF 원저자)
    slope          : 선형모델 기울기 (기본값: SeqFF 원저자)

    Returns
    -------
    FetalFractionResult (ff_chry 관련 필드는 None)
    """
    rows = read_tsv(normalized_tsv)
    if not rows:
        raise ValueError(f"빈 TSV: {normalized_tsv}")

    short_fractions: list[float] = []

    for r in rows:
        chrom = r["chrom"]
        if chrom in _EXCLUDED:
            continue

        total_gc = float(r.get("total_gc", 0) or 0)
        short_gc = float(r.get("short_gc", 0) or 0)

        # total_gc == 0 bin 제외 (mappability mask 이미 적용된 bin)
        if total_gc <= 0:
            continue

        short_fractions.append(short_gc / total_gc)

    if not short_fractions:
        log.warning("유효한 reference bin 없음 — SeqFF 계산 불가")
        return FetalFractionResult(
            sample_id           = sample_id,
            ff_chry             = None,
            chry_mean_cov       = None,
            autosome_mean_cov   = None,
            chry_bin_count      = None,
            ff_seqff            = None,
            short_fraction_mean = None,
            seqff_intercept     = intercept,
            seqff_slope         = slope,
            method_note         = "유효 bin 없음",
        )

    sf_mean = float(np.mean(short_fractions))
    ff      = intercept + slope * sf_mean

    # 물리적으로 불가능한 범위 클램프 (0~100%)
    ff_clamped = float(np.clip(ff, 0.0, 100.0))
    if ff != ff_clamped:
        log.warning("SeqFF 범위 초과 (raw=%.4f) → clamp → %.4f", ff, ff_clamped)

    log.info(
        "SeqFF FF: sample=%s ff=%.2f%% sf_mean=%.6f bins=%d",
        sample_id, ff_clamped, sf_mean, len(short_fractions),
    )

    return FetalFractionResult(
        sample_id           = sample_id,
        ff_chry             = None,
        chry_mean_cov       = None,
        autosome_mean_cov   = None,
        chry_bin_count      = None,
        ff_seqff            = ff_clamped,
        short_fraction_mean = sf_mean,
        seqff_intercept     = intercept,
        seqff_slope         = slope,
        method_note         = "SeqFF linear model",
    )