"""
fetal_fraction/chry.py
======================
chrY 평균 bin coverage / autosome 평균 bin coverage 기반 FF 추정

공식
----
ff_chry (%) = (chrY_mean_cov / autosome_mean_cov) * 2 * 100

- autosome: chr1~chr22 (chrX, chrY, chr15, chr18, chr21 제외)
- 입력: run_binning 출력 bincount.tsv
  컬럼: chrom, bin_start, bin_end, gc, mappability, total, short, long
- raw Y-ratio 그대로 반환 (sex 판정 없음)
"""
from __future__ import annotations

import logging
from typing import Optional

import numpy as np

from cnv.utils import read_tsv
from .utils import FetalFractionResult

log = logging.getLogger(__name__)

# autosome reference: aneuploidy 빈발 염색체 + sex chrom 제외
_EXCLUDED = {"chr15", "chr18", "chr21", "chrX", "chrY"}
_AUTOSOME_REF = {f"chr{i}" for i in range(1, 23)} - _EXCLUDED


def estimate(
    bincount_tsv: str,
    sample_id:    str = "",
) -> FetalFractionResult:
    """
    Parameters
    ----------
    bincount_tsv : run_binning 출력 TSV
    sample_id    : 결과 메타 표기용

    Returns
    -------
    FetalFractionResult (ff_seqff 필드는 None)
    """
    rows = read_tsv(bincount_tsv)
    if not rows:
        raise ValueError(f"빈 TSV: {bincount_tsv}")

    chry_covs:     list[float] = []
    autosome_covs: list[float] = []

    for r in rows:
        chrom = r["chrom"]
        total = float(r["total"])

        if chrom == "chrY":
            chry_covs.append(total)
        elif chrom in _AUTOSOME_REF:
            autosome_covs.append(total)

    if not chry_covs:
        log.warning("chrY bin 없음 — FF 계산 불가")
        return FetalFractionResult(
            sample_id         = sample_id,
            ff_chry           = None,
            chry_mean_cov     = None,
            autosome_mean_cov = None,
            chry_bin_count    = 0,
            ff_seqff          = None,
            short_fraction_mean = None,
            method_note       = "chrY bin 없음",
        )

    if not autosome_covs:
        log.warning("autosome bin 없음 — FF 계산 불가")
        return FetalFractionResult(
            sample_id         = sample_id,
            ff_chry           = None,
            chry_mean_cov     = float(np.mean(chry_covs)),
            autosome_mean_cov = None,
            chry_bin_count    = len(chry_covs),
            ff_seqff          = None,
            short_fraction_mean = None,
            method_note       = "autosome bin 없음",
        )

    chry_mean     = float(np.mean(chry_covs))
    autosome_mean = float(np.mean(autosome_covs))

    if autosome_mean == 0:
        log.warning("autosome 평균 coverage = 0 — FF 계산 불가")
        ff = None
        note = "autosome mean coverage = 0"
    else:
        ff   = (chry_mean / autosome_mean) * 2.0 * 100.0
        note = "chrY-based raw ratio"

    log.info(
        "chrY FF: sample=%s ff=%.2f%% chrY_mean=%.2f autosome_mean=%.2f bins=%d",
        sample_id, ff or 0.0, chry_mean, autosome_mean, len(chry_covs),
    )

    return FetalFractionResult(
        sample_id         = sample_id,
        ff_chry           = ff,
        chry_mean_cov     = chry_mean,
        autosome_mean_cov = autosome_mean,
        chry_bin_count    = len(chry_covs),
        ff_seqff          = None,
        short_fraction_mean = None,
        method_note       = note,
    )