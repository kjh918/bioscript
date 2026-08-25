"""
fetal_fraction/chry.py
======================
chrY coverage 기반 FF 추정 (male fetus만 유효)
+ Fragment size ratio 기반 FF proxy (남녀 공용)

공식
----
chrY 기반:
  ff_chry (%) = (chrY_mean_cov / autosome_mean_cov) × 2 × 100

Fragment size ratio (Illumina 특허 US11217330 기반):
  SR = P(100-150bp) / P(163-169bp)
  → FF와 Pearson r=0.831 (성별 무관)
  → 절대 FF% 추정 불가 (cohort 없음)
  → 샘플 간 상대 비교용 proxy로만 사용

chrY blacklist (hg38)
---------------------
1. PAR1 (chrY:10,000-2,781,479)  — chrX와 서열 동일, reads 누출
2. PAR2 (chrY:56,887,902-57,217,415)
3. Ampliconic DYZ1/DYZ2 (chrY:10,000,000-13,000,000)
4. Centromere 주변 (위 3번에 포함)
"""
from __future__ import annotations

import logging
from typing import Optional

import numpy as np

from cnv.utils import read_tsv
from .utils import FetalFractionResult

log = logging.getLogger(__name__)

_EXCLUDED     = {"chr15", "chr18", "chr21", "chrX", "chrY"}
_AUTOSOME_REF = {f"chr{i}" for i in range(1, 23)} - _EXCLUDED

# chrY blacklist (hg38, half-open intervals)
_CHRY_BLACKLIST: list[tuple[int, int]] = [
    (0,          3_000_000),   # PAR1 + buffer
    (10_000_000, 13_000_000),  # ampliconic DYZ1/DYZ2 + centromere
    (56_000_000, 57_300_000),  # PAR2 + buffer
]

# Fragment size ratio 계산용 구간 (bp)
_SHORT_LO  = 100
_SHORT_HI  = 150   # short fetal-enriched 구간
_LONG_LO   = 163
_LONG_HI   = 169   # nucleosome-protected peak (maternal dominant)


def _in_blacklist(start: int) -> bool:
    return any(lo <= start < hi for lo, hi in _CHRY_BLACKLIST)


def _compute_size_ratio(bincount_tsv: str) -> Optional[float]:
    """
    Fragment size ratio = P(100-150bp) / P(163-169bp)

    bincount.tsv의 short_raw / long_raw 컬럼으로 근사 계산.

    실제 구현 한계:
      binning.py에서 short = <150bp, long = ≥150bp 로만 나눠져 있어
      163-169bp 구간을 직접 분리할 수 없음.
      → 근사: short_raw / (total_raw - short_raw) 사용
        = P(<150bp) / P(≥150bp)
        이 비율도 FF와 양의 상관관계를 가짐 (Illumina 특허 근거)

    정확한 구현을 위해서는 binning.py에서 TLEN 분포를 직접 계산해야 함.
    현재는 proxy값으로만 사용.
    """
    rows = read_tsv(bincount_tsv)
    if not rows:
        return None

    total_short = 0.0
    total_long  = 0.0

    for r in rows:
        chrom = r.get("chrom", "")
        # autosome ref bin만 사용 (성별 영향 최소화)
        if chrom not in _AUTOSOME_REF:
            continue
        short = float(r.get("short", 0) or 0)
        total = float(r.get("total", 0) or 0)
        long  = total - short
        total_short += short
        total_long  += max(long, 0.0)

    total_reads = total_short + total_long
    if total_reads <= 0 or total_long <= 0:
        return None

    # SR = short_fraction / long_fraction (proxy)
    sr = (total_short / total_reads) / (total_long / total_reads)
    return float(sr)


def estimate(
    bincount_tsv: str,
    sample_id:    str = "",
) -> FetalFractionResult:
    """
    chrY coverage 기반 FF + fragment size ratio proxy 동시 계산.

    Returns
    -------
    FetalFractionResult
      ff_chry                : male fetus → FF (%), female → None
      ff_short_fragment_ratio: 남녀 공용 FF proxy (절대값 아님, 비교용)
    """
    rows = read_tsv(bincount_tsv)
    if not rows:
        raise ValueError(f"빈 TSV: {bincount_tsv}")

    chry_covs:          list[float] = []
    chry_covs_excluded: list[float] = []
    autosome_covs:      list[float] = []

    for r in rows:
        chrom = r["chrom"]
        total = float(r["total"])
        start = int(r["bin_start"])

        if chrom == "chrY":
            if _in_blacklist(start):
                chry_covs_excluded.append(total)
                log.debug("chrY blacklist 제외: %d (cov=%.1f)", start, total)
            else:
                chry_covs.append(total)
        elif chrom in _AUTOSOME_REF:
            autosome_covs.append(total)

    log.info("chrY bins: 유효=%d 제외(blacklist)=%d",
             len(chry_covs), len(chry_covs_excluded))

    # ── Fragment size ratio (남녀 공용) ───────────────────────────────
    sr = _compute_size_ratio(bincount_tsv)
    if sr is not None:
        log.info("fragment size ratio (SR): %.4f (FF proxy, 비교용)", sr)

    # ── chrY FF ───────────────────────────────────────────────────────
    if not chry_covs:
        log.warning("유효한 chrY bin 없음 — chrY FF 계산 불가")
        return FetalFractionResult(
            sample_id                = sample_id,
            ff_chry                  = None,
            chry_mean_cov            = None,
            autosome_mean_cov        = None,
            chry_bin_count           = 0,
            ff_short_fragment_ratio  = sr,
            method_note              = "chrY 유효 bin 없음 (blacklist 제외 후)",
        )

    if not autosome_covs:
        log.warning("autosome bin 없음 — chrY FF 계산 불가")
        return FetalFractionResult(
            sample_id                = sample_id,
            ff_chry                  = None,
            chry_mean_cov            = float(np.mean(chry_covs)),
            autosome_mean_cov        = None,
            chry_bin_count           = len(chry_covs),
            ff_short_fragment_ratio  = sr,
            method_note              = "autosome bin 없음",
        )

    chry_mean     = float(np.mean(chry_covs))
    autosome_mean = float(np.mean(autosome_covs))

    if autosome_mean <= 0:
        log.warning("autosome mean coverage = 0 — chrY FF 계산 불가")
        return FetalFractionResult(
            sample_id                = sample_id,
            ff_chry                  = None,
            chry_mean_cov            = chry_mean,
            autosome_mean_cov        = autosome_mean,
            chry_bin_count           = len(chry_covs),
            ff_short_fragment_ratio  = sr,
            method_note              = "autosome mean = 0",
        )

    raw_ratio = chry_mean / autosome_mean
    ff        = raw_ratio * 2.0 * 100.0

    # female fetus 판별: raw_ratio < 0.01
    # 정상 male FF=10% → raw_ratio ≈ 0.05
    # female / background noise → raw_ratio < 0.01
    is_female = raw_ratio < 0.01

    if is_female:
        ff_result = None
        note = (
            f"female fetus 추정 (chrY/autosome={raw_ratio:.4f} < 0.01)"
            f" — SR proxy={sr:.4f}" if sr else " — SR 계산 불가"
        )
        log.info("chrY FF: %s → female fetus 추정 (ratio=%.4f, SR=%.4f)",
                 sample_id, raw_ratio, sr or 0.0)
    else:
        ff_result = float(np.clip(ff, 0.0, 100.0))
        note = (
            f"chrY-based (blacklist 적용, bins={len(chry_covs)})"
            f", SR proxy={sr:.4f}" if sr else ""
        )
        log.info("chrY FF: %s ff=%.2f%% chrY=%.2f auto=%.2f bins=%d SR=%.4f",
                 sample_id, ff_result, chry_mean, autosome_mean,
                 len(chry_covs), sr or 0.0)

    return FetalFractionResult(
        sample_id                = sample_id,
        ff_chry                  = ff_result,
        chry_mean_cov            = chry_mean,
        autosome_mean_cov        = autosome_mean,
        chry_bin_count           = len(chry_covs),
        ff_short_fragments_ratio = sr,
        method_note              = note,
    )