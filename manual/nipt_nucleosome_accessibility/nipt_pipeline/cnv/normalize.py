"""
cnv/normalize.py

bin count TSV → normalized Z-score TSV

단계:
  1. LOESS GC correction  (total / short / long 각각)
  2. genome-wide Z-score  (reference chrom 기반)

reference chrom: chr15, chr18, chr21, chrX, chrY 제외한 autosome bin

출력 컬럼:
  chrom, bin_start, bin_end, gc, mappability,
  total_raw, short_raw, long_raw,
  total_gc, short_gc, long_gc,
  zscore_total, zscore_short, zscore_long
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess

from .utils import read_tsv, write_tsv

logger = logging.getLogger(__name__)

DEFAULT_EXCLUDED = {"chr15", "chr18", "chr21", "chrX", "chrY"}

FIELDNAMES = [
    "chrom", "bin_start", "bin_end", "gc", "mappability",
    "total_raw", "short_raw", "long_raw",
    "total_gc",  "short_gc",  "long_gc",
    "zscore_total", "zscore_short", "zscore_long",
]


# ---------------------------------------------------------------------------
# GC correction
# ---------------------------------------------------------------------------

def _loess_gc_correct(
    counts: np.ndarray,
    gc:     np.ndarray,
    frac:   float = 0.3,
    max_correction_ratio: float = 3.0,
    gc_trim_pct:          float = 0.01,
) -> np.ndarray:
    """
    LOESS GC bias 보정.
    corrected = count / predicted_by_gc * median(predicted)

    안전장치
    --------
    1. GC 극단 bin (상하 gc_trim_pct) 제외 후 LOESS 적합
    2. predicted 하한: median_pred * 0.1 클램프
    3. 보정 후 ratio > max_correction_ratio 인 bin → raw fallback
    """
    counts = counts.astype(np.float64)
    valid  = counts > 0

    if valid.sum() < 10:
        logger.warning("유효 bin이 너무 적어 GC correction 스킵")
        return counts.copy()

    # 1. GC 극단 제외 후 LOESS 적합용 mask
    gc_lo    = np.percentile(gc[valid], gc_trim_pct * 100)
    gc_hi    = np.percentile(gc[valid], (1 - gc_trim_pct) * 100)
    fit_mask = valid & (gc >= gc_lo) & (gc <= gc_hi)

    if fit_mask.sum() < 10:
        logger.warning("GC trim 후 적합 bin 부족 — GC correction 스킵")
        return counts.copy()

    # LOESS 적합 (fit_mask bin만, GC 오름차순 정렬)
    sorted_idx = np.argsort(gc[fit_mask])
    smoothed   = lowess(
        counts[fit_mask][sorted_idx],
        gc[fit_mask][sorted_idx],
        frac          = frac,
        return_sorted = False,
        it            = 3,
    )

    # LOESS 곡선 음수 클램프
    # GC 극저값 구간에서 smoothing 곡선이 음수로 내려가는 경우 차단
    smoothed = np.maximum(smoothed, 0.0)

    # 2. 전체 GC 범위 보간 (fit 범위 밖은 경계값 clamp — extrapolation 없음)
    predicted_all = np.interp(
        gc,
        gc[fit_mask][sorted_idx],
        smoothed,
    )

    median_pred = float(np.median(predicted_all[fit_mask]))
    if median_pred <= 0:
        logger.warning("median_pred <= 0 — GC correction 스킵")
        return counts.copy()

    # predicted 하한: median_pred * 0.1 (0 나누기 방지)
    predicted_safe = np.maximum(predicted_all, median_pred * 0.1)

    # 3. 보정값 계산 (np.where 대신 직접 배열 할당 — 후속 인덱싱 보장)
    corrected          = counts.copy()
    corrected[valid]   = counts[valid] / predicted_safe[valid] * median_pred

    # 4. 보정 배율 이상치 → raw fallback
    # np.where 반환값은 새 배열이라 인덱스 할당이 원본에 반영되지 않으므로
    # 반드시 corrected 배열에 직접 인덱싱으로 fallback
    ratio        = np.ones(len(counts), dtype=np.float64)
    safe_denom   = np.where(valid & (counts > 0), counts, 1.0)
    ratio[valid] = corrected[valid] / safe_denom[valid]

    outlier_mask = valid & (
        (ratio > max_correction_ratio) | (ratio < 1.0 / max_correction_ratio)
    )
    n_outlier = int(outlier_mask.sum())
    if n_outlier > 0:
        logger.warning(
            "GC correction 이상치 %d bin (ratio %.2f~%.2f 범위 초과) → raw fallback",
            n_outlier, 1.0 / max_correction_ratio, max_correction_ratio,
        )
        corrected[outlier_mask] = counts[outlier_mask]  # 직접 인덱싱으로 fallback

    return corrected


# ---------------------------------------------------------------------------
# Z-score
# ---------------------------------------------------------------------------

def _zscore(
    values:   np.ndarray,
    ref_mask: np.ndarray,
) -> np.ndarray:
    """
    ref_mask bin의 median / MAD 기반 robust Z-score 계산.

    mean/std 대신 median/MAD 사용 이유
    ------------------------------------
    - GC correction 후 잔여 이상치 bin이 mean/std를 오염시켜
      정상 bin의 Z-score 범위가 예상과 달라지는 문제 방지
    - MAD(Median Absolute Deviation) 기반 scale:
        z = (x - median) / (MAD * 1.4826)
      1.4826 = 정규분포 가정 하 MAD → std 환산 계수
    - MAD=0 (극단적으로 균일한 분포) 이면 std fallback
    """
    ref_vals = values[ref_mask]
    if len(ref_vals) < 2:
        logger.warning("reference bin이 부족해 Z-score 계산 불가")
        return np.full_like(values, np.nan, dtype=float)

    median = float(np.median(ref_vals))
    mad    = float(np.median(np.abs(ref_vals - median)))
    scale  = mad * 1.4826  # MAD → 정규분포 std 환산

    if scale == 0:
        # fallback: 일반 std
        scale = float(np.std(ref_vals, ddof=1))
        logger.warning("MAD=0 — std fallback (scale=%.6f)", scale)

    if scale == 0:
        logger.warning("scale=0 — Z-score 계산 불가")
        return np.full_like(values, np.nan, dtype=float)

    return (values - median) / scale


# ---------------------------------------------------------------------------
# 메인
# ---------------------------------------------------------------------------

def run_normalize(
    bin_tsv_path: str,
    output_path: str,
    excluded_chroms: set[str] | None = None,
    gc_correction: bool = True,
    gc_loess_span: float = 0.3,
) -> None:
    """
    Parameters
    ----------
    bin_tsv_path    : run_binning 출력 TSV
    output_path     : 정규화 결과 TSV
    excluded_chroms : Z-score reference에서 제외할 염색체 집합
    gc_correction   : LOESS GC correction 수행 여부
    gc_loess_span   : LOESS frac 파라미터
    """
    excluded = excluded_chroms or DEFAULT_EXCLUDED

    rows = read_tsv(bin_tsv_path)
    if not rows:
        raise ValueError(f"빈 TSV: {bin_tsv_path}")

    # numpy array 변환
    chroms     = np.array([r["chrom"]     for r in rows])
    bin_starts = np.array([int(r["bin_start"]) for r in rows])
    bin_ends   = np.array([int(r["bin_end"])   for r in rows])
    gc         = np.array([float(r["gc"])       for r in rows])
    mappability = np.array([float(r["mappability"]) for r in rows])
    total_raw  = np.array([int(r["total"]) for r in rows], dtype=float)
    short_raw  = np.array([int(r["short"]) for r in rows], dtype=float)
    long_raw   = np.array([int(r["long"])  for r in rows], dtype=float)

    # GC correction
    if gc_correction:
        logger.info("LOESS GC correction 수행 (span=%.2f)", gc_loess_span)
        total_gc = _loess_gc_correct(total_raw, gc, gc_loess_span)
        short_gc = _loess_gc_correct(short_raw, gc, gc_loess_span)
        long_gc  = _loess_gc_correct(long_raw,  gc, gc_loess_span)
    else:
        total_gc = total_raw.copy()
        short_gc = short_raw.copy()
        long_gc  = long_raw.copy()

    # chr19 GC 보정 후 임의 가중치 적용
    # chr19는 GC content ~48%로 LOESS correction 후 과도하게 낮아지는 경향
    # 30% 상향 가중치로 보정
    chr19_mask = chroms == "chr19"
    if chr19_mask.any():
        CHR19_WEIGHT = 1.30
        total_gc[chr19_mask] *= CHR19_WEIGHT
        short_gc[chr19_mask] *= CHR19_WEIGHT
        long_gc[chr19_mask]  *= CHR19_WEIGHT
        logger.info(
            "chr19 GC 가중치 적용 (×%.2f, bins=%d)",
            CHR19_WEIGHT, int(chr19_mask.sum()),
        )

    # reference bin mask: 제외 염색체 외의 bin
    ref_mask = np.array([c not in excluded for c in chroms])
    logger.info(
        "Z-score reference bin 수: %d / %d (제외 chroms: %s)",
        ref_mask.sum(), len(chroms), sorted(excluded),
    )

    # genome-wide Z-score
    zscore_total = _zscore(total_gc, ref_mask)
    zscore_short = _zscore(short_gc, ref_mask)
    zscore_long  = _zscore(long_gc,  ref_mask)

    # Z-score 이상치 clamp: GC correction 잔여 이상치 bin 제거
    # MAD 기반 Z-score 이므로 ±10 초과는 실질적으로 의미 없는 이상치
    _ZSCORE_CLAMP = 10.0
    for arr in (zscore_total, zscore_short, zscore_long):
        np.clip(arr, -_ZSCORE_CLAMP, _ZSCORE_CLAMP, out=arr)

    # 결과 조립
    out_rows = []
    for i in range(len(rows)):
        out_rows.append({
            "chrom":        chroms[i],
            "bin_start":    bin_starts[i],
            "bin_end":      bin_ends[i],
            "gc":           f"{gc[i]:.4f}",
            "mappability":  f"{mappability[i]:.4f}",
            "total_raw":    int(total_raw[i]),
            "short_raw":    int(short_raw[i]),
            "long_raw":     int(long_raw[i]),
            "total_gc":     f"{total_gc[i]:.4f}",
            "short_gc":     f"{short_gc[i]:.4f}",
            "long_gc":      f"{long_gc[i]:.4f}",
            "zscore_total": f"{zscore_total[i]:.4f}" if not np.isnan(zscore_total[i]) else "NA",
            "zscore_short": f"{zscore_short[i]:.4f}" if not np.isnan(zscore_short[i]) else "NA",
            "zscore_long":  f"{zscore_long[i]:.4f}"  if not np.isnan(zscore_long[i])  else "NA",
        })

    write_tsv(out_rows, output_path, FIELDNAMES)
    logger.info("normalization 완료 → %s", output_path)