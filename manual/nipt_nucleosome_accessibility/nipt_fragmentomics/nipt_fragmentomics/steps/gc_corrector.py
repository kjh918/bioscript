"""
gc_corrector.py
===============
bins_raw.parquet → bins_corrected.parquet

보정 전략
---------
[Count]
  CPM = count / library_total * 1e6
  library_total: 전체 autosome quality-pass bin count 합산 (전체 라이브러리 기준)
  LOWESS GC fit → correction term → log2_corrected

WPS 는 wps_compute.py 에서 bp 단위로 별도 계산합니다.
이 파일에서는 WPS 를 처리하지 않습니다.

출력 추가 컬럼
--------------
  mappability_pass
  log2_corrected_short / long
  gc_correction_short / long
  log2_ratio
"""

from __future__ import annotations

import logging
import os

import numpy as np
import pandas as pd
from statsmodels.nonparametric.smoothers_lowess import lowess

from nipt_fragmentomics.core.constants import (
    MIN_MAPPABILITY, GC_RANGE, GC_LOWESS_FRAC, GC_PSEUDOCOUNT, GC_CLIP,
)

log = logging.getLogger(__name__)
SEX_CHROMS = {"chrX", "chrY"}


# ─────────────────────────────────────────────────────────────────────
# Count LOWESS GC 보정
# ─────────────────────────────────────────────────────────────────────
def _correct_count(
    counts:    np.ndarray,
    gc:        np.ndarray,
    fit_mask:  np.ndarray,
    lib_total: float,
    pseudocount: float = GC_PSEUDOCOUNT,
    frac:        float = GC_LOWESS_FRAC,
    clip:        float = GC_CLIP,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    count → log2 CPM → LOWESS GC 보정.
    lib_total: 전체 autosome quality-pass bin count 합산.

    Returns: log2_corrected, gc_fit, correction
    """
    v   = counts.astype(np.float64)
    cpm = v / max(lib_total, 1.0) * 1_000_000
    y   = np.log2(cpm + pseudocount)

    active = fit_mask & np.isfinite(gc) & np.isfinite(y)
    if active.sum() < 50:
        log.warning("GC fit 가능 bin 부족 (%d) — 보정 생략", active.sum())
        return (y.astype(np.float32),
                np.full(len(y), np.nan, np.float32),
                np.zeros(len(y), np.float32))

    fit_vals = lowess(y[active], gc[active], frac=frac, return_sorted=False)
    order    = np.argsort(gc[active])
    fit_all  = np.interp(gc,
                         gc[active][order], fit_vals[order],
                         left=fit_vals[order][0], right=fit_vals[order][-1])

    baseline   = np.nanmedian(fit_vals)
    correction = np.clip(fit_all - baseline, -clip, clip).astype(np.float32)
    corrected  = (y - correction).astype(np.float32)

    return corrected, fit_all.astype(np.float32), correction


# ─────────────────────────────────────────────────────────────────────
# 공개 인터페이스
# ─────────────────────────────────────────────────────────────────────
def run(
    raw_path:        str,
    out_path:        str,
    min_mappability: float = MIN_MAPPABILITY,
    gc_frac:         float = GC_LOWESS_FRAC,
) -> pd.DataFrame:
    """bins_raw.parquet → bins_corrected.parquet"""
    df = pd.read_parquet(raw_path)
    log.info("bins_raw 로드: %d rows", len(df))

    gc = df["gc"].values.astype(np.float64)
    mp = df["mappability"].fillna(np.nan).values.astype(np.float64)

    # ── 1. Quality mask ──────────────────────────────────────────
    is_sex       = df["chrom"].isin(SEX_CHROMS).values
    gc_zero_mask = ~np.isfinite(gc) | (gc <= 0)

    mp_pass = (is_sex | np.isnan(mp) | (mp >= min_mappability)) & ~gc_zero_mask
    df["mappability_pass"] = mp_pass
    log.info("Quality 마스킹: %d / %d bins 제외", (~mp_pass).sum(), len(df))

    gc_valid      = (gc >= GC_RANGE[0]) & (gc <= GC_RANGE[1])
    fit_mask_base = mp_pass & gc_valid & ~is_sex

    # ── 2. Count GC 보정 (short / long 각각) ─────────────────────
    for frag in ("short", "long"):
        cnt      = df[f"{frag}_count"].values.astype(np.float64)
        cnt_m    = np.where(mp_pass, cnt, 0.0)

        # library total: 전체 autosome quality-pass bin 합산
        lib_total = float(np.nansum(cnt_m[mp_pass & ~is_sex]))
        log.info("[%s] library total = %.0f reads", frag, lib_total)

        fit_mask = fit_mask_base & (cnt >= 1)

        corrected, gc_fit, correction = _correct_count(
            cnt_m, gc, fit_mask,
            lib_total=lib_total,
            pseudocount=GC_PSEUDOCOUNT,
            frac=gc_frac,
            clip=GC_CLIP,
        )
        df[f"log2_corrected_{frag}"] = corrected
        df[f"gc_correction_{frag}"]  = correction
        df[f"gc_fit_{frag}"]         = gc_fit

    # ── 3. log2 ratio ─────────────────────────────────────────────
    df["log2_ratio"] = (
        df["log2_corrected_short"].values - df["log2_corrected_long"].values
    ).astype(np.float32)

    # ── 4. Quality 실패 bin → NaN ─────────────────────────────────
    bad = ~df["mappability_pass"]
    for col in ["log2_corrected_short", "log2_corrected_long", "log2_ratio"]:
        if col in df.columns:
            df.loc[bad, col] = np.nan

    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    df.to_parquet(out_path, engine="pyarrow", index=False)
    log.info("bins_corrected 저장: %s  (quality_pass=%d/%d)",
             out_path, int(mp_pass.sum()), len(df))
    return df