"""
wps_normalizer.py
=================
bins_wps_raw.parquet (1kb bin) → bins_wps_norm.parquet

Genome-wide WPS normalization.

정규화 공식 (Snyder et al. 2016 방식)
--------------------------------------
  adjusted_WPS = (WPS - local_median) / depth * 100

  WPS       : bin 내 spanning - endpoints 배열의 median (raw)
  local_median : 인접 bins 의 WPS median (1Mb 슬라이딩 윈도우, ±500 bins)
  depth     : bin 내 short 또는 long coverage median (depth proxy)
              = (short_count + long_count) / bin_len  [reads per bp]

  depth = 0 인 bin → NaN

추가 보정
---------
  GC bias: bins_corrected 의 gc_correction 값 재사용 가능하지만,
           WPS 는 count 와 GC 의존성이 달라 별도 적용하지 않음.
  mappability: mappability_pass=False bin → NaN

출력 컬럼
---------
  기존 bins_wps_raw 컬럼 모두 유지 +
  short_wps_L_norm, long_wps_L_norm,
  short_wps_S_norm, long_wps_S_norm
"""

from __future__ import annotations

import logging
import os

import numpy as np
import pandas as pd

from nipt_fragmentomics.core.constants import MIN_MAPPABILITY

log = logging.getLogger(__name__)

SEX_CHROMS = {"chrX", "chrY"}
LOCAL_WIN   = 500   # ±500 bins = 1kb bin 기준 ±500kb


def _local_median(arr: np.ndarray, half_win: int = LOCAL_WIN) -> np.ndarray:
    """
    각 bin 의 ±half_win 범위 median (NaN-safe).
    전체 배열에서 슬라이딩 윈도우 median 계산.
    """
    n   = len(arr)
    out = np.full(n, np.nan, dtype=np.float32)
    for i in range(n):
        lo  = max(0, i - half_win)
        hi  = min(n, i + half_win + 1)
        seg = arr[lo:hi]
        fin = seg[np.isfinite(seg)]
        if len(fin) > 0:
            out[i] = float(np.median(fin))
    return out


def _normalize_wps(
    wps:      np.ndarray,   # raw WPS (bin median)
    depth:    np.ndarray,   # reads per bp proxy
    half_win: int = LOCAL_WIN,
) -> np.ndarray:
    """
    adjusted_WPS = (WPS - local_median) / depth * 100

    depth = 0 → NaN
    """
    local_med = _local_median(wps, half_win)
    adjusted  = wps - local_med

    out = np.full(len(wps), np.nan, dtype=np.float32)
    mask = np.isfinite(depth) & (depth > 0) & np.isfinite(adjusted)
    out[mask] = (adjusted[mask] / depth[mask] * 100.0).astype(np.float32)

    return out


def run(
    raw_path:    str,
    out_path:    str,
    bin_len:     int   = 1000,
    half_win:    int   = LOCAL_WIN,
    min_mappability: float = MIN_MAPPABILITY,
) -> pd.DataFrame:
    """
    bins_wps_raw.parquet → bins_wps_norm.parquet

    Parameters
    ----------
    raw_path    : bin_extractor 로 생성한 1kb WPS parquet
    bin_len     : bin 크기 (bp), depth 계산에 사용
    half_win    : local median 윈도우 (bins)
    """
    df = pd.read_parquet(raw_path)
    log.info("bins_wps_raw 로드: %d bins (bin_len=%d bp)", len(df), bin_len)

    # mappability 마스킹
    bad = pd.Series(False, index=df.index)
    if "mappability_pass" in df.columns:
        bad |= ~df["mappability_pass"].astype(bool)
    elif "mappability" in df.columns:
        bad |= df["mappability"].fillna(0) < min_mappability
    if "gc" in df.columns:
        bad |= ~df["gc"].between(0.01, 0.99, inclusive="both").fillna(True)

    log.info("불량 bin: %d / %d", bad.sum(), len(df))

    # depth proxy: reads per bp
    # short_count + long_count / bin_len
    for frag in ("short", "long"):
        cnt_col = f"{frag}_count"
        if cnt_col not in df.columns:
            log.warning("%s 컬럼 없음", cnt_col)
            continue

        counts = df[cnt_col].values.astype(np.float64)
        depth  = counts / bin_len
        depth[bad.values] = np.nan

        for mode in ("L", "S"):
            wps_col = f"{frag}_wps_{mode}"
            if wps_col not in df.columns:
                continue

            # 염색체별 독립 normalization
            norm_vals = np.full(len(df), np.nan, dtype=np.float32)

            for chrom, cdf in df.groupby("chrom"):
                idx  = cdf.index
                wps  = df.loc[idx, wps_col].values.astype(np.float32)
                dep  = depth[idx.values if hasattr(idx, 'values') else idx]

                # 불량 bin → NaN 처리
                wps[bad.loc[idx].values] = np.nan

                norm = _normalize_wps(wps, dep, half_win=half_win)
                norm_vals[idx] = norm

            out_col = f"{frag}_wps_{mode}_norm"
            df[out_col] = norm_vals.astype(np.float32)
            log.info("  %s → %s  (finite: %d / %d)",
                     wps_col, out_col,
                     int(np.isfinite(norm_vals).sum()), len(df))

    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    df.to_parquet(out_path, engine="pyarrow", index=False)
    log.info("bins_wps_norm 저장: %s", out_path)
    return df
