"""
wps_normalizer.py
=================
bins_wps_raw.parquet → bins_wps_norm.parquet

Genome-wide WPS normalization.

정규화 방식
-----------
1. mappability 불량 bin → NaN
2. 염색체별 독립 z-score
   z = (WPS - median(WPS_chrom)) / MAD(WPS_chrom) * 1.4826
   - median 차감: 로컬 baseline 제거
   - MAD 나누기: 염색체 간 depth 차이 보정 (robust)
3. 전체 상하위 0.5% winsorization (centromere spike 제거)

z-score 가 0이면 해당 염색체 중앙값 수준,
양수면 뉴클레오솜 보호 증가, 음수면 감소를 의미합니다.

출력 추가 컬럼
--------------
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
WINSOR_Q   = 0.005   # 상하위 0.5% winsorization


def _zscore_chrom(arr: np.ndarray, quality_mask: np.ndarray) -> np.ndarray:
    """
    염색체 단위 robust z-score.
    quality_mask=True 인 bin 만 median/MAD 계산에 사용.
    """
    ref = arr[quality_mask & np.isfinite(arr)]
    if len(ref) < 5:
        return np.full(len(arr), np.nan, dtype=np.float32)

    med   = float(np.median(ref))
    mad   = float(np.median(np.abs(ref - med)))
    scale = mad * 1.4826 if mad > 1e-6 else float(ref.std() or 1.0)

    out = np.full(len(arr), np.nan, dtype=np.float32)
    finite = np.isfinite(arr)
    out[finite] = ((arr[finite] - med) / scale).astype(np.float32)
    return out


def run(
    raw_path:    str,
    out_path:    str,
    bin_len:     int   = 1000,
    min_mappability: float = MIN_MAPPABILITY,
) -> pd.DataFrame:
    """
    bins_wps_raw.parquet → bins_wps_norm.parquet

    WPS 컬럼: short_wps_L/S, long_wps_L/S → *_norm (z-score)
    """
    df = pd.read_parquet(raw_path)
    log.info("bins_wps_raw 로드: %d bins", len(df))

    wps_cols = [c for c in ["short_wps_L", "long_wps_L",
                              "short_wps_S", "long_wps_S"]
                if c in df.columns]

    if not wps_cols:
        log.error("WPS 컬럼 없음 (short_wps_L 등). "
                  "bin_extractor 가 WPS 를 계산하는 버전인지 확인하세요.")
        df.to_parquet(out_path, engine="pyarrow", index=False)
        return df

    # ── quality mask ─────────────────────────────────────────────
    gc = df["gc"].values.astype(np.float64) if "gc" in df.columns else np.ones(len(df))
    mp = df["mappability"].fillna(np.nan).values.astype(np.float64) \
         if "mappability" in df.columns else np.ones(len(df))
    is_sex = df["chrom"].isin(SEX_CHROMS).values

    # BED 필터 컬럼 우선
    if "is_filtered" in df.columns:
        bad = df["is_filtered"].astype(bool).values & ~is_sex
    elif "is_low_mappability" in df.columns:
        bad = df["is_low_mappability"].astype(bool).values & ~is_sex
    else:
        bad = (~np.isnan(mp) & (mp < min_mappability) & ~is_sex)

    bad |= (~np.isfinite(gc) | (gc <= 0))
    quality_mask = ~bad & ~is_sex

    log.info("quality_pass: %d / %d bins", quality_mask.sum(), len(df))

    # ── 염색체별 z-score ─────────────────────────────────────────
    for col in wps_cols:
        norm_vals = np.full(len(df), np.nan, dtype=np.float32)
        wps_all   = df[col].values.astype(np.float32)

        for chrom, cdf in df.groupby("chrom"):
            idx       = cdf.index
            wps_c     = wps_all[idx]
            qmask_c   = quality_mask[idx]
            # 불량 bin → NaN
            wps_c     = np.where(bad[idx], np.nan, wps_c)

            norm_c    = _zscore_chrom(wps_c, qmask_c)
            norm_vals[idx] = norm_c

        # winsorization (전체, 상하위 0.5%)
        finite = norm_vals[np.isfinite(norm_vals)]
        if len(finite) > 100:
            lo = float(np.quantile(finite, WINSOR_Q))
            hi = float(np.quantile(finite, 1.0 - WINSOR_Q))
            norm_vals = np.clip(norm_vals, lo, hi)

        out_col = f"{col}_norm"
        df[out_col] = norm_vals
        log.info("  %s → %s  finite=%d  range=[%.2f, %.2f]",
                 col, out_col,
                 int(np.isfinite(norm_vals).sum()),
                 float(np.nanmin(norm_vals)), float(np.nanmax(norm_vals)))

    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    df.to_parquet(out_path, engine="pyarrow", index=False)
    log.info("bins_wps_norm 저장: %s", out_path)
    return df