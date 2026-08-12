"""
cnv_caller.py
=============
bins_corrected.parquet 에서 short / long 각각의 CNV 를 호출합니다.

처리 순서
---------
1. LOO (Leave-One-Out) 정규화
     - 각 상염색체를 정규화할 때 자기 자신을 제외한 나머지 상염색체 median 사용
     - 거대 이수성 염색체(chr21 trisomy 등)가 전체 baseline 을 왜곡하는 것을 방지
     - 결과: copy_number_signal (정상=2.0), log2_chrom_norm (정상=0.0)
2. 성염색체 판별 (XX/XY)
     - copy_number_signal 기반 X/Y median 으로 판정
3. MAD 기반 robust z-score (LOO 정규화된 log2_chrom_norm 기준)
4. 성염색체 z-score 보정 (FF 반영)
5. CBS 근사 세그멘테이션 (short log2_chrom_norm 기준)
6. gain / loss / normal CNV call

출력 추가 컬럼
--------------
  short_copy_number / long_copy_number   : LOO copy number signal (정상=2.0)
  short_log2_norm   / long_log2_norm     : log2(copy_number/2), 정상=0.0
  short_zscore      / long_zscore
  short_cnv_call    / long_cnv_call
  segment_id
  sex_call                               : "XX" | "XY"
"""

from __future__ import annotations

import json
import logging
import os
from typing import Optional

import numpy as np
import pandas as pd

from nipt_fragmentomics.core.constants import (
    ZSCORE_GAIN, ZSCORE_LOSS, SEG_MIN_BINS, SEG_ALPHA,
)

log = logging.getLogger(__name__)

SEX_CHROMS = {"chrX", "chrY"}


# ─────────────────────────────────────────────────────────────────────
# LOO 정규화
# ─────────────────────────────────────────────────────────────────────
def loo_normalize(
    df:         pd.DataFrame,
    value_col:  str,
) -> tuple[pd.DataFrame, str, float, float]:
    """
    Leave-One-Out 정규화.

    mappability_pass=True 인 bin 만 baseline 계산에 사용합니다.
    mappability_pass=False 인 bin 은 copy_number / log2_norm 을 NaN 으로 설정합니다.
    """
    df_norm = df.copy()

    # log2 → linear 변환
    if "log2" in value_col:
        linear_vals = np.exp2(df_norm[value_col].values.astype(np.float64))
    else:
        linear_vals = df_norm[value_col].values.astype(np.float64)

    # mappability_pass + 유한값 마스크
    if "mappability_pass" in df_norm.columns:
        quality_mask = df_norm["mappability_pass"].values.astype(bool)
    else:
        quality_mask = np.ones(len(df_norm), dtype=bool)

    # gc=0 또는 count=0 bin 추가 제외
    if "gc" in df_norm.columns:
        gc_vals = df_norm["gc"].values.astype(np.float64)
        quality_mask &= np.isfinite(gc_vals) & (gc_vals > 0)

    quality_mask &= np.isfinite(linear_vals) & (linear_vals > 0)

    auto_mask    = ~df_norm["chrom"].isin(SEX_CHROMS)
    auto_quality = auto_mask & quality_mask   # 계산에 사용할 bin

    # 전체 autosome quality-pass 의 global median
    global_auto_med = float(np.nanmedian(linear_vals[auto_quality]))
    if not np.isfinite(global_auto_med) or global_auto_med == 0:
        global_auto_med = 1e-9

    loo_baseline = np.full(len(df_norm), global_auto_med, dtype=np.float64)

    # 상염색체별 LOO baseline (quality-pass bin 만)
    auto_chroms = df_norm.loc[auto_quality, "chrom"].unique()
    for chrom in auto_chroms:
        chrom_mask = df_norm["chrom"] == chrom
        loo_mask   = auto_quality & (~chrom_mask)
        loo_med    = float(np.nanmedian(linear_vals[loo_mask]))
        if not np.isfinite(loo_med) or loo_med == 0:
            loo_med = global_auto_med
        loo_baseline[chrom_mask.values] = loo_med

    # copy number signal
    cn_signal = 2.0 * (linear_vals / loo_baseline)

    # quality_pass=False bin → NaN
    cn_signal[~quality_mask] = np.nan

    df_norm[f"{value_col}_copy_number"] = cn_signal.astype(np.float32)

    with np.errstate(divide="ignore", invalid="ignore"):
        log2_norm = np.log2(cn_signal / 2.0)
    log2_norm = np.where(np.isfinite(log2_norm), log2_norm, np.nan)
    # mappability_pass=False 는 NaN 유지
    df_norm[f"{value_col}_log2_norm"] = log2_norm.astype(np.float32)

    # 성염색체 copy number median (quality_mask 무관하게 성염색체 전체 사용)
    x_mask = df_norm["chrom"] == "chrX"
    y_mask = df_norm["chrom"] == "chrY"
    x_cn   = float(np.nanmedian(cn_signal[x_mask.values])) if x_mask.any() else 0.0
    y_cn   = float(np.nanmedian(cn_signal[y_mask.values])) if y_mask.any() else 0.0
    sex_call = "XX" if x_cn >= 1.5 else "XY"

    n_total   = len(df_norm)
    n_quality = int(quality_mask.sum())
    log.info(
        "LOO 정규화 완료 [%s]  quality=%d/%d  sex=%s  X_CN=%.3f  Y_CN=%.3f",
        value_col, n_quality, n_total, sex_call, x_cn, y_cn,
    )
    return df_norm, sex_call, x_cn, y_cn


# ─────────────────────────────────────────────────────────────────────
# MAD 기반 robust z-score
# ─────────────────────────────────────────────────────────────────────
def _robust_zscore(
    values:   np.ndarray,
    ref_mask: np.ndarray,
) -> np.ndarray:
    """ref_mask 기준 median/MAD 로 z-score 계산."""
    ref = values[ref_mask & np.isfinite(values)]
    if len(ref) == 0:
        return np.full_like(values, np.nan, dtype=np.float32)
    med   = float(np.median(ref))
    mad   = float(np.median(np.abs(ref - med)))
    scale = mad * 1.4826 if mad > 1e-6 else (float(ref.std()) or 1.0)
    return ((values - med) / scale).astype(np.float32)


# ─────────────────────────────────────────────────────────────────────
# CBS 근사 세그멘테이션
# ─────────────────────────────────────────────────────────────────────
def _segment_chrom(
    zscore:   np.ndarray,
    min_bins: int   = SEG_MIN_BINS,
    alpha:    float = SEG_ALPHA,
) -> np.ndarray:
    """재귀 Welch t-test 기반 CBS 근사 세그멘테이션."""
    from scipy.stats import ttest_ind

    n       = len(zscore)
    seg_id  = np.zeros(n, dtype=np.int32)
    stack   = [(0, n, 0)]
    next_id = 1

    while stack:
        lo, hi, sid = stack.pop()
        if hi - lo < 2 * min_bins:
            seg_id[lo:hi] = sid
            continue

        chunk   = zscore[lo:hi]
        best_p  = 1.0
        best_cp = lo + min_bins

        for cp in range(lo + min_bins, hi - min_bins + 1):
            lf = chunk[:cp - lo];    lf = lf[np.isfinite(lf)]
            rf = chunk[cp - lo:];    rf = rf[np.isfinite(rf)]
            if len(lf) < 2 or len(rf) < 2:
                continue
            _, p = ttest_ind(lf, rf, equal_var=False)
            if p < best_p:
                best_p, best_cp = p, cp

        if best_p < alpha:
            seg_id[lo:best_cp] = sid
            stack.append((best_cp, hi, next_id))
            next_id += 1
        else:
            seg_id[lo:hi] = sid

    return seg_id


# ─────────────────────────────────────────────────────────────────────
# CNV call
# ─────────────────────────────────────────────────────────────────────
def _cnv_call(z: float, gain: float, loss: float) -> str:
    if not np.isfinite(z): return "unknown"
    if z >= gain:          return "gain"
    if z <= loss:          return "loss"
    return "normal"


# ─────────────────────────────────────────────────────────────────────
# 공개 인터페이스
# ─────────────────────────────────────────────────────────────────────
def run(
    corrected_path: str,
    out_path:       str,
    ff_json_path:   Optional[str] = None,
    zscore_gain:    float = ZSCORE_GAIN,
    zscore_loss:    float = ZSCORE_LOSS,
    seg_min_bins:   int   = SEG_MIN_BINS,
    seg_alpha:      float = SEG_ALPHA,
) -> pd.DataFrame:
    """
    bins_corrected.parquet → cnv_calls.parquet

    LOO 정규화 → z-score → 세그멘테이션 → CNV call
    """
    df = pd.read_parquet(corrected_path)
    log.info("CNV caller: %d bins", len(df))

    # FF 로드
    ff: Optional[float] = None
    if ff_json_path and os.path.exists(ff_json_path):
        with open(ff_json_path) as f:
            ff_data = json.load(f)
        ff = ff_data.get("consensus_ff")
        if ff:
            log.info("FF=%.4f 성염색체 보정 적용", ff)

    is_sex  = df["chrom"].isin(SEX_CHROMS).values
    is_auto = ~is_sex
    mp_pass = df["mappability_pass"].values
    ref_mask = is_auto & mp_pass   # z-score 기준 mask

    sex_calls = []

    for frag in ("short", "long"):
        log2_col = f"log2_corrected_{frag}"
        if log2_col not in df.columns:
            log.warning("%s 컬럼 없음 — 건너뜀", log2_col)
            df[f"{frag}_copy_number"] = np.nan
            df[f"{frag}_log2_norm"]   = np.nan
            df[f"{frag}_zscore"]      = np.nan
            df[f"{frag}_cnv_call"]    = "unknown"
            continue

        # ── Step 1: LOO 정규화 ────────────────────────────────────
        df, sex_call, x_cn, y_cn = loo_normalize(df, log2_col)

        # 컬럼명 정리: loo_normalize 가 {log2_col}_copy_number 등으로 생성
        df = df.rename(columns={
            f"{log2_col}_copy_number": f"{frag}_copy_number",
            f"{log2_col}_log2_norm":   f"{frag}_log2_norm",
        })
        sex_calls.append(sex_call)

        # ── Step 2: z-score (LOO log2_norm 기준) ─────────────────
        log2_norm = df[f"{frag}_log2_norm"].values.astype(np.float32)
        zscores   = _robust_zscore(log2_norm, ref_mask)

        # ── Step 3: 성염색체 z-score 보정 (FF 반영) ──────────────
        if ff is not None:
            ref_vals = log2_norm[ref_mask & np.isfinite(log2_norm)]
            med   = float(np.median(ref_vals)) if len(ref_vals) else 0.0
            mad   = float(np.median(np.abs(ref_vals - med))) if len(ref_vals) else 1.0
            scale = mad * 1.4826 if mad > 1e-6 else 1.0

            # chrX: FF 있으면 기대 copy number 가 2-FF (모체) + FF/2 (태아)
            # log2 기준으로 expected shift 계산
            for sex_chr, expected_cn in [
                ("chrX", 2.0 - ff / 2.0),   # 여아 태아: ~2.0, 남아 태아: ~1.5
                ("chrY", ff / 2.0),          # 남아 태아: ~FF/2
            ]:
                mask = df["chrom"] == sex_chr
                if not mask.any():
                    continue
                expected_log2_norm = np.log2(max(expected_cn, 1e-6) / 2.0)
                zscores[mask.values] = (
                    (log2_norm[mask.values] - expected_log2_norm - med) / scale
                ).astype(np.float32)

        df[f"{frag}_zscore"] = zscores

        # ── Step 4: CNV call ──────────────────────────────────────
        df[f"{frag}_cnv_call"] = [
            _cnv_call(float(z), zscore_gain, zscore_loss)
            for z in zscores
        ]

    # sex_call 통합 (short/long 일치하면 그대로, 불일치 시 short 우선)
    final_sex = sex_calls[0] if sex_calls else "Unknown"
    df["sex_call"] = final_sex
    log.info("최종 성별 판정: %s", final_sex)

    # ── Step 5: 세그멘테이션 (short log2_norm 기준) ───────────────
    seg_ids    = np.full(len(df), -1, dtype=np.int32)
    seg_offset = 0
    for chrom, cdf in df.groupby("chrom"):
        col = "short_log2_norm"
        if col not in cdf.columns:
            continue
        vals = np.where(
            np.isfinite(cdf[col].values),
            cdf[col].values, 0.0,
        )
        segs = _segment_chrom(vals, min_bins=seg_min_bins, alpha=seg_alpha)
        seg_ids[cdf.index] = segs + seg_offset
        seg_offset += int(segs.max()) + 1 if len(segs) else 0

    df["segment_id"] = seg_ids

    # ── 출력 컬럼 정리 ────────────────────────────────────────────
    keep = [
        "chrom", "start", "end",
        "short_count",          "long_count",
        "short_breadth",        "long_breadth",
        "short_wps_L_corrected","short_wps_S_corrected",
        "long_wps_L_corrected", "long_wps_S_corrected",
        # LOO 정규화 결과
        "short_copy_number",    "long_copy_number",
        "short_log2_norm",      "long_log2_norm",
        # 원본 GC 보정 log2 (참조용)
        "log2_corrected_short", "log2_corrected_long", "log2_ratio",
        # CNV
        "short_zscore",    "long_zscore",
        "short_cnv_call",  "long_cnv_call",
        "segment_id",      "sex_call",
        "gc", "mappability", "mappability_pass",
    ]
    out_df = df[[c for c in keep if c in df.columns]].copy()

    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    out_df.to_parquet(out_path, engine="pyarrow", index=False)
    log.info("CNV calls 저장: %s (%d bins)", out_path, len(out_df))
    return out_df