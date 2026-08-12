"""
wps_processor.py
================
wps_compute.py 가 생성한 npy 캐시 + marker BED 를 입력받아
marker 단위 WPS 요약 지표를 생성합니다.

wps_compute.py 역할 (bp 단위, 논문 방식):
  BAM → spanning/endpoints → WPS 배열
  → 1kb 블록 median 차감 (adjusted WPS)
  → Savitzky-Golay smoothing
  → bigWig + npy 저장

이 파일 역할:
  npy 캐시에서 marker 영역 프로파일 추출 → 요약 지표 저장

marker BED 포맷 (header 없음, 탭 구분):
  chrom   start   end   marker_id   marker_type
"""

from __future__ import annotations

import logging
import os
from typing import Optional

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

_MARKER_COLS = ["chrom", "start", "end", "marker_id", "marker_type"]


# ─────────────────────────────────────────────────────────────────────
# marker BED 로더
# ─────────────────────────────────────────────────────────────────────
def load_marker_bed(bed_path: str) -> pd.DataFrame:
    df = pd.read_csv(bed_path, sep="\t", header=None,
                     names=_MARKER_COLS[:5], comment="#")
    df["start"] = df["start"].astype(int)
    df["end"]   = df["end"].astype(int)
    if "marker_id" not in df.columns or df["marker_id"].isna().all():
        df["marker_id"] = (df["chrom"].astype(str) + ":" +
                           df["start"].astype(str) + "-" + df["end"].astype(str))
    if "marker_type" not in df.columns or df["marker_type"].isna().all():
        df["marker_type"] = "marker"
    log.info("marker BED: %d markers", len(df))
    return df


# ─────────────────────────────────────────────────────────────────────
# npy 캐시에서 marker 영역 WPS 추출 및 요약
# ─────────────────────────────────────────────────────────────────────
def _extract_and_summarize(
    chrom_arrays: dict[str, np.ndarray],   # chrom → bp 단위 WPS 배열
    marker_df:    pd.DataFrame,
    extend:       int = 0,                 # marker 중심 ±extend bp 확장 (0이면 marker 영역 그대로)
) -> pd.DataFrame:
    """
    marker 영역 (±extend 포함) 의 WPS 값을 집계합니다.

    집계 지표:
      wps_mean, wps_median, wps_std, wps_max, wps_min,
      wps_oscillation (max-min), n_positions
    """
    records = []
    for row in marker_df.itertuples(index=False):
        chrom     = row.chrom
        m_start   = int(row.start)
        m_end     = int(row.end)
        marker_id = str(row.marker_id)
        mtype     = str(getattr(row, "marker_type", "marker"))

        if chrom not in chrom_arrays:
            continue

        arr       = chrom_arrays[chrom]
        chrom_len = len(arr)

        if extend > 0:
            center  = (m_start + m_end) // 2
            lo      = max(0, center - extend)
            hi      = min(chrom_len, center + extend)
        else:
            lo = max(0, m_start)
            hi = min(chrom_len, m_end)

        if lo >= hi:
            continue

        vals = arr[lo:hi]
        finite = vals[np.isfinite(vals)]
        if len(finite) == 0:
            continue

        records.append({
            "marker_id":       marker_id,
            "marker_type":     mtype,
            "chrom":           chrom,
            "start":           m_start,
            "end":             m_end,
            "n_positions":     len(finite),
            "wps_mean":        float(np.mean(finite)),
            "wps_median":      float(np.median(finite)),
            "wps_std":         float(np.std(finite)),
            "wps_max":         float(np.max(finite)),
            "wps_min":         float(np.min(finite)),
            "wps_oscillation": float(np.max(finite) - np.min(finite)),
        })

    return pd.DataFrame(records)


# ─────────────────────────────────────────────────────────────────────
# 공개 인터페이스
# ─────────────────────────────────────────────────────────────────────
def run(
    wps_npy_path:        Optional[str],    # wps_compute.py 출력 npy
    out_marker_wps_path: str,
    marker_bed_path:     Optional[str] = None,
    wps_track:           str = "wps_norm", # npy dict 내 사용할 트랙
    extend:              int = 0,
    # 하위 호환 파라미터 (무시됨)
    corrected_path:      Optional[str] = None,
    out_peaks_path:      Optional[str] = None,
    sg_window:           int = 5,
    sg_polyorder:        int = 2,
) -> dict:
    """
    wps_compute.py npy + marker BED → marker_wps_summary.parquet

    wps_npy_path 없으면 marker WPS 생략.
    """
    if not marker_bed_path or not os.path.exists(marker_bed_path):
        log.info("marker BED 없음 — marker WPS 생략")
        return {"marker_wps": pd.DataFrame()}

    marker_df = load_marker_bed(marker_bed_path)

    if not wps_npy_path or not os.path.exists(wps_npy_path):
        log.warning("WPS npy 없음: %s — marker WPS 생략", wps_npy_path)
        log.warning("wps_compute.py 를 먼저 실행하세요.")
        return {"marker_wps": pd.DataFrame()}

    # npy 로드
    data = np.load(wps_npy_path, allow_pickle=True).item()
    chrom_arrays = data.get(wps_track, {})
    if not chrom_arrays:
        # frag 별로 분리된 경우 시도
        for key in data:
            if isinstance(data[key], dict):
                chrom_arrays = data[key]
                log.info("npy track '%s' 사용", key)
                break

    log.info("marker WPS 추출: %d markers, track=%s", len(marker_df), wps_track)
    summary_df = _extract_and_summarize(chrom_arrays, marker_df, extend=extend)

    os.makedirs(os.path.dirname(os.path.abspath(out_marker_wps_path)), exist_ok=True)
    summary_df.to_parquet(out_marker_wps_path, engine="pyarrow", index=False)
    log.info("marker WPS 저장: %s (%d markers)", out_marker_wps_path, len(summary_df))

    return {"marker_wps": summary_df}