"""
wps_processor.py
================
wps_compute.py 가 생성한 npy 에서 marker 단위 WPS 요약 지표를 생성합니다.

npy 구조 (원저자 방식)
----------------------
{
  "wps_norm": {marker_id: adjusted_WPS_array [2*extend+1]},
  "mode":     "L",
  "frag":     "long",
  "extend":   1000,
}

adjusted_WPS = (WPS - windowMedian) / Coverage * 100

요약 지표
---------
  wps_mean, wps_median, wps_center_mean (중심 ±200bp),
  wps_peak, wps_trough, wps_oscillation, n_positions
"""

from __future__ import annotations

import logging
import os
from typing import Optional

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

_MARKER_COLS = ["chrom", "start", "end", "marker_id", "marker_type"]


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
    return df


def run(
    wps_npy_path:        Optional[str],
    out_marker_wps_path: str,
    marker_bed_path:     Optional[str] = None,
    extend:              int = 1000,
    center_win:          int = 200,    # wps_center_mean 계산 범위 ±center_win bp
    # 하위 호환 파라미터
    corrected_path:      Optional[str] = None,
    out_peaks_path:      Optional[str] = None,
    sg_window:           int = 5,
    sg_polyorder:        int = 2,
) -> dict:
    """npy + marker BED → marker_wps_summary.parquet"""

    if not wps_npy_path or not os.path.exists(wps_npy_path):
        log.warning("WPS npy 없음 — marker WPS 생략: %s", wps_npy_path)
        return {"marker_wps": pd.DataFrame()}

    # npy 로드
    data = np.load(wps_npy_path, allow_pickle=True).item()
    wps_dict: dict[str, np.ndarray] = data.get("wps_norm", {})
    npy_extend = data.get("extend", extend)

    if not wps_dict:
        log.warning("npy 에 wps_norm 데이터 없음")
        return {"marker_wps": pd.DataFrame()}

    # marker BED (없으면 npy 키만으로 요약)
    if marker_bed_path and os.path.exists(marker_bed_path):
        marker_df = load_marker_bed(marker_bed_path)
        meta = {row.marker_id: row
                for row in marker_df.itertuples(index=False)}
    else:
        meta = {}

    # 요약 계산
    n_pos  = 2 * npy_extend + 1
    center = n_pos // 2
    c_lo   = max(0, center - center_win)
    c_hi   = min(n_pos, center + center_win + 1)

    records = []
    for marker_id, arr in wps_dict.items():
        if arr is None or len(arr) == 0:
            continue
        finite = arr[np.isfinite(arr)]
        if len(finite) == 0:
            continue

        center_vals = arr[c_lo:c_hi]
        center_finite = center_vals[np.isfinite(center_vals)]

        row = meta.get(marker_id)
        rec = {
            "marker_id":        marker_id,
            "marker_type":      getattr(row, "marker_type", "marker") if row else "marker",
            "chrom":            getattr(row, "chrom",  "") if row else "",
            "start":            getattr(row, "start",  -1) if row else -1,
            "end":              getattr(row, "end",    -1) if row else -1,
            "n_positions":      len(finite),
            "wps_mean":         float(np.mean(finite)),
            "wps_median":       float(np.median(finite)),
            "wps_center_mean":  float(np.mean(center_finite)) if len(center_finite) else np.nan,
            "wps_peak":         float(np.max(finite)),
            "wps_trough":       float(np.min(finite)),
            "wps_oscillation":  float(np.max(finite) - np.min(finite)),
        }
        records.append(rec)

    summary_df = pd.DataFrame(records)
    log.info("marker WPS 요약: %d markers", len(summary_df))

    os.makedirs(os.path.dirname(os.path.abspath(out_marker_wps_path)), exist_ok=True)
    summary_df.to_parquet(out_marker_wps_path, engine="pyarrow", index=False)
    log.info("marker WPS 저장: %s", out_marker_wps_path)

    return {"marker_wps": summary_df}
