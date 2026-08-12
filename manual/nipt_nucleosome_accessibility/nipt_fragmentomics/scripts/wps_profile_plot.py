"""
wps_profile_plot.py
===================
wps_compute.py 의 npy 캐시 + marker BED 를 입력받아
marker 중심 ±extend bp 범위의 WPS 평균 프로파일을 시각화합니다.

Snyder et al. 2015 Figure 2B 스타일:
  - 상하위 1% 제외 trimmed median 프로파일
  - 그룹 별 비교 오버레이 (y축 공유 고정)
  - WPS / endpoint / coverage 3-track 옵션

사용법
------
python wps_profile_plot.py \\
    --npy SID001.wps_L_long.npy \\
    --marker-bed markers.bed \\
    --extend 2000 \\
    --out SID001_profile.html

# 그룹 비교
python wps_profile_plot.py \\
    --npy low.npy medium.npy high.npy \\
    --labels Low Medium High \\
    --marker-bed markers.bed \\
    --extend 2000 --out compare.html

# 3-track
python wps_profile_plot.py \\
    --npy SID001.wps_L_long.npy \\
    --marker-bed markers.bed \\
    --extend 2000 \\
    --tracks wps_norm endpoint coverage \\
    --out 3track.html
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from typing import Optional

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.signal import savgol_filter

log = logging.getLogger(__name__)

_MARKER_COLS = ["chrom", "start", "end", "marker_id", "marker_type"]

_COLORS = [
    "rgba(220,60,60,0.9)",
    "rgba(50,120,220,0.9)",
    "rgba(50,180,100,0.9)",
    "rgba(180,80,220,0.9)",
    "rgba(220,150,40,0.9)",
    "rgba(40,190,210,0.9)",
]


# ─────────────────────────────────────────────────────────────────────
# 공통 유틸
# ─────────────────────────────────────────────────────────────────────
def _trimmed_median(arr: np.ndarray, q: float = 0.01) -> np.ndarray:
    """
    각 bp 위치(axis=0)에서 상하위 q% 제외 후 median.
    NaN 은 무시합니다.

    Parameters
    ----------
    arr : [n_markers, n_positions]
    q   : 제외 비율 (기본 0.01 = 상하위 1%)
    """
    result = np.full(arr.shape[1], np.nan, dtype=np.float32)
    for i in range(arr.shape[1]):
        col = arr[:, i]
        col = col[np.isfinite(col)]
        if len(col) == 0:
            continue
        lo = np.quantile(col, q)
        hi = np.quantile(col, 1.0 - q)
        trimmed = col[(col >= lo) & (col <= hi)]
        if len(trimmed) > 0:
            result[i] = float(np.median(trimmed))
    return result


def _shared_yrange(
    groups: list[tuple[str, np.ndarray]],
    pad:    float = 0.1,
) -> tuple[float, float]:
    """
    모든 그룹의 프로파일에서 공유 y 범위를 계산합니다.
    상하위 1% 제외 median 값 기준.
    """
    all_vals = []
    for _, profiles in groups:
        if profiles.shape[0] == 0:
            continue
        med = _trimmed_median(profiles)
        finite = med[np.isfinite(med)]
        if len(finite):
            all_vals.extend(finite.tolist())

    if not all_vals:
        return (-1.0, 1.0)

    y_min = float(np.min(all_vals))
    y_max = float(np.max(all_vals))
    margin = (y_max - y_min) * pad
    return (y_min - margin, y_max + margin)


# ─────────────────────────────────────────────────────────────────────
# 데이터 로더
# ─────────────────────────────────────────────────────────────────────
def load_npy(npy_path: str) -> dict:
    return np.load(npy_path, allow_pickle=True).item()


def load_marker_bed(bed_path: str) -> pd.DataFrame:
    df = pd.read_csv(bed_path, sep="\t", header=None,
                     names=_MARKER_COLS, comment="#")
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
# 프로파일 추출
# ─────────────────────────────────────────────────────────────────────
def extract_profiles(
    chrom_arrays: dict[str, np.ndarray],
    marker_df:    pd.DataFrame,
    extend:       int,
) -> np.ndarray:
    """
    marker 중심 ±extend 범위의 배열을 추출합니다.

    Returns
    -------
    profiles : float32 [n_markers, 2*extend+1]
    """
    profile_len = 2 * extend + 1
    profiles    = []

    for row in marker_df.itertuples(index=False):
        chrom  = row.chrom
        center = (int(row.start) + int(row.end)) // 2

        if chrom not in chrom_arrays:
            continue

        arr       = chrom_arrays[chrom]
        chrom_len = len(arr)
        reg_lo    = center - extend
        reg_hi    = center + extend + 1

        src_lo  = max(reg_lo, 0)
        src_hi  = min(reg_hi, chrom_len)
        profile = np.full(profile_len, np.nan, dtype=np.float32)
        dst_lo  = src_lo - reg_lo
        dst_hi  = dst_lo + (src_hi - src_lo)
        profile[dst_lo:dst_hi] = arr[src_lo:src_hi]
        profiles.append(profile)

    if not profiles:
        return np.empty((0, profile_len), dtype=np.float32)
    return np.array(profiles, dtype=np.float32)


# ─────────────────────────────────────────────────────────────────────
# 프로파일 집계: trimmed median + SEM + S-G smoothing
# ─────────────────────────────────────────────────────────────────────
def summarize_profile(
    profiles:  np.ndarray,
    extend:    int,
    sg_window: int = 21,
    sg_poly:   int = 2,
    trim_q:    float = 0.01,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    상하위 1% 제외 trimmed median + SEM + S-G smoothing.

    Returns
    -------
    x    : [-extend, ..., +extend]
    med  : smoothed trimmed median
    sem  : standard error (smoothed)
    """
    n   = profiles.shape[1]
    x   = np.arange(n) - n // 2

    med_raw = _trimmed_median(profiles, q=trim_q)

    # SEM: trimmed 기준 std / sqrt(n_valid)
    sem_raw = np.full(n, np.nan, dtype=np.float32)
    for i in range(n):
        col = profiles[:, i]
        col = col[np.isfinite(col)]
        if len(col) < 2:
            continue
        lo = np.quantile(col, trim_q)
        hi = np.quantile(col, 1.0 - trim_q)
        tr = col[(col >= lo) & (col <= hi)]
        if len(tr) > 1:
            sem_raw[i] = float(tr.std() / np.sqrt(len(tr)))

    # S-G smoothing
    def _sg(arr):
        filled = np.where(np.isfinite(arr), arr, 0.0)
        wl = min(sg_window, len(arr))
        if wl % 2 == 0: wl -= 1
        if wl >= 3 and wl > sg_poly:
            return savgol_filter(filled, window_length=wl,
                                 polyorder=sg_poly).astype(np.float32)
        return arr

    return x, _sg(med_raw), _sg(sem_raw)


# ─────────────────────────────────────────────────────────────────────
# 플롯: 그룹별 오버레이 (y축 공유 고정)
# ─────────────────────────────────────────────────────────────────────
def plot_profiles(
    groups:    list[tuple[str, np.ndarray]],
    extend:    int,
    track:     str  = "wps_norm",
    sg_window: int  = 21,
    title:     str  = "WPS Profile",
    height:    int  = 500,
    show_sem:  bool = True,
    trim_q:    float = 0.01,
) -> go.Figure:
    """
    그룹별 trimmed median 프로파일을 y축 공유 고정으로 오버레이합니다.
    """
    fig = go.Figure()

    # 공유 y 범위 사전 계산
    y_min, y_max = _shared_yrange(groups)

    for i, (label, profiles) in enumerate(groups):
        if profiles.shape[0] == 0:
            log.warning("%s: 프로파일 없음", label)
            continue

        x, med, sem = summarize_profile(
            profiles, extend=extend, sg_window=sg_window, trim_q=trim_q
        )
        color      = _COLORS[i % len(_COLORS)]
        fill_color = color.replace("0.9)", "0.15)")

        if show_sem:
            fig.add_trace(go.Scatter(
                x=np.concatenate([x, x[::-1]]),
                y=np.concatenate([med + sem, (med - sem)[::-1]]),
                fill="toself", fillcolor=fill_color,
                line=dict(color="rgba(0,0,0,0)"),
                showlegend=False, hoverinfo="skip",
            ))

        fig.add_trace(go.Scatter(
            x=x, y=med,
            mode="lines",
            line=dict(color=color, width=2),
            name=f"{label} (n={profiles.shape[0]})",
            hovertemplate=f"{label}<br>pos=%{{x}}<br>{track}=%{{y:.4f}}",
        ))

    fig.add_vline(x=0, line_dash="dash", line_color="black", line_width=1.0,
                  annotation_text="center", annotation_position="top")
    fig.add_hline(y=0, line_dash="dot", line_color="gray", line_width=0.7)

    fig.update_layout(
        title=title, height=height,
        xaxis_title="Distance from marker center (bp)",
        yaxis=dict(
            title=f"Trimmed median WPS ({track})",
            range=[y_min, y_max],   # y축 고정
        ),
        plot_bgcolor="white", paper_bgcolor="white",
        font=dict(size=12),
        legend=dict(x=0.01, y=0.99, bordercolor="lightgray", borderwidth=1),
    )
    fig.update_xaxes(showgrid=True, gridcolor="rgba(200,200,200,0.3)")
    fig.update_yaxes(showgrid=True, gridcolor="rgba(200,200,200,0.3)")

    return fig


# ─────────────────────────────────────────────────────────────────────
# 플롯: 3-track (y축 track별 공유 고정)
# ─────────────────────────────────────────────────────────────────────
def plot_3track(
    npy_data:  dict,
    marker_df: pd.DataFrame,
    extend:    int,
    tracks:    list[str] = None,
    sg_window: int  = 21,
    title:     str  = "WPS 3-Track Profile",
    height:    int  = 800,
    trim_q:    float = 0.01,
) -> go.Figure:
    """
    WPS / endpoint / coverage 3-track. 각 row y축 고정.
    short / long 을 같은 row에 함께 그려서 직접 비교 가능.
    """
    if tracks is None:
        tracks = ["wps_norm", "endpoint", "coverage"]

    track_labels = {
        "wps_raw":  "WPS (raw)",
        "wps_norm": "WPS (normalized)",
        "endpoint": "Endpoint Density",
        "coverage": "Coverage",
    }
    # short/long 컬럼이 npy dict 에 frag 키로 들어있음
    # npy_data = {"wps_raw": {chrom: arr}, "wps_norm": {...}, ...}
    # 각 트랙에서 short/long 분리 처리

    n_rows = len(tracks)
    fig = make_subplots(
        rows=n_rows, cols=1, shared_xaxes=True,
        vertical_spacing=0.06,
        subplot_titles=[track_labels.get(t, t) for t in tracks],
    )

    short_colors = {"wps_norm": "rgba(50,120,220,0.9)",
                    "wps_raw":  "rgba(50,120,220,0.9)",
                    "endpoint": "rgba(220,80,80,0.9)",
                    "coverage": "rgba(50,180,80,0.9)"}
    long_colors  = {"wps_norm": "rgba(180,80,220,0.9)",
                    "wps_raw":  "rgba(180,80,220,0.9)",
                    "endpoint": "rgba(220,140,40,0.9)",
                    "coverage": "rgba(20,130,130,0.9)"}

    for row_idx, track in enumerate(tracks, start=1):
        row_profiles = {}   # "short" / "long" → profiles

        for frag in ("short", "long"):
            # npy 키: "{track}_{frag}" 또는 "{track}" (frag 구분 없는 경우)
            key = f"{track}_{frag}" if f"{track}_{frag}" in npy_data else track
            if key not in npy_data:
                continue
            chrom_arrs = npy_data[key]
            profiles   = extract_profiles(chrom_arrs, marker_df, extend)
            if profiles.shape[0] > 0:
                row_profiles[frag] = profiles

        if not row_profiles:
            # frag 구분 없이 단일 트랙
            if track in npy_data:
                chrom_arrs = npy_data[track]
                profiles   = extract_profiles(chrom_arrs, marker_df, extend)
                row_profiles["all"] = profiles

        if not row_profiles:
            log.warning("'%s' 트랙 데이터 없음", track)
            continue

        # 이 row 의 y 범위: short + long 합쳐서 계산
        all_profiles_list = [(k, v) for k, v in row_profiles.items()]
        y_min_r, y_max_r  = _shared_yrange(all_profiles_list)

        for frag, profiles in row_profiles.items():
            x, med, sem = summarize_profile(
                profiles, extend=extend, sg_window=sg_window, trim_q=trim_q
            )
            if frag == "long":
                color = long_colors.get(track, "rgba(180,80,220,0.9)")
                dash  = "dot"
            else:
                color = short_colors.get(track, "rgba(50,120,220,0.9)")
                dash  = "solid"
            fill_color = color.replace("0.9)", "0.12)")

            fig.add_trace(go.Scatter(
                x=np.concatenate([x, x[::-1]]),
                y=np.concatenate([med + sem, (med - sem)[::-1]]),
                fill="toself", fillcolor=fill_color,
                line=dict(color="rgba(0,0,0,0)"),
                showlegend=False, hoverinfo="skip",
            ), row=row_idx, col=1)

            fig.add_trace(go.Scatter(
                x=x, y=med, mode="lines",
                line=dict(color=color, width=1.5, dash=dash),
                name=f"{track} ({frag})",
                showlegend=(row_idx == 1),
            ), row=row_idx, col=1)

        # y축 고정
        fig.update_yaxes(range=[y_min_r, y_max_r], row=row_idx, col=1)
        fig.add_vline(x=0, line_dash="dash", line_color="black",
                      line_width=0.8, row=row_idx, col=1)
        if track in ("wps_raw", "wps_norm", "endpoint"):
            fig.add_hline(y=0, line_dash="dot", line_color="gray",
                          line_width=0.6, row=row_idx, col=1)

    fig.update_xaxes(title_text="Distance from marker center (bp)",
                     row=n_rows, col=1)
    fig.update_yaxes(showgrid=True, gridcolor="rgba(200,200,200,0.3)")
    fig.update_layout(
        title=title, height=height,
        plot_bgcolor="white", paper_bgcolor="white",
        font=dict(size=11), showlegend=True,
    )
    return fig


# ─────────────────────────────────────────────────────────────────────
# 공개 인터페이스
# ─────────────────────────────────────────────────────────────────────
def run(
    npy_paths:  list[str],
    marker_bed: str,
    out_html:   str,
    labels:     Optional[list[str]] = None,
    extend:     int   = 2000,
    track:      str   = "wps_norm",
    tracks:     Optional[list[str]] = None,
    group_col:  Optional[str] = None,
    sg_window:  int   = 21,
    title:      str   = "WPS Profile",
    height:     int   = 500,
    show_sem:   bool  = True,
    trim_q:     float = 0.01,
) -> go.Figure:
    marker_df = load_marker_bed(marker_bed)

    if labels is None:
        labels = [os.path.basename(p).replace(".npy", "") for p in npy_paths]

    # 3-track 모드
    if tracks is not None:
        assert len(npy_paths) == 1, "3-track 모드는 npy 1개만 지원"
        npy_data = load_npy(npy_paths[0])
        fig = plot_3track(
            npy_data, marker_df, extend=extend,
            tracks=tracks, sg_window=sg_window,
            title=title, height=height, trim_q=trim_q,
        )
        fig.write_html(out_html)
        log.info("저장: %s", out_html)
        return fig

    # 단일/다중 샘플 오버레이 모드
    groups: list[tuple[str, np.ndarray]] = []

    for label, npy_path in zip(labels, npy_paths):
        npy_data   = load_npy(npy_path)
        chrom_arrs = npy_data.get(track, {})

        if group_col and group_col in marker_df.columns:
            for grp_val, grp_df in marker_df.groupby(group_col):
                profiles = extract_profiles(chrom_arrs, grp_df, extend)
                groups.append((f"{label} / {grp_val}", profiles))
        else:
            profiles = extract_profiles(chrom_arrs, marker_df, extend)
            groups.append((label, profiles))

    fig = plot_profiles(
        groups, extend=extend, track=track,
        sg_window=sg_window, title=title,
        height=height, show_sem=show_sem, trim_q=trim_q,
    )
    fig.write_html(out_html)
    log.info("저장: %s", out_html)
    return fig


# ─────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────
def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="wps_profile_plot.py",
        description="WPS npy + marker BED → 평균 프로파일 HTML",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--npy",        nargs="+", required=True)
    p.add_argument("--marker-bed", required=True, dest="marker_bed")
    p.add_argument("--out",        required=True)
    p.add_argument("--labels",     nargs="*", default=None)
    p.add_argument("--extend",     type=int, default=2000)
    p.add_argument("--track",      default="wps_norm",
                   choices=["wps_raw","wps_norm","endpoint","coverage"])
    p.add_argument("--tracks",     nargs="*", default=None,
                   choices=["wps_raw","wps_norm","endpoint","coverage"])
    p.add_argument("--group-col",  default=None, dest="group_col")
    p.add_argument("--sg-window",  type=int, default=21, dest="sg_window")
    p.add_argument("--trim-q",     type=float, default=0.01, dest="trim_q",
                   help="outlier 제외 비율 (기본 0.01 = 상하위 1%%)")
    p.add_argument("--no-sem",     action="store_true", dest="no_sem")
    p.add_argument("--title",      default="WPS Profile")
    p.add_argument("--height",     type=int, default=500)
    p.add_argument("--log-level",  default="INFO", dest="log_level",
                   choices=["DEBUG","INFO","WARNING","ERROR"])
    return p


def main():
    args = _build_parser().parse_args()
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] — %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )
    for p in args.npy:
        if not os.path.exists(p):
            sys.exit(f"npy 없음: {p}")
    if not os.path.exists(args.marker_bed):
        sys.exit(f"marker BED 없음: {args.marker_bed}")
    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)

    run(
        npy_paths  = args.npy,
        marker_bed = args.marker_bed,
        out_html   = args.out,
        labels     = args.labels,
        extend     = args.extend,
        track      = args.track,
        tracks     = args.tracks,
        group_col  = args.group_col,
        sg_window  = args.sg_window,
        title      = args.title,
        height     = args.height,
        show_sem   = not args.no_sem,
        trim_q     = args.trim_q,
    )


if __name__ == "__main__":
    main()