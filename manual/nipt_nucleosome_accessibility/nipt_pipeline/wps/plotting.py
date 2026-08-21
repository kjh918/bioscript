"""
wps/plotting.py
===============
WPS 시각화 2종

1. plot_genome_wide()
   - 전체 염색체 연속 X축
   - raw / adjusted(norm) / coverage 3-track subplot (1 PNG per mode)

2. plot_marker_aggregate()
   - marker_stats.parquet / _profiles.npy 기반
   - L × (raw / adjusted / cov) + S × (raw / adjusted / cov) = 6-panel
"""
from __future__ import annotations

import logging
import os
from typing import Optional

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

from .utils import (
    STANDARD_CHROMS, HG38_CHROM_SIZES, WPS_PARAMS,
    load_track, read_manifest,
)

log = logging.getLogger(__name__)

# 염색체 배경 교대색
_BG = ["#F7F7F7", "#EFEFEF"]

# 트랙별 색상
_TRACK_COLOR = {
    "raw":      "#455A64",  # blue-grey
    "norm":     "#1E88E5",  # blue
    "cov":      "#43A047",  # green
    "frag_cov": "#8BC34A",  # light green
}


# ─────────────────────────────────────────────────────────────────────
# 1. Genome-wide line plot
# ─────────────────────────────────────────────────────────────────────
def _subsample(arr: np.ndarray, max_pts: int = 500_000) -> np.ndarray:
    """너무 많은 포인트는 균등 subsample하여 렌더링 속도 확보"""
    if arr.size <= max_pts:
        return arr.astype(np.float32)
    step = arr.size // max_pts
    return arr[::step].astype(np.float32)


def _rolling_mean(arr: np.ndarray, window: int) -> np.ndarray:
    """
    시각화용 rolling mean — cumsum 기반 O(n), 빠름.
    window는 원본 해상도(bp) 기준.
    경계는 available 포인트 수로 나눠 edge bias 방지.
    """
    if window <= 1 or arr.size < window:
        return arr.astype(np.float32)

    a      = arr.astype(np.float64)
    # NaN → 0, count 배열로 나누기
    valid  = np.isfinite(a).astype(np.float64)
    filled = np.where(np.isfinite(a), a, 0.0)

    half = window // 2

    # cumsum으로 sliding sum
    cs      = np.concatenate([[0.0], np.cumsum(filled)])
    cv      = np.concatenate([[0.0], np.cumsum(valid)])

    # 각 포인트에 대해 [i-half, i+half] 범위 합산
    lo = np.maximum(np.arange(len(a)) - half, 0)
    hi = np.minimum(np.arange(len(a)) + half + 1, len(a))

    sums   = cs[hi] - cs[lo]
    counts = cv[hi] - cv[lo]

    result = np.where(counts > 0, sums / counts, 0.0)
    return result.astype(np.float32)


def _robust_ylim(
    chunks: list[np.ndarray],
    positive_only: bool,
    pad_ratio: float = 0.15,
) -> tuple[float, float]:
    """
    IQR 기반 robust ylim 계산.
    이상치 1~2개로 Y축 스케일이 압축되는 문제 방지.
    """
    if not chunks:
        return (0.0, 1.0)
    all_y = np.concatenate(chunks)
    all_y = all_y[np.isfinite(all_y)]
    if all_y.size == 0:
        return (0.0, 1.0)

    if positive_only:
        p99 = float(np.percentile(all_y, 99))
        return (0.0, max(p99 * 1.2, 1.0))

    q25, med, q75 = np.percentile(all_y, [25, 50, 75])
    iqr  = q75 - q25 or 1.0
    ymin = med - 3.5 * iqr
    ymax = med + 3.5 * iqr
    pad  = (ymax - ymin) * pad_ratio
    return (ymin - pad, ymax + pad)


def _collect_mode_data(
    wps_dir:    str,
    mode:       str,
    prefix:     str,
    chroms:     list[str],
    downsample: int,
    win_smooth: int,
) -> tuple[dict[str, list[np.ndarray]], list[np.ndarray], list[float], dict[str, float]]:
    """
    단일 mode(L or S) 의 track 데이터를 수집.
    Returns (track_data, x_data, chrom_boundaries, chrom_centers)
    """
    manifest = read_manifest(wps_dir, mode)

    track_data: dict[str, list[np.ndarray]] = {
        "raw": [], "norm": [], "frag_cov": []
    }
    x_data:           list[np.ndarray] = []
    chrom_boundaries: list[float]      = [0.0]
    chrom_centers:    dict[str, float] = {}
    x_cursor = 0.0

    for chrom in chroms:
        chrom_len = HG38_CHROM_SIZES.get(chrom, 0)
        if chrom_len == 0:
            continue

        arrays: dict[str, Optional[np.ndarray]] = {}
        for metric in ("raw", "norm", "frag_cov"):
            arrays[metric] = load_track(
                wps_dir  = wps_dir,
                mode     = mode,
                chrom    = chrom,
                metric   = metric,
                prefix   = prefix,
                manifest = manifest,
            )

        ref_arr = arrays.get("raw")
        n = int(ref_arr.size) if ref_arr is not None else chrom_len
        chrom_centers[chrom] = x_cursor + n / 2

        for metric in ("raw", "norm", "frag_cov"):
            arr = arrays.get(metric)
            if arr is not None:
                arr_f = arr.astype(np.float32)
                if win_smooth > 1:
                    arr_f = _rolling_mean(arr_f, win_smooth)
                y = _subsample(arr_f, downsample)
            else:
                y = np.full(min(n, downsample), np.nan, dtype=np.float32)
            track_data[metric].append(y)

        step = n / max(len(track_data["raw"][-1]), 1)
        xs   = np.arange(len(track_data["raw"][-1])) * step + x_cursor
        x_data.append(xs)

        x_cursor += n
        chrom_boundaries.append(x_cursor)

    return track_data, x_data, chrom_boundaries, chrom_centers


def _draw_genome_panels(
    axes:             list,
    panel_specs:      list[tuple],
    track_data:       dict[str, list[np.ndarray]],
    x_data:           list[np.ndarray],
    chrom_boundaries: list[float],
    chrom_centers:    dict[str, float],
    target:           list[str],
    x_cursor:         float,
    last_ax:          bool = False,
) -> None:
    """공통 패널 드로잉 헬퍼"""
    for ax, (metric, ylabel, positive_only, color) in zip(axes, panel_specs):
        ax.set_facecolor("white")

        for idx, chrom in enumerate(target):
            if idx >= len(chrom_boundaries) - 1:
                break
            ax.axvspan(chrom_boundaries[idx], chrom_boundaries[idx + 1],
                       color=_BG[idx % 2], alpha=0.6, zorder=0)

        for xs, ys in zip(x_data, track_data[metric]):
            finite = np.isfinite(ys)
            if not finite.any():
                continue
            ax.plot(xs, ys, linewidth=0.4, color=color, alpha=0.85, rasterized=True)

        ax.axhline(0, linewidth=0.5, color="#9E9E9E", linestyle="--", alpha=0.5)
        for xb in chrom_boundaries[1:-1]:
            ax.axvline(xb, color="#BDBDBD", linewidth=0.4, zorder=1)

        ax.set_ylabel(ylabel, fontsize=8)
        ax.tick_params(axis="y", labelsize=7)
        ax.spines[["top", "right"]].set_visible(False)

        _chunks = [y[np.isfinite(y)] for y in track_data[metric] if np.any(np.isfinite(y))]
        ymin, ymax = _robust_ylim(_chunks, positive_only)
        ax.set_ylim(ymin, ymax)

    if last_ax:
        axes[-1].set_xticks([chrom_centers[c] for c in target if c in chrom_centers])
        axes[-1].set_xticklabels(
            [c.replace("chr", "") for c in target if c in chrom_centers],
            fontsize=7, rotation=0,
        )
        axes[-1].set_xlabel("Chromosome", fontsize=9)
    axes[-1].set_xlim(0, x_cursor)


def plot_genome_wide(
    wps_dir:     str,
    output_path: str,
    sample_id:   str   = "",
    prefix:      str   = "",
    chroms:      Optional[list[str]] = None,
    downsample:  int   = 500_000,
    win_smooth:  int   = 10000,
) -> None:
    """
    L + S WPS genome-wide 통합 plot (6패널, PNG 1개)

    패널 (위→아래):
      L Raw WPS / L Adjusted WPS / L Fragment coverage
      S Raw WPS / S Adjusted WPS / S Fragment coverage

    Parameters
    ----------
    wps_dir     : compute.run() out_dir (L/ S/ 상위)
    output_path : PNG 저장 경로
    sample_id   : 타이틀 표시용
    prefix      : compute prefix (sample_id)
    chroms      : 표시할 염색체 (None이면 standard 전체)
    downsample  : 염색체당 최대 표시 포인트
    win_smooth  : 시각화용 rolling mean window (bp)
    """
    target = chroms or STANDARD_CHROMS

    log.info("genome-wide plot 데이터 수집 중 (L)...")
    l_data, l_xs, l_bounds, l_centers = _collect_mode_data(
        wps_dir, "L", prefix, target, downsample, win_smooth
    )
    log.info("genome-wide plot 데이터 수집 중 (S)...")
    s_data, s_xs, s_bounds, s_centers = _collect_mode_data(
        wps_dir, "S", prefix, target, downsample, win_smooth
    )

    if not l_xs and not s_xs:
        log.warning("genome-wide plot: 로드된 데이터 없음")
        return

    # x_cursor: L/S는 같은 chrom 크기 기준이므로 동일
    x_cursor = l_bounds[-1] if l_bounds else s_bounds[-1]

    l_prm = WPS_PARAMS["L"]
    s_prm = WPS_PARAMS["S"]

    panel_specs_L = [
        ("raw",      f"L Raw WPS  ({l_prm['frag_min']}-{l_prm['frag_max']}bp)",      False, _TRACK_COLOR["raw"]),
        ("norm",     f"L Adjusted WPS",                                               False, _TRACK_COLOR["norm"]),
        ("frag_cov", f"L Frag coverage",                                              True,  _TRACK_COLOR["frag_cov"]),
    ]
    panel_specs_S = [
        ("raw",      f"S Raw WPS  ({s_prm['frag_min']}-{s_prm['frag_max']}bp)",      False, "#E53935"),
        ("norm",     f"S Adjusted WPS",                                               False, "#FB8C00"),
        ("frag_cov", f"S Frag coverage",                                              True,  "#A5D6A7"),
    ]

    fig, axes = plt.subplots(6, 1, figsize=(22, 14), sharex=True, facecolor="white")
    fig.subplots_adjust(hspace=0.06)
    fig.suptitle(
        f"{sample_id} — WPS genome-wide  L({l_prm['frag_min']}-{l_prm['frag_max']}bp)"
        f"  S({s_prm['frag_min']}-{s_prm['frag_max']}bp)",
        fontsize=11, fontweight="bold", y=0.995,
    )

    _draw_genome_panels(
        axes[:3], panel_specs_L, l_data, l_xs, l_bounds, l_centers,
        target, x_cursor, last_ax=False,
    )
    _draw_genome_panels(
        axes[3:], panel_specs_S, s_data, s_xs, s_bounds, s_centers,
        target, x_cursor, last_ax=True,
    )

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    fig.savefig(output_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    log.info("genome-wide plot 저장: %s", output_path)

def plot_genome_wide_by_chrom(
    wps_dir:    str,
    output_dir: str,
    sample_id:  str  = "",
    prefix:     str  = "",
    chroms:     Optional[list[str]] = None,
    downsample: int  = 500_000,
    win_smooth: int  = 10000,
) -> list[str]:
    """
    chromosome별 L+S WPS 6패널 PNG 저장.

    Parameters
    ----------
    wps_dir    : compute.run() out_dir
    output_dir : PNG 저장 디렉토리
    sample_id  : 타이틀 표기용
    prefix     : compute prefix
    chroms     : 대상 염색체 (None이면 standard 전체)
    downsample : 염색체당 최대 포인트
    win_smooth : rolling mean window (bp)

    Returns
    -------
    저장된 PNG 경로 리스트
    """
    target = chroms or STANDARD_CHROMS
    os.makedirs(output_dir, exist_ok=True)

    manifests = {mode: read_manifest(wps_dir, mode) for mode in ("L", "S")}
    l_prm = WPS_PARAMS["L"]
    s_prm = WPS_PARAMS["S"]

    panel_specs_L = [
        ("raw",      f"L Raw WPS  ({l_prm['frag_min']}-{l_prm['frag_max']}bp)", False, _TRACK_COLOR["raw"]),
        ("norm",     "L Adjusted WPS",                                           False, _TRACK_COLOR["norm"]),
        ("frag_cov", "L Frag coverage",                                          True,  _TRACK_COLOR["frag_cov"]),
    ]
    panel_specs_S = [
        ("raw",      f"S Raw WPS  ({s_prm['frag_min']}-{s_prm['frag_max']}bp)", False, "#E53935"),
        ("norm",     "S Adjusted WPS",                                           False, "#FB8C00"),
        ("frag_cov", "S Frag coverage",                                          True,  "#A5D6A7"),
    ]

    saved: list[str] = []

    for chrom in target:
        chrom_len = HG38_CHROM_SIZES.get(chrom, 0)
        if chrom_len == 0:
            continue

        # L / S track 로드
        mode_arrays: dict[str, dict[str, Optional[np.ndarray]]] = {}
        any_data = False
        for mode in ("L", "S"):
            arrays = {}
            for metric in ("raw", "norm", "frag_cov"):
                arr = load_track(
                    wps_dir  = wps_dir,
                    mode     = mode,
                    chrom    = chrom,
                    metric   = metric,
                    prefix   = prefix,
                    manifest = manifests[mode],
                )
                arrays[metric] = arr
                if arr is not None:
                    any_data = True
            mode_arrays[mode] = arrays

        if not any_data:
            log.warning("[%s] 모든 track 없음 — skip", chrom)
            continue

        # 데이터 준비 (rolling mean → subsample)
        def _prep(arr: Optional[np.ndarray], n: int) -> np.ndarray:
            if arr is None:
                return np.full(min(n, downsample), np.nan, dtype=np.float32)
            arr_f = arr.astype(np.float32)
            if win_smooth > 1:
                arr_f = _rolling_mean(arr_f, win_smooth)
            return _subsample(arr_f, downsample)

        n = chrom_len
        l_arrays = mode_arrays["L"]
        s_arrays = mode_arrays["S"]

        # x축
        ref = next((a for a in l_arrays.values() if a is not None), None)
        n   = int(ref.size) if ref is not None else chrom_len
        xs  = np.linspace(0, chrom_len, min(n, downsample))

        fig, axes = plt.subplots(6, 1, figsize=(18, 10),
                                  sharex=True, facecolor="white")
        fig.subplots_adjust(hspace=0.06)
        fig.suptitle(
            f"{sample_id} — {chrom}  WPS  "
            f"L({l_prm['frag_min']}-{l_prm['frag_max']}bp) "
            f"S({s_prm['frag_min']}-{s_prm['frag_max']}bp)",
            fontsize=10, fontweight="bold", y=0.995,
        )

        for ax, (metric, ylabel, positive_only, color) in zip(axes[:3], panel_specs_L):
            y = _prep(l_arrays.get(metric), n)
            _x = xs[:len(y)]
            ax.set_facecolor("white")
            finite = np.isfinite(y)
            if finite.any():
                ax.plot(_x, y, linewidth=0.5, color=color, alpha=0.85, rasterized=True)
            ax.axhline(0, linewidth=0.5, color="#9E9E9E", linestyle="--", alpha=0.5)
            ax.set_ylabel(ylabel, fontsize=8)
            ax.tick_params(axis="y", labelsize=7)
            ax.spines[["top", "right"]].set_visible(False)
            chunks = [y[finite]] if finite.any() else []
            ymin, ymax = _robust_ylim(chunks, positive_only)
            ax.set_ylim(ymin, ymax)

        for ax, (metric, ylabel, positive_only, color) in zip(axes[3:], panel_specs_S):
            y = _prep(s_arrays.get(metric), n)
            _x = xs[:len(y)]
            ax.set_facecolor("white")
            finite = np.isfinite(y)
            if finite.any():
                ax.plot(_x, y, linewidth=0.5, color=color, alpha=0.85, rasterized=True)
            ax.axhline(0, linewidth=0.5, color="#9E9E9E", linestyle="--", alpha=0.5)
            ax.set_ylabel(ylabel, fontsize=8)
            ax.tick_params(axis="y", labelsize=7)
            ax.spines[["top", "right"]].set_visible(False)
            chunks = [y[finite]] if finite.any() else []
            ymin, ymax = _robust_ylim(chunks, positive_only)
            ax.set_ylim(ymin, ymax)

        # X축: Mbp 단위
        axes[-1].set_xlabel(f"{chrom} position (Mbp)", fontsize=9)
        axes[-1].xaxis.set_major_formatter(
            mticker.FuncFormatter(lambda x, _: f"{x/1e6:.0f}")
        )
        axes[-1].set_xlim(0, chrom_len)
        axes[-1].tick_params(axis="x", labelsize=7)

        out_png = os.path.join(output_dir, f"{sample_id}.{chrom}.wps.png")
        fig.savefig(out_png, dpi=120, bbox_inches="tight")
        plt.close(fig)
        saved.append(out_png)
        log.info("  ✓ %s → %s", chrom, out_png)

    log.info("chrom별 plot 완료: %d개 → %s", len(saved), output_dir)
    return saved


# ─────────────────────────────────────────────────────────────────────
# 2. Marker aggregate plot
# ─────────────────────────────────────────────────────────────────────
def plot_marker_aggregate(
    profiles_npy: str,
    output_path:  str,
    sample_id:    str = "",
    title_suffix: str = "",
) -> None:
    """
    aggregate.run() 출력 _profiles.npy 기반 6-panel aggregate plot

    패널 순서 (위→아래):
      L raw WPS / L adjusted WPS / L fragment coverage
      S raw WPS / S adjusted WPS / S fragment coverage

    Parameters
    ----------
    profiles_npy  : aggregate._profiles.npy 경로
    output_path   : PNG 저장 경로
    sample_id     : 타이틀 표시용
    title_suffix  : 추가 타이틀 문자열 (marker_type 등)
    """
    if not os.path.isfile(profiles_npy):
        log.warning("profiles.npy 없음: %s", profiles_npy)
        return

    data = np.load(profiles_npy, allow_pickle=True).item()
    pos          = data.get("pos",          np.array([]))
    marker_count = data.get("marker_count", 0)

    panel_specs = [
        # (data_key,   counts_key,    ylabel,                       positive_only, color,   mode_label)
        ("L_raw",      "L_raw_n",      "L Raw WPS",                  False, _TRACK_COLOR["raw"],      "L"),
        ("L_norm",     "L_norm_n",     "L Adjusted WPS",             False, _TRACK_COLOR["norm"],     "L"),
        ("L_frag_cov", "L_frag_cov_n", "L Fragment coverage",        True,  _TRACK_COLOR["frag_cov"],"L"),
        ("S_raw",      "S_raw_n",      "S Raw WPS",                  False, _TRACK_COLOR["raw"],      "S"),
        ("S_norm",     "S_norm_n",     "S Adjusted WPS",             False, _TRACK_COLOR["norm"],     "S"),
        ("S_frag_cov", "S_frag_cov_n", "S Fragment coverage",        True,  _TRACK_COLOR["frag_cov"],"S"),
    ]

    fig, axes = plt.subplots(6, 1, figsize=(14, 13), sharex=True, facecolor="white")
    fig.subplots_adjust(hspace=0.10)

    title = f"{sample_id} — Marker aggregate WPS | n={marker_count:,}"
    if title_suffix:
        title += f" [{title_suffix}]"
    fig.suptitle(title, fontsize=11, fontweight="bold", y=0.99)

    for ax, (dk, ck, ylabel, positive_only, color, mode_label) in zip(axes, panel_specs):
        ax.set_facecolor("white")
        arr    = data.get(dk, np.array([]))
        n_arr  = data.get(ck, np.zeros(len(pos), dtype=np.int64))
        finite = np.isfinite(arr)

        if not finite.any() or len(pos) == 0:
            ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                    ha="center", va="center", fontsize=9, color="#9E9E9E")
            ax.set_ylabel(ylabel, fontsize=8)
            continue

        ax.plot(pos, arr, linewidth=1.0, color=color)
        ax.axvline(0, linewidth=0.9, linestyle="--", alpha=0.6, color="#F44336")
        ax.axhline(0, linewidth=0.5, linestyle="--", alpha=0.4, color="#9E9E9E")

        vals = arr[finite]
        vmin, vmax = float(vals.min()), float(vals.max())

        if positive_only:
            ax.set_ylim(0, max(vmax * 1.15, 1.0))
        else:
            pad = max((vmax - vmin) * 0.12, 1e-6)
            ax.set_ylim(vmin - pad, vmax + pad)

        n_valid = int(n_arr.max()) if n_arr.size else 0
        ax.set_ylabel(f"[{mode_label}] {ylabel}\n(n≤{n_valid})", fontsize=8)
        ax.grid(axis="y", linewidth=0.3, alpha=0.25)
        ax.tick_params(axis="y", labelsize=7)
        ax.spines[["top", "right"]].set_visible(False)

    axes[-1].set_xlabel("Distance from marker center (bp)", fontsize=9)
    axes[-1].xaxis.set_major_formatter(mticker.FuncFormatter(
        lambda x, _: f"{int(x):+d}" if x != 0 else "0"
    ))
    axes[-1].tick_params(axis="x", labelsize=7)

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info("aggregate plot 저장: %s", output_path)


def plot_marker_aggregate_from_parquet(
    parquet_path: str,
    output_path:  str,
    sample_id:    str = "",
    title_suffix: str = "",
) -> None:
    """
    parquet 직접 입력 버전 (profiles.npy 없을 때 fallback)
    """
    import pandas as pd

    df = pd.read_parquet(parquet_path, engine="pyarrow")
    pos          = df["relative_position"].values
    marker_count = int(df.get("L_raw_n", pd.Series([0])).max())

    # npy 포맷으로 변환 후 재활용
    data = {
        "pos":           pos,
        "marker_count":  marker_count,
        "L_raw":         df["L_raw_wps"].values.astype(np.float32),
        "L_norm":        df["L_adjusted_wps"].values.astype(np.float32),
        "L_frag_cov":    df["L_frag_coverage"].values.astype(np.float32),
        "S_raw":         df["S_raw_wps"].values.astype(np.float32),
        "S_norm":        df["S_adjusted_wps"].values.astype(np.float32),
        "S_frag_cov":    df["S_frag_coverage"].values.astype(np.float32),
        **{f"{k}_n": df[f"{k}_n"].values if f"{k}_n" in df.columns
           else np.zeros(len(pos), dtype=np.int64)
           for k in ("L_raw", "L_norm", "L_frag_cov", "S_raw", "S_norm", "S_frag_cov")},
    }

    # 임시 npy 저장 후 재활용
    tmp_npy = parquet_path.replace(".parquet", "_tmp_profiles.npy")
    np.save(tmp_npy, data, allow_pickle=True)
    plot_marker_aggregate(tmp_npy, output_path, sample_id, title_suffix)
    os.remove(tmp_npy)


# ─────────────────────────────────────────────────────────────────────
# 3. Marker group overlay plot
# ─────────────────────────────────────────────────────────────────────
# 팔레트: marker_type 수에 따라 자동 할당
_GROUP_PALETTE = [
    "#1E88E5", "#E53935", "#43A047", "#FB8C00",
    "#8E24AA", "#00ACC1", "#F4511E", "#6D4C41",
    "#039BE5", "#3949AB", "#00897B", "#C0CA33",
]


def _draw_group_panel(
    axes:          list,
    group_data:    dict,
    profile_flank: int,
    color:         str,
) -> None:
    """단일 그룹 데이터를 6개 패널에 그림 (재사용 헬퍼)"""
    panel_specs = [
        ("L_raw",      "L Raw WPS",          False),
        ("L_norm",     "L Adjusted WPS",     False),
        ("L_frag_cov", "L Fragment coverage", True),
        ("S_raw",      "S Raw WPS",          False),
        ("S_norm",     "S Adjusted WPS",     False),
        ("S_frag_cov", "S Fragment coverage", True),
    ]
    pos = group_data.get("pos", np.arange(-profile_flank, profile_flank + 1))
    n   = int(group_data.get("marker_count", 0))

    for ax, (track_key, ylabel, positive_only) in zip(axes, panel_specs):
        ax.set_facecolor("white")
        ax.axvline(0, linewidth=0.9, linestyle="--", color="#F44336", alpha=0.7, zorder=3)
        ax.axhline(0, linewidth=0.5, linestyle="--", color="#9E9E9E", alpha=0.4, zorder=1)

        arr    = group_data.get(track_key)
        label  = f"{group_data.get('_label', '')} (n={n:,})"

        if arr is not None:
            finite = np.isfinite(arr)
            if finite.any():
                ax.plot(pos, arr, linewidth=1.3, color=color,
                        label=label, alpha=0.9, zorder=2)
                vals = arr[finite]
                if positive_only:
                    ax.set_ylim(0, float(np.percentile(vals, 99)) * 1.2)
                else:
                    p1, p99 = np.percentile(vals, [0.5, 99.5])
                    pad = (p99 - p1) * 0.15 or 1.0
                    ax.set_ylim(p1 - pad, p99 + pad)
            else:
                ax.text(0.5, 0.5, "No finite data",
                        transform=ax.transAxes, ha="center", va="center",
                        fontsize=8, color="#9E9E9E")
        else:
            ax.text(0.5, 0.5, "No data",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=8, color="#9E9E9E")

        ax.set_ylabel(ylabel, fontsize=8)
        ax.tick_params(axis="y", labelsize=7)
        ax.legend(fontsize=7, loc="upper right", framealpha=0.7)
        ax.grid(axis="y", linewidth=0.3, alpha=0.25)
        ax.spines[["top", "right"]].set_visible(False)


def plot_group_overlay(
    group_dir:     str,
    output_dir:    str,
    sample_id:     str = "",
    profile_flank: int = 1000,
) -> list[str]:
    """
    by_group/*_profiles.npy 를 읽어 그룹별로 1개씩 PNG 저장.

    패널 (위→아래):
      L Raw WPS / L Adjusted WPS / L Fragment coverage
      S Raw WPS / S Adjusted WPS / S Fragment coverage

    Parameters
    ----------
    group_dir     : aggregate.run() 이 저장한 by_group/ 디렉토리
    output_dir    : PNG 저장 디렉토리 (그룹별 파일 생성)
    sample_id     : 타이틀 표기용
    profile_flank : ±flank bp

    Returns
    -------
    저장된 PNG 경로 리스트
    """
    import glob

    npy_files = sorted(glob.glob(os.path.join(group_dir, "*_profiles.npy")))
    if not npy_files:
        log.warning("by_group/ 에 npy 파일 없음: %s", group_dir)
        return []

    os.makedirs(output_dir, exist_ok=True)
    saved: list[str] = []

    for idx, f in enumerate(npy_files):
        try:
            d = np.load(f, allow_pickle=True).item()
        except Exception as e:
            log.warning("npy 로드 실패 %s: %s", f, e)
            continue

        mt    = os.path.basename(f).replace("_profiles.npy", "")
        d["_label"] = mt
        color = _GROUP_PALETTE[idx % len(_GROUP_PALETTE)]
        n     = int(d.get("marker_count", 0))

        fig, axes = plt.subplots(6, 1, figsize=(14, 14),
                                  sharex=True, facecolor="white")
        fig.subplots_adjust(hspace=0.08)
        fig.suptitle(
            f"{sample_id} — [{mt}]  WPS aggregate profile  "
            f"(n={n:,} markers, ±{profile_flank}bp)",
            fontsize=11, fontweight="bold", y=0.995,
        )

        _draw_group_panel(axes, d, profile_flank, color)

        # X축
        axes[-1].set_xlabel("Distance from marker center (bp)", fontsize=9)
        axes[-1].xaxis.set_major_formatter(
            mticker.FuncFormatter(lambda x, _: f"{int(x):+d}" if x != 0 else "0")
        )
        axes[-1].tick_params(axis="x", labelsize=7)
        axes[-1].set_xlim(-profile_flank, profile_flank)

        out_png = os.path.join(output_dir, f"{sample_id}.{mt}.wps_profile.png")
        fig.savefig(out_png, dpi=150, bbox_inches="tight")
        plt.close(fig)
        saved.append(out_png)
        log.info("group plot 저장: %s (n=%d)", out_png, n)

    log.info("총 %d개 group plot 완료 → %s", len(saved), output_dir)
    return saved