"""
marker_extractor.py
===================
wps_compute 출력 NPY에서 marker 영역을 slice하여
L/S-WPS aggregate profile을 생성합니다.

입력 WPS 구조 (wps_compute 출력)
---------------------------------
{wps_dir}/L/chr1.raw.npy
{wps_dir}/L/chr1.cov.npy
{wps_dir}/L/chr1.norm.npy
{wps_dir}/S/chr1.raw.npy
{wps_dir}/S/chr1.cov.npy
{wps_dir}/S/chr1.norm.npy

출력
----
marker_stats.parquet          : position별 aggregate (relative_position, L/S track)
marker_stats_profiles.npy     : aggregate array dict (pipeline 호환)
{plot_dir}/aggregate_marker_wps.png
"""
from __future__ import annotations

import logging
import os
import re
from typing import Optional

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

BUILD_TAG = "marker-aggregate-20260820-v1"


# ─────────────────────────────────────────────────────────────────────
# BED 로더
# ─────────────────────────────────────────────────────────────────────
def load_marker_bed(bed_path: str) -> pd.DataFrame:
    df = pd.read_csv(
        bed_path, sep="\t", header=None, comment="#",
        names=["chrom", "start", "end", "marker_id", "marker_type"],
    )
    df["start"] = pd.to_numeric(df["start"], errors="raise").astype("int64")
    df["end"]   = pd.to_numeric(df["end"],   errors="raise").astype("int64")

    if df["marker_id"].isna().all():
        df["marker_id"] = (df["chrom"].astype(str) + ":" +
                           df["start"].astype(str) + "-" + df["end"].astype(str))
    else:
        miss = df["marker_id"].isna()
        df.loc[miss, "marker_id"] = (
            df.loc[miss, "chrom"].astype(str) + ":" +
            df.loc[miss, "start"].astype(str) + "-" + df.loc[miss, "end"].astype(str)
        )

    df["marker_type"] = df["marker_type"].fillna("marker").astype(str)
    df["marker_id"]   = df["marker_id"].astype(str)
    df = df[df["end"] > df["start"]].copy()
    log.info("marker BED: %d markers [%s]", len(df), bed_path)
    return df


# ─────────────────────────────────────────────────────────────────────
# track 로더
# ─────────────────────────────────────────────────────────────────────
def _load_track(
    wps_dir: str, mode: str, chrom: str, metric: str,
    prefix: str = "",
) -> Optional[np.ndarray]:
    """
    {wps_dir}/{mode}/{prefix}.{chrom}.{metric}.npy 로드.
    prefix 없으면 {chrom}.{metric}.npy.
    없으면 None 반환.
    """
    suffix     = {"raw": "raw.npy", "cov": "cov.npy", "norm": "norm.npy"}[metric]
    fname_base = f"{prefix}.{chrom}" if prefix else chrom
    path       = os.path.join(wps_dir, mode, f"{fname_base}.{suffix}")
    if not os.path.isfile(path):
        return None
    arr = np.load(path, mmap_mode="r", allow_pickle=False)
    if arr.ndim != 1:
        raise ValueError(f"WPS NPY must be 1-D: {path}")
    return arr


# ─────────────────────────────────────────────────────────────────────
# aggregate plot
# ─────────────────────────────────────────────────────────────────────
def _plot_aggregate(
    relative_position: np.ndarray,
    aggregate: dict[str, np.ndarray],
    counts: dict[str, np.ndarray],
    marker_count: int,
    out_png: str,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    panel_order = [
        ("L_cov",  "L coverage (120-180bp)",  True,  "#2a78d6"),
        ("L_raw",  "L raw WPS",               False, "#2a78d6"),
        ("L_norm", "L adjusted WPS",           False, "#2a78d6"),
        ("S_cov",  "S coverage (35-80bp)",     True,  "#eb6834"),
        ("S_raw",  "S raw WPS",                False, "#eb6834"),
        ("S_norm", "S adjusted WPS",           False, "#eb6834"),
    ]

    os.makedirs(os.path.dirname(os.path.abspath(out_png)), exist_ok=True)
    fig, axes = plt.subplots(6, 1, figsize=(14, 12), sharex=True, facecolor="white")

    for ax, (key, ylabel, positive_only, color) in zip(axes, panel_order):
        arr    = aggregate.get(key, np.array([]))
        finite = np.isfinite(arr)

        if not finite.any():
            ax.text(0.5, 0.5, "No data",
                    transform=ax.transAxes, ha="center", va="center")
            ax.set_ylabel(ylabel, fontsize=9)
            continue

        ax.plot(relative_position, arr, linewidth=1.0, color=color)
        ax.axvline(0, linewidth=0.8, linestyle="--", alpha=0.6, color="black")

        vals = arr[finite]
        vmin, vmax = float(vals.min()), float(vals.max())

        if positive_only:
            ax.set_ylim(0, max(vmax * 1.1, 1.0))
        else:
            pad = max((vmax - vmin) * 0.1, 0.01) if not np.isclose(vmin, vmax) else 0.01
            ax.set_ylim(vmin - pad, vmax + pad)
            ax.axhline(0, linewidth=0.5, linestyle="--", alpha=0.4, color="gray")

        n_valid = int(counts.get(key, np.array([0])).max())
        ax.set_ylabel(f"{ylabel}\n(n≤{n_valid})", fontsize=9)
        ax.grid(axis="y", linewidth=0.3, alpha=0.25)
        ax.set_facecolor("white")
        for sp in ax.spines.values(): sp.set_linewidth(0.4)

    axes[-1].set_xlabel("Distance from marker center (bp)", fontsize=10)
    fig.suptitle(f"Aggregate marker-centered WPS | total markers={marker_count:,}",
                 fontsize=13)
    plt.tight_layout()
    plt.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)
    log.info("aggregate plot 저장: %s", out_png)


# ─────────────────────────────────────────────────────────────────────
# 공개 인터페이스
# ─────────────────────────────────────────────────────────────────────
def run(
    marker_bed:    str,
    out_path:      str,
    wps_dir:       str,
    bam_path:      Optional[str] = None,   # 하위 호환 (사용 안 함)
    fasta_path:    Optional[str] = None,   # 하위 호환 (사용 안 함)
    bw_path:       Optional[str] = None,   # 하위 호환 (사용 안 함)
    min_mapq:      int  = 20,              # 하위 호환 (사용 안 함)
    n_jobs:        int  = 1,               # 하위 호환 (사용 안 함)
    save_profiles: bool = True,
    make_plots:    bool = True,
    plot_dir:      Optional[str] = None,
    profile_flank: int  = 1000,
    max_plot_markers: Optional[int] = None,  # 하위 호환 (사용 안 함)
    prefix:        str  = "",              # wps_compute prefix (sample_id)
) -> pd.DataFrame:
    """
    wps_compute NPY에서 marker 영역을 slice하여 L/S-WPS aggregate profile을 생성합니다.

    Parameters
    ----------
    marker_bed    : BED 파일 (chrom start end marker_id marker_type)
    out_path      : marker_stats.parquet 저장 경로
    wps_dir       : wps_compute out_dir (L/, S/ 하위 디렉토리 포함)
    profile_flank : marker center 기준 ±flank bp
    prefix        : sample_id prefix (wps_compute의 prefix와 동일해야 함)
                    manifest.json 에서 실제 경로를 읽어 prefix를 자동 감지합니다.
    """
    marker_df = load_marker_bed(marker_bed)
    if marker_df.empty:
        log.warning("marker 없음")
        return pd.DataFrame()

    wps_dir  = os.path.abspath(wps_dir)
    out_path = os.path.abspath(out_path)

    if plot_dir is None:
        plot_dir = os.path.join(os.path.dirname(out_path), "marker_plots")

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    if make_plots:
        os.makedirs(plot_dir, exist_ok=True)

    # manifest.json에서 실제 파일 경로 로드 → prefix 자동 감지
    # {wps_dir}/{mode}/manifest.json → chromosomes: {chrom: {raw, cov, norm}}
    import json as _json
    real_paths: dict[str, dict[str, dict[str, str]]] = {}  # {mode: {chrom: {metric: path}}}
    for mode in ("L", "S"):
        manifest_path = os.path.join(wps_dir, mode, "manifest.json")
        if not os.path.isfile(manifest_path):
            log.warning("manifest 없음: %s", manifest_path)
            continue
        with open(manifest_path) as _f:
            manifest = _json.load(_f)
        real_paths[mode] = {}
        for chrom, entry in manifest.get("chromosomes", {}).items():
            if isinstance(entry, dict):
                real_paths[mode][chrom] = {
                    m: entry.get(m, "") for m in ("raw", "cov", "norm")
                }
            elif isinstance(entry, str):
                real_paths[mode][chrom] = {"norm": entry, "raw": "", "cov": ""}
        log.info("manifest 로드: mode=%s chroms=%d", mode, len(real_paths[mode]))

    log.info("marker build=%s flank=±%d bp total=%d",
             BUILD_TAG, profile_flank, len(marker_df))

    profile_len       = profile_flank * 2 + 1
    relative_position = np.arange(-profile_flank, profile_flank + 1, dtype=np.int32)

    TRACK_KEYS = ("L_cov", "L_raw", "L_norm", "S_cov", "S_raw", "S_norm")
    sums   = {k: np.zeros(profile_len, dtype=np.float64) for k in TRACK_KEYS}
    counts = {k: np.zeros(profile_len, dtype=np.int64)   for k in TRACK_KEYS}

    marker_count = 0

    for chrom, cdf in marker_df.groupby("chrom", sort=False):
        chrom = str(chrom)
        log.info("[Marker] chrom=%s markers=%d", chrom, len(cdf))

        # 트랙 로드 (manifest 실제 경로 우선, fallback: prefix 기반 경로)
        tracks: dict[str, Optional[np.ndarray]] = {}
        for mode in ("L", "S"):
            for metric in ("raw", "cov", "norm"):
                key  = f"{mode}_{metric}"
                path = ""
                # manifest에서 실제 경로 확인
                if mode in real_paths and chrom in real_paths[mode]:
                    path = real_paths[mode][chrom].get(metric, "")
                # fallback: prefix 기반 경로 구성
                if not path or not os.path.isfile(path):
                    arr = _load_track(wps_dir, mode, chrom, metric, prefix=prefix)
                    tracks[key] = arr
                    if arr is None:
                        log.warning("  ! track 없음: %s/%s/%s.%s.npy",
                                    wps_dir, mode, chrom, metric)
                    continue
                try:
                    arr = np.load(path, mmap_mode="r", allow_pickle=False)
                    if arr.ndim != 1:
                        raise ValueError(f"1-D expected: {path}")
                    tracks[key] = arr
                except Exception as e:
                    log.warning("  ! track 로드 실패 %s: %s", path, e)
                    tracks[key] = None

        valid = [a for a in tracks.values() if a is not None]
        if not valid:
            log.warning("  ! %s: 모든 track 없음 — 건너뜀", chrom)
            continue

        chrom_len = min(int(a.size) for a in valid)
        n_chrom   = 0

        for r in cdf.itertuples(index=False):
            center = (int(r.start) + int(r.end)) // 2
            if center < 0 or center >= chrom_len:
                continue

            sl_start = max(0,          center - profile_flank)
            sl_end   = min(chrom_len,  center + profile_flank + 1)
            if sl_end <= sl_start:
                continue

            ag_start = (sl_start - center) + profile_flank
            ag_end   = (sl_end   - center) + profile_flank

            for key in TRACK_KEYS:
                arr = tracks.get(key)
                if arr is None:
                    continue
                seg    = np.asarray(arr[sl_start:sl_end], dtype=np.float32)
                finite = np.isfinite(seg)
                if not finite.any():
                    continue
                sums[key][ag_start:ag_end][finite]   += seg[finite]
                counts[key][ag_start:ag_end][finite] += 1

            marker_count += 1
            n_chrom      += 1

        log.info("  ✓ %s aggregated=%d", chrom, n_chrom)
        tracks.clear()

    if marker_count == 0:
        log.warning("aggregate 가능한 marker 없음")
        return pd.DataFrame()

    # aggregate 계산
    aggregate: dict[str, np.ndarray] = {}
    for key in TRACK_KEYS:
        arr   = np.full(profile_len, np.nan, dtype=np.float32)
        valid = counts[key] > 0
        arr[valid] = (sums[key][valid] / counts[key][valid]).astype(np.float32)
        aggregate[key] = arr
        if valid.any():
            log.info("[AGGREGATE] key=%s markers=%d valid=%d min=%.4f max=%.4f",
                     key, marker_count, int(valid.sum()),
                     float(np.nanmin(arr)), float(np.nanmax(arr)))
        else:
            log.warning("[AGGREGATE] key=%s valid position 없음", key)

    # marker_stats.parquet 저장
    result_df = pd.DataFrame({
        "relative_position": relative_position,
        "L_coverage":   aggregate["L_cov"],
        "L_raw_wps":    aggregate["L_raw"],
        "L_adjusted_wps": aggregate["L_norm"],
        "S_coverage":   aggregate["S_cov"],
        "S_raw_wps":    aggregate["S_raw"],
        "S_adjusted_wps": aggregate["S_norm"],
        "L_cov_n":   counts["L_cov"],
        "L_raw_n":   counts["L_raw"],
        "L_norm_n":  counts["L_norm"],
        "S_cov_n":   counts["S_cov"],
        "S_raw_n":   counts["S_raw"],
        "S_norm_n":  counts["S_norm"],
    })
    result_df.to_parquet(out_path, engine="pyarrow", index=False)
    log.info("marker_stats 저장: %s (positions=%d, markers=%d)",
             out_path, len(result_df), marker_count)

    # profiles.npy 저장 — pipeline._plot_marker_profiles 호환 구조
    if save_profiles:
        npy_path = out_path.replace(".parquet", "_profiles.npy")
        # pipeline은 {marker_id: {...}} 구조를 기대하므로
        # aggregate 결과를 단일 "aggregate" 키로 저장
        # _plot_marker_profiles는 marker_type별 iterate하므로
        # type별 dummy key로 래핑
        profiles_dict = {
            "__aggregate__": {
                "pos":        relative_position,
                "short_wps_L": aggregate["L_norm"],
                "long_wps_L":  aggregate["L_norm"],
                "short_wps_S": aggregate["S_norm"],
                "long_wps_S":  aggregate["S_norm"],
            }
        }
        # 실제 aggregate dict도 함께 저장
        save_dict = {
            "marker_count":  marker_count,
            "profile_flank": profile_flank,
            "pos":           relative_position,
            "L_cov":  aggregate["L_cov"],
            "L_raw":  aggregate["L_raw"],
            "L_norm": aggregate["L_norm"],
            "S_cov":  aggregate["S_cov"],
            "S_raw":  aggregate["S_raw"],
            "S_norm": aggregate["S_norm"],
            "L_cov_n":  counts["L_cov"],
            "L_raw_n":  counts["L_raw"],
            "L_norm_n": counts["L_norm"],
            "S_cov_n":  counts["S_cov"],
            "S_raw_n":  counts["S_raw"],
            "S_norm_n": counts["S_norm"],
            # pipeline 호환 profiles 키
            "profiles": profiles_dict,
        }
        np.save(npy_path, save_dict, allow_pickle=True)
        log.info("profiles.npy 저장: %s", npy_path)

    # aggregate plot (전체)
    if make_plots:
        agg_png = os.path.join(plot_dir, "aggregate_marker_wps.png")
        _plot_aggregate(relative_position, aggregate, counts, marker_count, agg_png)

        # marker_type별 개별 aggregate plot
        for mtype in marker_df["marker_type"].unique():
            mtype_df = marker_df[marker_df["marker_type"] == mtype]
            type_sums   = {k: np.zeros(profile_len, dtype=np.float64) for k in TRACK_KEYS}
            type_counts = {k: np.zeros(profile_len, dtype=np.int64)   for k in TRACK_KEYS}
            type_count  = 0

            for chrom, cdf in mtype_df.groupby("chrom", sort=False):
                chrom = str(chrom)
                type_tracks: dict[str, Optional[np.ndarray]] = {}
                for mode in ("L", "S"):
                    for metric in ("raw", "cov", "norm"):
                        key = f"{mode}_{metric}"
                        # manifest 실제 경로 우선
                        _tp = ""
                        if mode in real_paths and chrom in real_paths[mode]:
                            _tp = real_paths[mode][chrom].get(metric, "")
                        if _tp and os.path.isfile(_tp):
                            try:
                                _arr = np.load(_tp, mmap_mode="r", allow_pickle=False)
                                type_tracks[key] = _arr if _arr.ndim == 1 else None
                            except Exception:
                                type_tracks[key] = None
                        else:
                            type_tracks[key] = _load_track(wps_dir, mode, chrom, metric, prefix=prefix)

                valid = [a for a in type_tracks.values() if a is not None]
                if not valid:
                    continue
                chrom_len = min(int(a.size) for a in valid)

                for r in cdf.itertuples(index=False):
                    center   = (int(r.start) + int(r.end)) // 2
                    if center < 0 or center >= chrom_len:
                        continue
                    sl_start = max(0, center - profile_flank)
                    sl_end   = min(chrom_len, center + profile_flank + 1)
                    ag_start = (sl_start - center) + profile_flank
                    ag_end   = (sl_end   - center) + profile_flank
                    for key in TRACK_KEYS:
                        arr = type_tracks.get(key)
                        if arr is None:
                            continue
                        seg    = np.asarray(arr[sl_start:sl_end], dtype=np.float32)
                        finite = np.isfinite(seg)
                        if not finite.any():
                            continue
                        type_sums[key][ag_start:ag_end][finite]   += seg[finite]
                        type_counts[key][ag_start:ag_end][finite] += 1
                    type_count += 1

            if type_count == 0:
                continue

            type_agg = {}
            for key in TRACK_KEYS:
                arr   = np.full(profile_len, np.nan, dtype=np.float32)
                valid = type_counts[key] > 0
                arr[valid] = (type_sums[key][valid] / type_counts[key][valid]).astype(np.float32)
                type_agg[key] = arr

            safe_type = re.sub(r"[^A-Za-z0-9._-]+", "_", str(mtype)).strip("._") or "marker"
            type_png  = os.path.join(plot_dir, f"{safe_type}_marker_wps.png")
            _plot_aggregate(relative_position, type_agg, type_counts, type_count, type_png)
            log.info("  marker_type=%s n=%d -> %s", mtype, type_count, type_png)

    return result_df