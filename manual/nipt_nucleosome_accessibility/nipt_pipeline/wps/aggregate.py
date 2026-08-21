"""
wps/aggregate.py
================
manifest NPY에서 marker BED 영역을 slice → aggregate WPS profile

출력
----
{out_path}                                marker_stats.parquet  (전체 aggregate)
{out_path/_profiles.npy}                  전체 aggregate npy
{out_dir}/by_group/{marker_type}_profiles.npy   그룹별 aggregate npy
"""
from __future__ import annotations

import logging
import os
from typing import Optional

import numpy as np
import pandas as pd

from .utils import STANDARD_CHROMS, load_track, track_path, read_manifest

log = logging.getLogger(__name__)

TRACK_KEYS = (
    "L_raw", "L_cov", "L_frag_cov", "L_norm",
    "S_raw", "S_cov", "S_frag_cov", "S_norm",
)


# ── BED 로더 ─────────────────────────────────────────────────────────
def load_marker_bed(bed_path: str) -> pd.DataFrame:
    df = pd.read_csv(
        bed_path, sep="\t", header=None, comment="#",
        names=["chrom", "start", "end", "marker_id", "marker_type"],
    )
    df["start"] = pd.to_numeric(df["start"], errors="raise").astype("int64")
    df["end"]   = pd.to_numeric(df["end"],   errors="raise").astype("int64")

    miss = df["marker_id"].isna()
    df.loc[miss, "marker_id"] = (
        df.loc[miss, "chrom"].astype(str) + ":" +
        df.loc[miss, "start"].astype(str) + "-" +
        df.loc[miss, "end"].astype(str)
    )
    df["marker_type"] = df["marker_type"].fillna("marker").astype(str)
    df["marker_id"]   = df["marker_id"].astype(str)
    df = df[df["end"] > df["start"]].copy()
    log.info("marker BED 로드: %d markers, %d types [%s]",
             len(df), df["marker_type"].nunique(), bed_path)
    return df


# ── chrom 트랙 로드 ───────────────────────────────────────────────────
def _load_chrom_tracks(
    wps_dir:   str,
    chrom:     str,
    prefix:    str,
    manifests: dict[str, dict],
) -> dict[str, Optional[np.ndarray]]:
    tracks: dict[str, Optional[np.ndarray]] = {}
    for key in TRACK_KEYS:
        mode, metric = key.split("_", 1)
        tracks[key] = load_track(
            wps_dir  = wps_dir,
            mode     = mode,
            chrom    = chrom,
            metric   = metric,
            prefix   = prefix,
            manifest = manifests.get(mode),
        )
    return tracks


# ── 누적기 (sums / counts 딕셔너리 쌍) ────────────────────────────────
def _make_accum(profile_len: int) -> tuple[dict, dict]:
    sums   = {k: np.zeros(profile_len, dtype=np.float64) for k in TRACK_KEYS}
    counts = {k: np.zeros(profile_len, dtype=np.int64)   for k in TRACK_KEYS}
    return sums, counts


def _accumulate(
    sums:          dict[str, np.ndarray],
    counts:        dict[str, np.ndarray],
    tracks:        dict[str, Optional[np.ndarray]],
    center:        int,
    chrom_len:     int,
    profile_flank: int,
    profile_len:   int,
) -> bool:
    if center < 0 or center >= chrom_len:
        return False

    sl_start = max(0,         center - profile_flank)
    sl_end   = min(chrom_len, center + profile_flank + 1)
    if sl_end <= sl_start:
        return False

    ag_start = (sl_start - center) + profile_flank
    ag_end   = (sl_end   - center) + profile_flank

    any_valid = False
    for key in TRACK_KEYS:
        arr = tracks.get(key)
        if arr is None:
            continue
        seg    = np.asarray(arr[sl_start:sl_end], dtype=np.float64)
        finite = np.isfinite(seg)
        if not finite.any():
            continue
        sums[key][ag_start:ag_end][finite]   += seg[finite]
        counts[key][ag_start:ag_end][finite] += 1
        any_valid = True

    return any_valid


def _baseline_mean(arr: np.ndarray, profile_flank: int, baseline_pct: float = 0.2) -> float:
    """
    profile 양 끝 baseline_pct 구간의 평균으로 baseline 추정.

    center 근처 신호에 영향받지 않고 flank 끝 구간만 사용.
    예) flank=1000, baseline_pct=0.2 → 양쪽 끝 200bp 구간 평균

    TSS/CTCF 같이 center에 강한 신호가 있는 경우
    전체 trimmed mean을 쓰면 신호 방향이 반전되는 문제 방지.
    """
    n         = len(arr)
    flank_n   = int(n * baseline_pct / 2)  # 양쪽 각각
    flank_n   = max(flank_n, 10)

    left_end  = arr[:flank_n]
    right_end = arr[n - flank_n:]
    baseline  = np.concatenate([left_end, right_end])
    vals      = baseline[np.isfinite(baseline)]

    if vals.size == 0:
        return 0.0

    # 극단값 제거 (25~75% 내)
    lo  = np.percentile(vals, 25)
    hi  = np.percentile(vals, 75)
    mid = vals[(vals >= lo) & (vals <= hi)]
    return float(np.mean(mid)) if mid.size > 0 else float(np.mean(vals))


def _finalize(
    sums:              dict[str, np.ndarray],
    counts:            dict[str, np.ndarray],
    profile_len:       int,
    marker_count:      int,
    profile_flank:     int,
    relative_position: np.ndarray,
    baseline_pct:      float = 0.2,
) -> dict:
    """
    sums/counts → nanmean → baseline subtraction → save_dict

    baseline: profile 양 끝 20% 구간 평균 (flank 끝 = signal 없는 구간)
    → center 신호에 영향받지 않아 TSS/CTCF 패턴이 올바른 방향으로 표시됨
    → WPS / norm 트랙만 적용, coverage 트랙은 원본 유지
    """
    aggregate: dict[str, np.ndarray] = {}
    for key in TRACK_KEYS:
        arr   = np.full(profile_len, np.nan, dtype=np.float32)
        valid = counts[key] > 0
        arr[valid] = (sums[key][valid] / counts[key][valid]).astype(np.float32)

        # baseline subtraction: WPS / norm 트랙만 (coverage는 유지)
        if not key.endswith("_cov") and not key.endswith("_frag_cov"):
            baseline = _baseline_mean(arr, profile_flank, baseline_pct)
            arr = np.where(np.isfinite(arr), arr - baseline, np.nan).astype(np.float32)

        aggregate[key] = arr

    return {
        "marker_count":  marker_count,
        "profile_flank": profile_flank,
        "pos":           relative_position,
        **{k: aggregate[k]   for k in TRACK_KEYS},
        **{f"{k}_n": counts[k] for k in TRACK_KEYS},
    }


# ── 공개 인터페이스 ────────────────────────────────────────────────────
def run(
    marker_bed:    str,
    out_path:      str,
    wps_dir:       str,
    profile_flank: int  = 1000,
    prefix:        str  = "",
    save_profiles: bool = True,
) -> pd.DataFrame:
    """
    Parameters
    ----------
    marker_bed    : BED (chrom start end marker_id marker_type)
    out_path      : marker_stats.parquet 저장 경로
    wps_dir       : compute.run() out_dir 루트 (L/ S/ 상위)
    profile_flank : marker center ±flank bp  (default 1000 → 2001bp 프로파일)
    prefix        : compute prefix (sample_id)
    save_profiles : npy 저장 여부

    출력
    ----
    {out_path}                           전체 aggregate parquet
    {out_dir}/_profiles.npy              전체 aggregate npy
    {out_dir}/by_group/{type}_profiles.npy  marker_type별 aggregate npy
    """
    marker_df = load_marker_bed(marker_bed)
    if marker_df.empty:
        log.warning("marker 없음 — aggregate 건너뜀")
        return pd.DataFrame()

    wps_dir  = os.path.abspath(wps_dir)
    out_path = os.path.abspath(out_path)
    out_dir  = os.path.dirname(out_path)
    os.makedirs(out_dir, exist_ok=True)

    manifests: dict[str, dict] = {
        mode: read_manifest(wps_dir, mode) for mode in ("L", "S")
    }

    profile_len       = profile_flank * 2 + 1
    relative_position = np.arange(-profile_flank, profile_flank + 1, dtype=np.int32)

    # 전체 누적기
    all_sums, all_counts = _make_accum(profile_len)
    all_total = 0

    # marker_type별 누적기
    mtypes      = marker_df["marker_type"].unique().tolist()
    group_sums  = {mt: _make_accum(profile_len) for mt in mtypes}
    group_total = {mt: 0 for mt in mtypes}

    # marker를 STANDARD_CHROMS 순서로 정렬
    # → 동일 chrom NPY를 한 번만 로드하고 순서대로 처리
    chrom_order = {c: i for i, c in enumerate(STANDARD_CHROMS)}
    marker_df = marker_df.copy()
    marker_df["_chrom_order"] = marker_df["chrom"].map(
        lambda c: chrom_order.get(c, 9999)
    )
    marker_df = marker_df.sort_values(["_chrom_order", "start"]).drop(
        columns="_chrom_order"
    )

    # chrom별 순서대로 NPY 1회 로드 후 처리
    for chrom, cdf in marker_df.groupby("chrom", sort=False):
        chrom = str(chrom)
        log.info("[aggregate] chrom=%s markers=%d", chrom, len(cdf))

        tracks   = _load_chrom_tracks(wps_dir, chrom, prefix, manifests)
        valid_tr = [a for a in tracks.values() if a is not None]

        if not valid_tr:
            log.warning("  ! %s: 모든 track 없음 — 건너뜀", chrom)
            for key in TRACK_KEYS:
                mode, metric = key.split("_", 1)
                p = track_path(wps_dir, mode, chrom, metric, prefix)
                log.debug("    tried: %s (exists=%s)", p, os.path.isfile(p))
            continue

        chrom_len = min(int(a.size) for a in valid_tr)
        n_chrom   = 0

        for row in cdf.itertuples(index=False):
            center = (int(row.start) + int(row.end)) // 2
            mtype  = str(row.marker_type)

            # 전체 누적
            hit = _accumulate(
                all_sums, all_counts, tracks,
                center, chrom_len, profile_flank, profile_len,
            )
            if hit:
                all_total += 1
                n_chrom   += 1
                # 그룹별 누적
                gs, gc = group_sums[mtype]
                _accumulate(gs, gc, tracks, center, chrom_len, profile_flank, profile_len)
                group_total[mtype] += 1

        log.info("  ✓ %s aggregated=%d", chrom, n_chrom)
        tracks.clear()

    if all_total == 0:
        log.warning("aggregate 가능한 marker 없음")
        return pd.DataFrame()

    # ── 전체 aggregate 저장 ───────────────────────────────────────────
    all_dict = _finalize(all_sums, all_counts, profile_len, all_total,
                         profile_flank, relative_position)

    # parquet
    result_df = pd.DataFrame({
        "relative_position": relative_position,
        "L_raw_wps":       all_dict["L_raw"],
        "L_coverage":      all_dict["L_cov"],
        "L_frag_coverage": all_dict["L_frag_cov"],
        "L_adjusted_wps":  all_dict["L_norm"],
        "S_raw_wps":       all_dict["S_raw"],
        "S_coverage":      all_dict["S_cov"],
        "S_frag_coverage": all_dict["S_frag_cov"],
        "S_adjusted_wps":  all_dict["S_norm"],
        **{f"{k}_n": all_dict[f"{k}_n"] for k in TRACK_KEYS},
    })
    result_df.to_parquet(out_path, engine="pyarrow", index=False)
    log.info("전체 marker_stats 저장: %s (markers=%d)", out_path, all_total)

    if save_profiles:
        npy_path = out_path.replace(".parquet", "_profiles.npy")
        np.save(npy_path, all_dict, allow_pickle=True)
        log.info("전체 profiles.npy 저장: %s", npy_path)

    # ── marker_type별 aggregate 저장 ──────────────────────────────────
    if save_profiles:
        group_dir = os.path.join(out_dir, "by_group")
        os.makedirs(group_dir, exist_ok=True)

        for mt in mtypes:
            n = group_total[mt]
            if n == 0:
                log.warning("marker_type=%s aggregated=0 — 저장 스킵", mt)
                continue
            gs, gc = group_sums[mt]
            g_dict = _finalize(gs, gc, profile_len, n, profile_flank, relative_position)
            # 파일명에 사용 불가 문자 치환
            safe_mt  = "".join(c if c.isalnum() or c in "-_." else "_" for c in mt)
            npy_path = os.path.join(group_dir, f"{safe_mt}_profiles.npy")
            np.save(npy_path, g_dict, allow_pickle=True)
            log.info("  group [%s] markers=%d → %s", mt, n, npy_path)

    return result_df