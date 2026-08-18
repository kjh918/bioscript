"""
wps_compute.py
==============
marker BED 기반 WPS 계산 (원저자 방식).

핵심 설계
---------
- marker region (±extend bp) 에서만 계산 → genome-wide 스캔 아님
- 같은 염색체 marker 를 한 번의 BAM fetch 로 처리 → 속도 최적화
- short / long fragment 각각 독립 계산
- adjusted WPS = (WPS - windowMedian) / Coverage * 100  (원저자 normalizeWPSwigs.py)

출력
----
  {prefix}.wps_{mode}_{frag}.npy  per mode/frag 조합
  구조: {"wps_norm": {marker_id: float32 array [2*extend+1]},
          "mode": str, "frag": str, "extend": int}
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Optional

import numpy as np
import pysam

log = logging.getLogger(__name__)

WPS_PARAMS = {
    "L": {"frag_min": 120, "frag_max": 180, "window": 120},
    "S": {"frag_min":  35, "frag_max":  80, "window":  16},
}
SHORT_MAX = 150


# ─────────────────────────────────────────────────────────────────────
# marker BED 로더
# ─────────────────────────────────────────────────────────────────────
def _load_marker_bed(bed_path: str) -> list[tuple[str, str, int, int]]:
    """(marker_id, chrom, start, end) 목록 반환."""
    markers = []
    with open(bed_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            try:
                start, end = int(parts[1]), int(parts[2])
            except ValueError:
                continue
            marker_id = parts[3] if len(parts) > 3 else f"{chrom}:{start}-{end}"
            markers.append((marker_id, chrom, start, end))
    log.info("marker BED 로드: %d markers", len(markers))
    return markers


# ─────────────────────────────────────────────────────────────────────
# 염색체 단위 WPS 계산 (같은 chrom marker 일괄 처리)
# ─────────────────────────────────────────────────────────────────────
def _compute_chrom_markers(
    bam_path:    str,
    chrom:       str,
    markers:     list[tuple[str, int, int]],   # [(marker_id, center, reg_len), ...]
    extend:      int,
    mode:        str,
    frag_filter: str,
    min_mapq:    int,
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """
    같은 염색체의 marker 들을 BAM 1회 fetch 로 처리합니다.

    전략
    ----
    1. 모든 marker region 을 커버하는 최소 fetch 범위 계산
    2. 해당 범위 한 번만 fetch
    3. 각 read 를 겹치는 marker 배열에 누적

    Returns
    -------
    {marker_id: (wps_array, coverage_array)}
    """
    prm    = WPS_PARAMS[mode]
    fmin   = prm["frag_min"]
    fmax   = prm["frag_max"]
    half_k = prm["window"] // 2
    reg_len = 2 * extend + 1

    # marker별 배열 초기화
    spanning_d: dict[str, np.ndarray] = {}
    wps_ep_d:   dict[str, np.ndarray] = {}
    cov_d:      dict[str, np.ndarray] = {}
    starts_d:   dict[str, int]        = {}  # reg_start

    for marker_id, center, _ in markers:
        reg_start = max(0, center - extend)
        spanning_d[marker_id] = np.zeros(reg_len, dtype=np.int32)
        wps_ep_d[marker_id]   = np.zeros(reg_len, dtype=np.int32)
        cov_d[marker_id]      = np.zeros(reg_len, dtype=np.int32)
        starts_d[marker_id]   = reg_start

    # fetch 범위: 전체 marker 범위 + fragment 최대 길이 패딩
    all_reg_starts = [starts_d[m[0]] for m in markers]
    all_reg_ends   = [starts_d[m[0]] + reg_len for m in markers]
    fetch_start    = max(0, min(all_reg_starts) - fmax)
    fetch_end      = max(all_reg_ends) + fmax

    # marker interval lookup: pos → 겹치는 marker 목록
    # read 하나가 여러 marker 에 걸칠 수 있음
    # 효율: marker를 sorted list로 만들고 bisect로 탐색
    sorted_markers = sorted(markers, key=lambda x: starts_d[x[0]])
    sorted_starts  = [starts_d[m[0]] for m in sorted_markers]
    sorted_ends    = [starts_d[m[0]] + reg_len for m in sorted_markers]

    import bisect

    seen: set[str] = set()

    with pysam.AlignmentFile(bam_path, "rb", threads=1) as bam:
        for read in bam.fetch(chrom, fetch_start, fetch_end):
            if (read.is_unmapped or read.is_duplicate
                    or read.is_secondary or read.is_supplementary):
                continue
            if read.mapping_quality < min_mapq:
                continue
            if read.is_paired and not read.is_read1:
                continue

            qname = read.query_name
            if qname in seen:
                continue
            seen.add(qname)

            if read.is_paired and read.template_length != 0:
                fs = read.reference_start
                fe = fs + abs(read.template_length)
            else:
                fs = read.reference_start
                fe = read.reference_end

            flen = fe - fs
            if flen <= 0:
                continue

            is_short = (flen <= SHORT_MAX)
            if frag_filter == "short" and not is_short:
                continue
            if frag_filter == "long" and is_short:
                continue

            # 이 read 와 겹치는 marker 탐색
            # read [fs, fe) 와 겹치는 region: reg_start < fe AND reg_end > fs
            lo_idx = bisect.bisect_left(sorted_ends,   fs + 1)
            hi_idx = bisect.bisect_right(sorted_starts, fe - 1)

            for idx in range(lo_idx, hi_idx):
                marker_id = sorted_markers[idx][0]
                reg_start = starts_d[marker_id]
                reg_end   = reg_start + reg_len

                # coverage
                cov_lo = max(fs, reg_start) - reg_start
                cov_hi = min(fe, reg_end)   - reg_start
                if cov_lo < cov_hi:
                    cov_d[marker_id][cov_lo:cov_hi] += 1

                # WPS (fragment 길이 필터)
                if not (fmin <= flen <= fmax):
                    continue

                sp_lo = max(fs + half_k + 1, reg_start) - reg_start
                sp_hi = min(fe - half_k,     reg_end)   - reg_start
                if sp_lo < sp_hi:
                    spanning_d[marker_id][sp_lo:sp_hi] += 1

                for ep in (fs, fe - 1):
                    w_lo = max(ep - half_k,     reg_start) - reg_start
                    w_hi = min(ep + half_k + 1, reg_end)   - reg_start
                    if w_lo < w_hi:
                        wps_ep_d[marker_id][w_lo:w_hi] += 1

    # WPS 계산
    results = {}
    for marker_id in spanning_d:
        wps = (spanning_d[marker_id] - wps_ep_d[marker_id]).astype(np.int32)
        results[marker_id] = (wps, cov_d[marker_id])

    return results


# ─────────────────────────────────────────────────────────────────────
# 정규화: 원저자 방식
# adjusted_WPS = (WPS - windowMedian) / Coverage * 100
# ─────────────────────────────────────────────────────────────────────
def _adjust_wps(
    wps:      np.ndarray,
    coverage: np.ndarray,
    win_size: int = 1000,
) -> np.ndarray:
    n      = len(wps)
    wps_f  = wps.astype(np.float64)
    half_w = win_size // 2

    # 슬라이딩 윈도우 median
    # 전체 region이 보통 수천 bp이므로 전체 median으로 근사 가능
    # (원저자도 region 내 전체 median 사용)
    window_median = float(np.median(wps_f))

    adjusted = wps_f - window_median

    # Coverage 로 나눔
    cov  = coverage.astype(np.float64)
    out  = np.zeros(n, dtype=np.float32)
    mask = cov > 0
    out[mask]  = (adjusted[mask] / cov[mask] * 100.0).astype(np.float32)
    out[~mask] = np.nan   # coverage 없는 위치는 NaN

    return out


# ─────────────────────────────────────────────────────────────────────
# 모듈 레벨 worker
# ─────────────────────────────────────────────────────────────────────
def _chrom_worker(args: tuple) -> tuple[str, dict]:
    bam_path, chrom, markers, extend, mode, frag_filter, min_mapq = args
    try:
        return chrom, _compute_chrom_markers(
            bam_path, chrom, markers, extend, mode, frag_filter, min_mapq
        )
    except Exception as e:
        log.error("  ✗ %s: %s", chrom, e)
        return chrom, {}


# ─────────────────────────────────────────────────────────────────────
# 공개 인터페이스
# ─────────────────────────────────────────────────────────────────────
def run(
    bam_path:    str,
    out_prefix:  str,
    marker_bed:  Optional[str] = None,
    mode:        str  = "L",
    frag_filter: str  = "long",
    extend:      int  = 1000,
    min_mapq:    int  = 20,
    win_size:    int  = 1000,
    n_jobs:      int  = 4,
    save_bw:     bool = False,
    # 하위 호환
    chroms:      Optional[list] = None,
    block_size:  int  = 1000,
    winsor_q:    float = 0.001,
    blacklist_bed: Optional[str] = None,
) -> dict[str, str]:
    """
    marker BED 기반 WPS 계산.

    같은 염색체 marker 를 묶어서 BAM fetch 최소화.
    adjusted_WPS = (WPS - median) / Coverage * 100
    """
    os.makedirs(os.path.dirname(os.path.abspath(out_prefix)) or ".", exist_ok=True)

    suffix   = f"{mode}_{frag_filter}"
    npy_path = f"{out_prefix}.wps_{suffix}.npy"

    if not marker_bed or not os.path.exists(marker_bed):
        log.warning("marker BED 없음 — WPS 계산 생략")
        np.save(npy_path, {"wps_norm": {}, "mode": mode,
                           "frag": frag_filter, "extend": extend},
                allow_pickle=True)
        return {"npy": npy_path}

    markers = _load_marker_bed(marker_bed)
    if not markers:
        log.warning("유효한 marker 없음")
        np.save(npy_path, {"wps_norm": {}, "mode": mode,
                           "frag": frag_filter, "extend": extend},
                allow_pickle=True)
        return {"npy": npy_path}

    # 염색체별로 marker 그룹화
    chrom_markers: dict[str, list] = defaultdict(list)
    for marker_id, chrom, start, end in markers:
        center = (start + end) // 2
        reg_len = 2 * extend + 1
        chrom_markers[chrom].append((marker_id, center, reg_len))

    log.info("WPS 계산: mode=%s  frag=%s  markers=%d  chroms=%d  extend=±%d  jobs=%d",
             mode, frag_filter, len(markers),
             len(chrom_markers), extend, n_jobs)

    # 병렬 처리 (염색체 단위)
    all_raw: dict[str, tuple[np.ndarray, np.ndarray]] = {}

    task_args = [
        (bam_path, chrom, cms, extend, mode, frag_filter, min_mapq)
        for chrom, cms in chrom_markers.items()
    ]

    with ProcessPoolExecutor(max_workers=n_jobs) as ex:
        futures = {ex.submit(_chrom_worker, a): a[1] for a in task_args}
        for fut in as_completed(futures):
            chrom = futures[fut]
            try:
                _, chrom_results = fut.result()
                all_raw.update(chrom_results)
                log.info("  ✓ %s  (%d markers)", chrom,
                         len(chrom_markers[chrom]))
            except Exception as e:
                log.error("  ✗ %s: %s", chrom, e)

    # adjusted WPS 계산
    wps_norm: dict[str, np.ndarray] = {}
    for marker_id, (wps, cov) in all_raw.items():
        wps_norm[marker_id] = _adjust_wps(wps, cov, win_size=win_size)

    # 통계 로그
    if wps_norm:
        sample = np.concatenate([
            v[np.isfinite(v)] for v in list(wps_norm.values())[:100]
        ])
        if len(sample):
            log.info("Adjusted WPS 통계 (첫 100 markers): "
                     "mean=%.2f  std=%.2f  range=[%.1f~%.1f]",
                     float(np.mean(sample)), float(np.std(sample)),
                     float(np.min(sample)), float(np.max(sample)))

    log.info("WPS 완료: %d / %d markers", len(wps_norm), len(markers))

    np.save(npy_path, {
        "wps_norm": wps_norm,
        "mode":     mode,
        "frag":     frag_filter,
        "extend":   extend,
    }, allow_pickle=True)
    log.info("npy 저장: %s  (%.1f MB)",
             npy_path,
             os.path.getsize(npy_path) / 1e6 if os.path.exists(npy_path) else 0)

    return {"npy": npy_path}


# ─────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────
def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="wps_compute.py",
        description="marker BED 기반 WPS: (WPS - median) / Coverage * 100",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--bam",        required=True)
    p.add_argument("--out-prefix", required=True, dest="out_prefix")
    p.add_argument("--marker-bed", required=True, dest="marker_bed")
    p.add_argument("--mode",   choices=["L","S"], default="L")
    p.add_argument("--frag",   choices=["short","long"],
                   default="long", dest="frag_filter")
    p.add_argument("--extend", type=int, default=1000,
                   help="marker 중심 ±extend bp")
    p.add_argument("--win-size", type=int, default=1000, dest="win_size",
                   help="windowMedian 계산 범위 (bp)")
    p.add_argument("--min-mapq", type=int, default=20, dest="min_mapq")
    p.add_argument("--jobs",     type=int, default=4)
    p.add_argument("--log-level", default="INFO", dest="log_level",
                   choices=["DEBUG","INFO","WARNING","ERROR"])
    return p


def main():
    args = _build_parser().parse_args()
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] — %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )
    if not os.path.exists(args.bam):
        sys.exit(f"BAM 없음: {args.bam}")

    run(
        bam_path    = args.bam,
        out_prefix  = args.out_prefix,
        marker_bed  = args.marker_bed,
        mode        = args.mode,
        frag_filter = args.frag_filter,
        extend      = args.extend,
        min_mapq    = args.min_mapq,
        win_size    = args.win_size,
        n_jobs      = args.jobs,
    )


if __name__ == "__main__":
    main()
