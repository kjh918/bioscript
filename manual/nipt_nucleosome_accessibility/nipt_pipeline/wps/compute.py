"""
wps/compute.py
==============
BAM → chromosome 단위 WPS NPY (multiprocess)

출력 파일 (chromosome당 4종)
----------------------------
{out_dir}/{mode}/{prefix}.{chrom}.raw.npy      int32   spanning - endpoints
{out_dir}/{mode}/{prefix}.{chrom}.cov.npy      int32   전체 fragment coverage (norm 분모)
{out_dir}/{mode}/{prefix}.{chrom}.frag_cov.npy int32   mode별 fragment coverage (시각화용)
{out_dir}/{mode}/{prefix}.{chrom}.norm.npy     float32 adjusted WPS = (WPS - blockMedian) / cov * 100

L-WPS: frag 120-180bp, k=120
S-WPS: frag 35-80bp,   k=16
"""
from __future__ import annotations

import logging
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Optional

import numpy as np
import pysam

from .utils import (
    HG38_CHROM_SIZES, STANDARD_CHROMS, WPS_PARAMS, MAX_FRAG,
    track_path, write_manifest,
)

log = logging.getLogger(__name__)

# pysam FLAG
_FLAG_DUP         = 0x400
_FLAG_SECONDARY   = 0x100
_FLAG_SUPPLEMENTARY = 0x800


# ── chromosome 단위 WPS 계산 ─────────────────────────────────────────
def _calc_chrom(
    chrom:     str,
    chrom_len: int,
    bam_path:  str,
    mode:      str,
    min_mapq:  int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Returns
    -------
    wps      : int32[chrom_len]  spanning - endpoints (mode별 frag 범위)
    cov_all  : int32[chrom_len]  전체 fragment coverage (norm 분모)
    cov_frag : int32[chrom_len]  mode별 fragment coverage (시각화용)
    """
    prm    = WPS_PARAMS[mode]
    fmin   = prm["frag_min"]
    fmax   = prm["frag_max"]
    half_k = prm["window"] // 2

    spanning  = np.zeros(chrom_len, dtype=np.int32)
    endpoints = np.zeros(chrom_len, dtype=np.int32)
    cov_all   = np.zeros(chrom_len, dtype=np.int32)
    cov_frag  = np.zeros(chrom_len, dtype=np.int32)

    seen: set[str] = set()

    with pysam.AlignmentFile(bam_path, "rb", threads=1) as bam:
        for read in bam.fetch(chrom):
            # 기본 필터
            if (read.is_unmapped
                    or read.is_duplicate
                    or read.is_secondary
                    or read.is_supplementary):
                continue
            if read.mapping_quality < min_mapq:
                continue
            # paired-end: read1만 처리 (pair당 1회 카운트)
            if read.is_paired and read.template_length <= 0:
                continue

            # qname 중복 방지
            qname = read.query_name
            if qname in seen:
                continue
            seen.add(qname)

            # fragment 좌표 결정
            if read.is_paired and read.template_length != 0:
                tlen = abs(int(read.template_length))
                if 0 < tlen <= MAX_FRAG:
                    fs = int(read.reference_start)
                    fe = fs + tlen
                else:
                    fs = int(read.reference_start)
                    fe = int(read.reference_end) if read.reference_end else fs
            else:
                fs = int(read.reference_start)
                fe = int(read.reference_end) if read.reference_end else fs

            if fs < 0 or fe <= fs:
                continue
            flen = fe - fs
            if flen <= 0 or flen > MAX_FRAG:
                continue

            clo = max(fs, 0)
            chi = min(fe, chrom_len)
            if clo >= chi:
                continue

            # 전체 coverage 누적 (norm 분모)
            cov_all[clo:chi] += 1

            # mode별 frag 범위에 해당하는 경우에만 WPS 계산
            if fmin <= flen <= fmax:
                cov_frag[clo:chi] += 1

                # spanning: fragment 내부 window 제외 영역
                sp_lo = max(fs + half_k + 1, 0)
                sp_hi = min(fe - half_k,     chrom_len)
                if sp_lo < sp_hi:
                    spanning[sp_lo:sp_hi] += 1

                # endpoints: 양 끝 ± half_k window
                for ep in (fs, fe - 1):
                    w_lo = max(ep - half_k,     0)
                    w_hi = min(ep + half_k + 1, chrom_len)
                    if w_lo < w_hi:
                        endpoints[w_lo:w_hi] += 1

    wps = (spanning - endpoints).astype(np.int32)
    return wps, cov_all, cov_frag


# ── 정규화: WPS - blockMedian + Savitzky-Golay smoothing ─────────────
# Snyder et al. 2016:
#   1. 1kb window running median subtraction (locally zero-adjusted)
#   2. Savitzky-Golay filter (window=21, poly=2) smoothing
# coverage 나누기는 적용하지 않음
# aggregate 단계에서 trimmed mean subtraction으로 샘플간 정규화
def _normalize(
    wps:       np.ndarray,
    coverage:  np.ndarray,   # coverage 보존 (cov.npy용, norm 계산엔 미사용)
    win_size:  int = 1000,
    sg_window: int = 21,
    sg_poly:   int = 2,
) -> np.ndarray:
    from scipy.signal import savgol_filter

    n    = int(wps.size)
    norm = np.zeros(n, dtype=np.float32)  # NaN 대신 0으로 초기화

    # step 1: 1kb block median subtraction
    # coverage > 0 위치만 median 계산에 사용하되,
    # WPS 값은 전체 위치에 대해 기록 (coverage=0 → WPS=0으로 처리)
    # → Snyder 2016: locally zero-adjusted, 위치 제한 없음
    for s in range(0, n, win_size):
        e         = min(s + win_size, n)
        wps_block = wps[s:e].astype(np.float32)
        cov_block = coverage[s:e]
        cov_mask  = cov_block > 0

        if not cov_mask.any():
            # 이 block은 coverage가 전혀 없음 → 0으로 채움
            norm[s:e] = 0.0
            continue

        # median은 coverage 있는 위치만으로 계산
        med       = float(np.median(wps_block[cov_mask]))
        # subtraction은 전체 위치에 적용
        norm[s:e] = wps_block - med

    # step 2: Savitzky-Golay smoothing
    # norm이 0으로 초기화되어 있으므로 전체 배열에 한 번에 적용 가능
    # coverage 없는 구간(값=0)도 smoothing에 포함 → 경계 artifact 최소화
    if n >= sg_window:
        norm = savgol_filter(
            norm.astype(np.float64),
            window_length = sg_window,
            polyorder     = sg_poly,
        ).astype(np.float32)

    return norm


# ── ProcessPoolExecutor worker (module-level, pickle 호환) ─────────────
def _chrom_worker(args: tuple) -> tuple[str, dict]:
    """
    args = (chrom, chrom_len, bam_path, mode, min_mapq, win_size,
            out_dir, prefix, resume)
    Returns (chrom, stats_dict)
    """
    chrom, chrom_len, bam_path, mode, min_mapq, win_size, out_dir, prefix, resume = args

    raw_path      = track_path(out_dir, mode, chrom, "raw",      prefix)
    cov_path      = track_path(out_dir, mode, chrom, "cov",      prefix)
    frag_cov_path = track_path(out_dir, mode, chrom, "frag_cov", prefix)
    norm_path     = track_path(out_dir, mode, chrom, "norm",     prefix)

    all_paths = (raw_path, cov_path, frag_cov_path, norm_path)

    # resume: 4종 모두 존재하면 건너뜀
    if resume and all(os.path.isfile(p) for p in all_paths):
        raw  = np.load(raw_path,  mmap_mode="r", allow_pickle=False)
        norm = np.load(norm_path, mmap_mode="r", allow_pickle=False)
        finite = norm[np.isfinite(norm)]
        return chrom, {
            "raw":      raw_path,
            "cov":      cov_path,
            "frag_cov": frag_cov_path,
            "norm":     norm_path,
            "wps_min":  int(raw.min())       if raw.size    else 0,
            "wps_max":  int(raw.max())       if raw.size    else 0,
            "norm_min": float(finite.min())  if finite.size else float("nan"),
            "norm_max": float(finite.max())  if finite.size else float("nan"),
            "resumed":  True,
        }

    wps, cov_all, cov_frag = _calc_chrom(
        chrom, chrom_len, bam_path, mode, min_mapq
    )
    norm = _normalize(wps, cov_all, win_size=win_size)

    # 원자적 저장 (tmp → rename)
    # worker 프로세스에서 디렉토리가 없을 수 있으므로 makedirs 보장
    os.makedirs(os.path.dirname(raw_path), exist_ok=True)
    for arr, path in [
        (wps,      raw_path),
        (cov_all,  cov_path),
        (cov_frag, frag_cov_path),
        (norm,     norm_path),
    ]:
        # np.save는 확장자가 없으면 .npy를 자동 추가하므로
        # tmp 경로도 .npy로 끝나야 os.replace 대상과 일치함
        # path = "...chrX.raw.npy" → tmp = "...chrX.raw.tmp.npy"
        tmp = path.replace(".npy", ".tmp.npy")
        np.save(tmp, arr, allow_pickle=False)
        os.replace(tmp, path)

    finite = norm[np.isfinite(norm)]
    return chrom, {
        "raw":      raw_path,
        "cov":      cov_path,
        "frag_cov": frag_cov_path,
        "norm":     norm_path,
        "wps_min":  int(wps.min())       if wps.size    else 0,
        "wps_max":  int(wps.max())       if wps.size    else 0,
        "norm_min": float(finite.min())  if finite.size else float("nan"),
        "norm_max": float(finite.max())  if finite.size else float("nan"),
        "finite_fraction": float(finite.size / norm.size) if norm.size else 0.0,
    }


# ── 공개 인터페이스 ────────────────────────────────────────────────────
def run(
    bam_path: str,
    out_dir:  str,
    mode:     str,
    chroms:   Optional[list[str]] = None,
    min_mapq: int  = 20,
    win_size: int  = 1000,
    n_jobs:   int  = 4,
    resume:   bool = False,
    prefix:   str  = "",
) -> dict:
    """
    BAM → {out_dir}/{mode}/ 에 chromosome별 4종 NPY 저장.

    Parameters
    ----------
    bam_path : sorted + indexed BAM
    out_dir  : 출력 루트 (하위에 mode 디렉토리 생성)
    mode     : "L" (long, 120-180bp) | "S" (short, 35-80bp)
    chroms   : 처리할 염색체 목록 (None이면 BAM 헤더 기반 standard)
    min_mapq : 최소 MAPQ
    win_size : blockMedian 계산 window 크기 (bp)
    n_jobs   : 병렬 프로세스 수
    resume   : 기존 NPY 있으면 건너뜀
    prefix   : 출력 파일명 prefix (예: sample_id)

    Returns
    -------
    {
      "mode": "L",
      "dir":  "{out_dir}/L",
      "manifest": "...",
      "chromosomes": {chrom: {raw, cov, frag_cov, norm, ...}, ...}
    }
    """
    if mode not in WPS_PARAMS:
        raise ValueError(f"지원하지 않는 mode: {mode!r}. L 또는 S 중 선택")
    if not os.path.isfile(bam_path):
        raise FileNotFoundError(f"BAM 없음: {bam_path}")

    track_dir = os.path.join(out_dir, mode)
    os.makedirs(track_dir, exist_ok=True)

    # BAM 헤더에서 chrom 크기 로드
    chrom_sizes: dict[str, int] = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for sq in bam.header.get("SQ", []):
            chrom_sizes[str(sq["SN"])] = int(sq["LN"])

    # 처리 대상 염색체 결정
    canonical = [c for c in STANDARD_CHROMS if c in chrom_sizes]
    target    = chroms if chroms else canonical
    target    = [c for c in target if c in chrom_sizes]

    prm = WPS_PARAMS[mode]
    log.info(
        "WPS compute mode=%s frag=%d-%dbp k=%d chroms=%d jobs=%d",
        mode, prm["frag_min"], prm["frag_max"], prm["window"], len(target), n_jobs,
    )

    task_args = [
        (chrom, chrom_sizes[chrom], bam_path, mode,
         min_mapq, win_size, out_dir, prefix, resume)
        for chrom in target
    ]

    result: dict[str, dict] = {}

    with ProcessPoolExecutor(max_workers=n_jobs) as ex:
        futures = {ex.submit(_chrom_worker, a): a[0] for a in task_args}
        for fut in as_completed(futures):
            chrom = futures[fut]
            try:
                chrom_out, stats = fut.result()
                result[chrom_out] = stats
                nm_min = f"{stats['norm_min']:.2f}" if np.isfinite(stats["norm_min"]) else "NA"
                nm_max = f"{stats['norm_max']:.2f}" if np.isfinite(stats["norm_max"]) else "NA"
                log.info(
                    "  ✓ %s wps=[%d~%d] norm=[%s~%s]%s",
                    chrom_out, stats["wps_min"], stats["wps_max"], nm_min, nm_max,
                    " (resumed)" if stats.get("resumed") else "",
                )
            except Exception as exc:
                log.exception("  ✗ %s: %s", chrom, exc)

    # 염색체 순서 유지
    ordered = {c: result[c] for c in target if c in result}

    manifest_path = write_manifest(track_dir, mode, ordered)
    log.info("compute 완료 mode=%s chroms=%d manifest=%s", mode, len(ordered), manifest_path)

    return {
        "mode":        mode,
        "dir":         track_dir,
        "manifest":    manifest_path,
        "chromosomes": ordered,
    }