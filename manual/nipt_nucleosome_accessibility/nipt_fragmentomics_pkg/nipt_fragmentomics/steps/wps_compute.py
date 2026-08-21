"""
wps_compute.py
==============
BAM → 염색체별 WPS NPY (raw / cov / norm 3종)

출력 구조
---------
{out_dir}/{mode}/
    chr1.raw.npy   : int32  spanning - endpoints (raw WPS)
    chr1.cov.npy   : int32  coverage (WPS fragment 범위 기준)
    chr1.norm.npy  : float32 adjusted WPS = (WPS - blockMedian) / cov * 100
    manifest.json

L-WPS: frag 120-180bp, k=120
S-WPS: frag 35-80bp,  k=16
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Optional

import numpy as np
import pysam

log = logging.getLogger(__name__)

BUILD_TAG = "LS-chrom-npy-3track-20260820-v1"

WPS_PARAMS = {
    "L": {"frag_min": 120, "frag_max": 180, "window": 120},
    "S": {"frag_min":  35, "frag_max":  80, "window":  16},
}
MAX_FRAG = 1000
SHORT_MAX = 150


# ─────────────────────────────────────────────────────────────────────
# 염색체 단위 WPS 계산
# ─────────────────────────────────────────────────────────────────────
def _calc_chrom(
    chrom:    str,
    chrom_len: int,
    bam_path: str,
    mode:     str,
    min_mapq: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Returns
    -------
    wps      : int32 [chrom_len]  spanning - endpoints (mode별 frag 범위)
    cov_all  : int32 [chrom_len]  전체 fragment coverage (원저자 normalization 분모)
    cov_frag : int32 [chrom_len]  mode별 frag 범위 coverage (L: 120-180bp / S: 35-80bp)
    """
    prm    = WPS_PARAMS[mode]
    fmin   = prm["frag_min"]
    fmax   = prm["frag_max"]
    half_k = prm["window"] // 2

    spanning  = np.zeros(chrom_len, dtype=np.int32)
    endpoints = np.zeros(chrom_len, dtype=np.int32)
    cov_all   = np.zeros(chrom_len, dtype=np.int32)   # 전체 coverage
    cov_frag  = np.zeros(chrom_len, dtype=np.int32)   # mode별 fragment coverage

    seen: set[str] = set()

    with pysam.AlignmentFile(bam_path, "rb", threads=1) as bam:
        for read in bam.fetch(chrom):
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

            # fragment 좌표
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

            clo = max(fs, 0);  chi = min(fe, chrom_len)

            # 전체 coverage (원저자 normalization 분모)
            if clo < chi:
                cov_all[clo:chi] += 1

            # mode별 fragment coverage
            if fmin <= flen <= fmax:
                if clo < chi:
                    cov_frag[clo:chi] += 1

                # spanning
                sp_lo = max(fs + half_k + 1, 0)
                sp_hi = min(fe - half_k,     chrom_len)
                if sp_lo < sp_hi:
                    spanning[sp_lo:sp_hi] += 1

                # endpoints
                for ep in (fs, fe - 1):
                    w_lo = max(ep - half_k,     0)
                    w_hi = min(ep + half_k + 1, chrom_len)
                    if w_lo < w_hi:
                        endpoints[w_lo:w_hi] += 1

    wps = (spanning - endpoints).astype(np.int32)
    return wps, cov_all, cov_frag


# ─────────────────────────────────────────────────────────────────────
# 정규화: (WPS - blockMedian) / coverage * 100
# ─────────────────────────────────────────────────────────────────────
def _normalize(
    wps:      np.ndarray,
    coverage: np.ndarray,
    win_size: int = 1000,
) -> np.ndarray:
    n    = int(wps.size)
    norm = np.full(n, np.nan, dtype=np.float32)

    for s in range(0, n, win_size):
        e         = min(s + win_size, n)
        wps_block = wps[s:e].astype(np.float32)
        cov_block = coverage[s:e]
        med       = float(np.median(wps_block))
        mask      = cov_block > 0
        if not mask.any():
            continue
        adjusted       = wps_block[mask] - med
        norm_block     = norm[s:e]
        norm_block[mask] = (
            adjusted / cov_block[mask].astype(np.float32) * 100.0
        ).astype(np.float32)

    return norm


# ─────────────────────────────────────────────────────────────────────
# 모듈 레벨 worker (ProcessPoolExecutor pickle 호환)
# ─────────────────────────────────────────────────────────────────────
def _chrom_worker(args: tuple) -> tuple[str, dict]:
    """
    (chrom, chrom_len, bam_path, mode, min_mapq, win_size, track_dir, resume)
    → (chrom, stats_dict)

    raw / cov / norm 3종을 {track_dir}/{chrom}.{metric}.npy 로 저장합니다.
    """
    chrom, chrom_len, bam_path, mode, min_mapq, win_size, track_dir, resume, prefix = args

    fname_base   = f"{prefix}.{chrom}" if prefix else chrom
    raw_path     = os.path.join(track_dir, f"{fname_base}.raw.npy")
    cov_path     = os.path.join(track_dir, f"{fname_base}.cov.npy")
    norm_path    = os.path.join(track_dir, f"{fname_base}.norm.npy")
    frag_cov_path = os.path.join(track_dir, f"{fname_base}.frag_cov.npy")

    # resume: 4종 모두 있으면 건너뜀
    if resume and all(os.path.isfile(p) for p in (raw_path, cov_path, norm_path, frag_cov_path)):
        raw  = np.load(raw_path,  mmap_mode="r", allow_pickle=False)
        norm = np.load(norm_path, mmap_mode="r", allow_pickle=False)
        finite = norm[np.isfinite(norm)]
        return chrom, {
            "raw": raw_path, "cov": cov_path, "norm": norm_path,
            "frag_cov": frag_cov_path,
            "wps_min": int(raw.min()) if raw.size else 0,
            "wps_max": int(raw.max()) if raw.size else 0,
            "norm_min": float(finite.min()) if finite.size else float("nan"),
            "norm_max": float(finite.max()) if finite.size else float("nan"),
            "resumed": True,
        }

    wps, cov_all, cov_frag = _calc_chrom(chrom, chrom_len, bam_path, mode, min_mapq)
    # normalization 분모: 원저자 방식 (전체 coverage)
    norm = _normalize(wps, cov_all, win_size=win_size)

    # 원자적 저장
    for arr, path in [
        (wps,      raw_path),
        (cov_all,  cov_path),      # 전체 coverage (normalization 분모)
        (norm,     norm_path),
        (cov_frag, frag_cov_path), # mode별 fragment coverage (시각화용)
    ]:
        tmp = path + ".tmp.npy"
        np.save(tmp, arr, allow_pickle=False)
        os.replace(tmp, path)

    finite = norm[np.isfinite(norm)]
    stats = {
        "raw":       raw_path,
        "cov":       cov_path,
        "norm":      norm_path,
        "frag_cov":  frag_cov_path,
        "wps_min":   int(wps.min())       if wps.size    else 0,
        "wps_max":   int(wps.max())       if wps.size    else 0,
        "norm_min":  float(finite.min())  if finite.size else float("nan"),
        "norm_max":  float(finite.max())  if finite.size else float("nan"),
        "finite_fraction": float(finite.size / norm.size) if norm.size else 0.0,
    }
    return chrom, stats


# ─────────────────────────────────────────────────────────────────────
# 공개 인터페이스
# ─────────────────────────────────────────────────────────────────────
def run(
    bam_path: str,
    out_dir:  str,
    mode:     str            = "L",
    chroms:   Optional[list[str]] = None,
    min_mapq: int            = 20,
    win_size: int            = 1000,
    n_jobs:   int            = 4,
    resume:   bool           = False,
    prefix:   str            = "",     # 파일명 prefix (예: sample_id)
    # 하위 호환 파라미터 (무시됨)
    **kwargs,
) -> dict:
    """
    BAM → {out_dir}/{mode}/chr*.raw/cov/norm.npy + manifest.json

    Returns
    -------
    {
      "mode": "L",
      "dir":  "{out_dir}/L",
      "manifest": "{out_dir}/L/manifest.json",
      "chromosomes": {
          "chr1": {"raw": "...", "cov": "...", "norm": "..."},
          ...
      }
    }
    """
    if mode not in WPS_PARAMS:
        raise ValueError(f"지원하지 않는 WPS mode: {mode}")
    if not os.path.exists(bam_path):
        raise FileNotFoundError(f"BAM 없음: {bam_path}")

    prm       = WPS_PARAMS[mode]
    track_dir = os.path.join(out_dir, mode)
    os.makedirs(track_dir, exist_ok=True)

    # 염색체 크기 로드
    chrom_sizes: dict[str, int] = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for sq in bam.header.get("SQ", []):
            chrom_sizes[str(sq["SN"])] = int(sq["LN"])

    canonical = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    if chroms is None:
        chroms = [c for c in canonical if c in chrom_sizes]
    else:
        chroms = [c for c in chroms if c in chrom_sizes]

    log.info(
        "WPS build=%s mode=%s frag=%d-%dbp k=%d chroms=%d jobs=%d dir=%s",
        BUILD_TAG, mode, prm["frag_min"], prm["frag_max"],
        prm["window"], len(chroms), n_jobs, track_dir,
    )

    task_args = [
        (chrom, chrom_sizes[chrom], bam_path, mode,
         min_mapq, win_size, track_dir, resume, prefix)
        for chrom in chroms
    ]

    result_paths: dict[str, dict] = {}

    with ProcessPoolExecutor(max_workers=n_jobs) as ex:
        futures = {ex.submit(_chrom_worker, a): a[0] for a in task_args}
        for fut in as_completed(futures):
            chrom = futures[fut]
            try:
                chrom, stats = fut.result()
                result_paths[chrom] = stats

                norm_min = f"{stats['norm_min']:.2f}" if np.isfinite(stats['norm_min']) else "NA"
                norm_max = f"{stats['norm_max']:.2f}" if np.isfinite(stats['norm_max']) else "NA"
                log.info(
                    "  ✓ %s wps=[%d~%d] norm=[%s~%s]%s",
                    chrom, stats["wps_min"], stats["wps_max"],
                    norm_min, norm_max,
                    " (resumed)" if stats.get("resumed") else "",
                )
            except Exception as exc:
                log.exception("  ✗ %s: %s", chrom, exc)

    # manifest 저장
    ordered = {c: result_paths[c] for c in chroms if c in result_paths}
    manifest_path = os.path.join(track_dir, "manifest.json")
    with open(manifest_path, "w") as f:
        json.dump({
            "build": BUILD_TAG, "mode": mode,
            "frag_min": prm["frag_min"], "frag_max": prm["frag_max"],
            "window": prm["window"], "win_size": win_size,
            "chromosomes": ordered,
        }, f, indent=2)

    log.info("WPS 완료 mode=%s chroms=%d manifest=%s",
             mode, len(ordered), manifest_path)

    return {
        "mode":        mode,
        "dir":         track_dir,
        "manifest":    manifest_path,
        "chromosomes": ordered,
    }


# ─────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────
def _build_parser():
    p = argparse.ArgumentParser(
        prog="wps_compute.py",
        description="BAM → L/S WPS NPY (raw/cov/norm 3종)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--bam",      required=True)
    p.add_argument("--out-dir",  required=True, dest="out_dir")
    p.add_argument("--mode",     choices=["L","S"], required=True)
    p.add_argument("--chroms",   nargs="*", default=None)
    p.add_argument("--min-mapq", type=int, default=20, dest="min_mapq")
    p.add_argument("--win-size", type=int, default=1000, dest="win_size")
    p.add_argument("--jobs",     type=int, default=4)
    p.add_argument("--resume",   action="store_true")
    p.add_argument("--log-level", default="INFO", dest="log_level",
                   choices=["DEBUG","INFO","WARNING","ERROR"])
    return p


def main():
    args = _build_parser().parse_args()
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(name)s — %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )
    if not os.path.exists(args.bam):
        sys.exit(f"BAM 없음: {args.bam}")
    run(
        bam_path=args.bam, out_dir=args.out_dir,
        mode=args.mode, chroms=args.chroms,
        min_mapq=args.min_mapq, win_size=args.win_size,
        n_jobs=args.jobs, resume=args.resume,
    )


if __name__ == "__main__":
    main()