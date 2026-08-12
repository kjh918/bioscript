"""
wps_compute.py
==============
BAM 파일에서 bp 단위 WPS 를 계산하고 bigWig / npy 로 저장합니다.

WPS(i, k) = N_spanning - N_endpoints

출력
----
  {prefix}.wps_{mode}_{frag}.bw       : raw WPS bigWig
  {prefix}.wps_{mode}_{frag}.norm.bw  : 1000bp 블록 median 차감 정규화
  {prefix}.endpoint_{frag}.bw         : endpoint density bigWig
  {prefix}.coverage_{frag}.bw         : fragment coverage bigWig
  {prefix}.wps_{mode}_{frag}.npy      : numpy dict 캐시 (프로파일 플롯용)

사용법
------
python wps_compute.py \\
    --bam sample.bam \\
    --out-prefix ./results/SID001/wps/SID001 \\
    --mode L --frag long --jobs 4

python wps_compute.py \\
    --bam sample.bam \\
    --out-prefix ./results/SID001/wps/SID001 \\
    --mode S --frag short --chroms chr1 chr21
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
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

try:
    import pyBigWig as _pbw
    _HAS_BW = True
except ImportError:
    _HAS_BW = False


# ─────────────────────────────────────────────────────────────────────
# 염색체 단위 WPS 계산
# ─────────────────────────────────────────────────────────────────────
def _compute_chrom(
    chrom:       str,
    chrom_len:   int,
    bam_path:    str,
    mode:        str,
    frag_filter: str,
    min_mapq:    int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Returns
    -------
    wps      : int32 [chrom_len]  WPS = spanning - endpoints
    ep_arr   : int32 [chrom_len]  endpoint density
    cov_arr  : int32 [chrom_len]  coverage depth
    """
    prm    = WPS_PARAMS[mode]
    fmin   = prm["frag_min"]
    fmax   = prm["frag_max"]
    half_k = prm["window"] // 2

    spanning = np.zeros(chrom_len, dtype=np.int32)
    ep_arr   = np.zeros(chrom_len, dtype=np.int32)
    cov_arr  = np.zeros(chrom_len, dtype=np.int32)

    # WPS endpoints 와 endpoint density 를 분리
    wps_ep   = np.zeros(chrom_len, dtype=np.int32)

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
            if frag_filter == "long"  and is_short:
                continue

            # coverage
            lo = max(fs, 0)
            hi = min(fe, chrom_len)
            if lo < hi:
                cov_arr[lo:hi] += 1

            # endpoint density (단순 말단 위치 빈도)
            for ep in (fs, fe - 1):
                if 0 <= ep < chrom_len:
                    ep_arr[ep] += 1

            # WPS (fragment 길이 필터 적용)
            if not (fmin <= flen <= fmax):
                continue

            # spanning
            sp_lo = max(fs + half_k + 1, 0)
            sp_hi = min(fe - half_k,     chrom_len)
            if sp_lo < sp_hi:
                spanning[sp_lo:sp_hi] += 1

            # WPS endpoints (윈도우 내 말단)
            for ep in (fs, fe - 1):
                w_lo = max(ep - half_k,     0)
                w_hi = min(ep + half_k + 1, chrom_len)
                if w_lo < w_hi:
                    wps_ep[w_lo:w_hi] += 1

    wps = (spanning - wps_ep).astype(np.int32)
    return wps, ep_arr, cov_arr


def _normalize(
    arr:        np.ndarray,
    block_size: int   = 1000,
    winsor_q:   float = 0.001,
) -> np.ndarray:
    """
    논문(Snyder et al. 2016) 방식 WPS 정규화.

    1. Winsorization: 상하위 0.1% 값을 clamp
       → centromere/telomere/blacklist 의 극단값 스파이크 제거
    2. 1000bp 비중첩 블록 median 차감 (adjusted WPS)
       → 로컬 coverage 편향 제거, baseline → 0

    winsor_q=0.001 : 상하위 0.1% clamp (전장 유전체 ~300만 bp 기준 3000개)
    """
    out = arr.astype(np.float32).copy()

    # 1. Winsorization — 극단 outlier clamp
    finite = out[np.isfinite(out)]
    if len(finite) > 100:
        lo = float(np.quantile(finite, winsor_q))
        hi = float(np.quantile(finite, 1.0 - winsor_q))
        np.clip(out, lo, hi, out=out)

    # 2. 1000bp 블록 median 차감
    for s in range(0, len(out), block_size):
        e   = min(s + block_size, len(out))
        med = np.nanmedian(out[s:e])
        if np.isfinite(med):
            out[s:e] -= med

    return out


# ─────────────────────────────────────────────────────────────────────
# 모듈 레벨 worker (ProcessPoolExecutor pickle 호환)
# ─────────────────────────────────────────────────────────────────────
def _chrom_worker(
    chrom:       str,
    chrom_len:   int,
    bam_path:    str,
    mode:        str,
    frag_filter: str,
    min_mapq:    int,
):
    """ProcessPoolExecutor 에서 직접 호출되는 모듈 레벨 함수."""
    return chrom, _compute_chrom(chrom, chrom_len, bam_path,
                                 mode, frag_filter, min_mapq)


# ─────────────────────────────────────────────────────────────────────
# bigWig 저장
# ─────────────────────────────────────────────────────────────────────
def _write_bw(
    path:        str,
    chrom_arrs:  dict[str, np.ndarray],
    chrom_sizes: dict[str, int],
) -> None:
    if not _HAS_BW:
        log.warning("pyBigWig 없음 — %s 생략", path)
        return

    # 데이터가 실제로 있는 염색체만 포함
    chrom_order = [
        c for c in chrom_sizes
        if c in chrom_arrs and chrom_arrs[c] is not None and len(chrom_arrs[c]) > 0
    ]

    if not chrom_order:
        log.warning("bigWig 저장 생략 — 유효한 염색체 데이터 없음: %s", path)
        return

    bw = _pbw.open(path, "w")
    bw.addHeader([(c, chrom_sizes[c]) for c in chrom_order])

    for chrom in chrom_order:
        arr = chrom_arrs[chrom].astype(np.float32)

        # NaN 은 0.0 으로 채움 (bigWig 은 NaN 미지원)
        arr = np.where(np.isfinite(arr), arr, 0.0).astype(np.float32)

        # 연속 유효 구간을 묶어서 addEntries (span+step 방식)
        # → 항목 수를 크게 줄여 API 제한 회피
        n = len(arr)
        if n == 0:
            continue

        # 연속 구간 탐지: 값이 0이 아닌 구간만 저장 (0은 skip)
        nonzero = arr != 0.0
        if not nonzero.any():
            continue

        # 구간 경계 계산
        edges   = np.diff(nonzero.astype(np.int8), prepend=0, append=0)
        starts  = np.where(edges == 1)[0]
        ends_   = np.where(edges == -1)[0]

        for s, e in zip(starts.tolist(), ends_.tolist()):
            seg_vals = arr[s:e].tolist()
            if not seg_vals:
                continue
            try:
                bw.addEntries(
                    chrom,
                    int(s),
                    values=seg_vals,
                    span=1,
                    step=1,
                )
            except Exception as bw_err:
                log.warning("bigWig addEntries 실패 [%s:%d-%d]: %s",
                            chrom, s, e, bw_err)

    bw.close()
    log.info("저장: %s", path)


# ─────────────────────────────────────────────────────────────────────
# 공개 인터페이스
# ─────────────────────────────────────────────────────────────────────
def run(
    bam_path:    str,
    out_prefix:  str,
    mode:        str  = "L",
    frag_filter: str  = "long",
    chroms:      Optional[list[str]] = None,
    min_mapq:    int  = 20,
    block_size:  int  = 1000,
    n_jobs:      int  = 4,
) -> dict[str, str]:
    """
    BAM → WPS / endpoint / coverage bigWig + npy 캐시

    Returns
    -------
    out_paths : {name: filepath}
    """
    os.makedirs(os.path.dirname(os.path.abspath(out_prefix)) or ".", exist_ok=True)

    # 염색체 크기
    chrom_sizes: dict[str, int] = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for sq in bam.header["SQ"]:
            chrom_sizes[sq["SN"]] = sq["LN"]

    if chroms is None:
        chroms = [c for c in chrom_sizes if c not in ("chrM","MT","M")]
    else:
        chroms = [c for c in chroms if c in chrom_sizes]

    log.info("WPS 계산 시작: mode=%s  frag=%s  chroms=%d  jobs=%d",
             mode, frag_filter, len(chroms), n_jobs)

    wps_raw:  dict[str, np.ndarray] = {}
    wps_norm: dict[str, np.ndarray] = {}
    ep_dict:  dict[str, np.ndarray] = {}
    cov_dict: dict[str, np.ndarray] = {}

    with ProcessPoolExecutor(max_workers=n_jobs) as ex:
        futures = {
            ex.submit(_chrom_worker, c, chrom_sizes[c],
                      bam_path, mode, frag_filter, min_mapq): c
            for c in chroms
        }
        for fut in as_completed(futures):
            chrom = futures[fut]
            try:
                _, (wps, ep, cov) = fut.result()
                raw_f32         = wps.astype(np.float32)
                norm_arr        = _normalize(wps, block_size)
                wps_raw[chrom]  = raw_f32
                wps_norm[chrom] = norm_arr
                ep_dict[chrom]  = ep.astype(np.float32)
                cov_dict[chrom] = cov.astype(np.float32)
                log.info(
                    "  ✓ %s  raw=[%.1f~%.1f]  norm=[%.3f~%.3f]",
                    chrom,
                    float(wps.min()), float(wps.max()),
                    float(norm_arr.min()), float(norm_arr.max()),
                )
            except Exception as exc:
                log.error("  ✗ %s: %s", chrom, exc)

    if not wps_raw:
        log.error("WPS 계산 결과 없음 — 모든 염색체 worker 실패. BAM/파라미터 확인 필요.")
        return {}

    log.info("WPS 계산 완료: %d / %d 염색체", len(wps_raw), len(chroms))

    suffix = f"{mode}_{frag_filter}"
    out_paths = {
        "wps_raw":  f"{out_prefix}.wps_{suffix}.bw",
        "wps_norm": f"{out_prefix}.wps_{suffix}.norm.bw",
        "endpoint": f"{out_prefix}.endpoint_{frag_filter}.bw",
        "coverage": f"{out_prefix}.coverage_{frag_filter}.bw",
        "npy":      f"{out_prefix}.wps_{suffix}.npy",
    }

    _write_bw(out_paths["wps_raw"],  wps_raw,  chrom_sizes)
    _write_bw(out_paths["wps_norm"], wps_norm, chrom_sizes)
    _write_bw(out_paths["endpoint"], ep_dict,  chrom_sizes)
    _write_bw(out_paths["coverage"], cov_dict, chrom_sizes)

    # npy 캐시: dict {chrom: array} 저장 (allow_pickle=True 필요)
    np.save(out_paths["npy"], {
        "wps_raw":  wps_raw,
        "wps_norm": wps_norm,
        "endpoint": ep_dict,
        "coverage": cov_dict,
        "mode":     mode,
        "frag":     frag_filter,
    }, allow_pickle=True)
    log.info("npy 캐시 저장: %s", out_paths["npy"])

    return out_paths


# ─────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────
def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="wps_compute.py",
        description="BAM → bp 단위 WPS / endpoint / coverage bigWig",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--bam",        required=True)
    p.add_argument("--out-prefix", required=True, dest="out_prefix",
                   help="출력 prefix (디렉토리/파일명 앞부분)")
    p.add_argument("--mode",   choices=["L","S"], default="L",
                   help="L=뉴클레오솜(k=120,frag 120-180bp) "
                        "S=TF footprint(k=16,frag 35-80bp)")
    p.add_argument("--frag",   choices=["short","long","all"],
                   default="long", dest="frag_filter",
                   help="short≤150bp | long≥151bp | all")
    p.add_argument("--chroms", nargs="*", default=None,
                   help="처리 염색체 (미지정=전체, chrM 제외)")
    p.add_argument("--min-mapq",   type=int, default=20,  dest="min_mapq")
    p.add_argument("--block-size", type=int, default=1000, dest="block_size",
                   help="정규화 블록 크기 (bp)")
    p.add_argument("--jobs",       type=int, default=4)
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
    if not os.path.exists(args.bam):
        sys.exit(f"BAM 없음: {args.bam}")

    run(
        bam_path    = args.bam,
        out_prefix  = args.out_prefix,
        mode        = args.mode,
        frag_filter = args.frag_filter,
        chroms      = args.chroms,
        min_mapq    = args.min_mapq,
        block_size  = args.block_size,
        n_jobs      = args.jobs,
    )


if __name__ == "__main__":
    main()