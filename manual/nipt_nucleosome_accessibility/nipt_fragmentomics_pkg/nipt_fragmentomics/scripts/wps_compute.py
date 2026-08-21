"""
wps_compute.py
==============
BAM → bp 단위 WPS bigWig (원저자 방식)

원저자(Snyder et al. 2016) 방식:
  WPS(i) = N_spanning(i) - N_endpoints(i)

  N_spanning(i) : window [i-k/2, i+k/2] 를 완전히 가로지르는 fragment 수
                  frag_start < i-k/2  AND  frag_end > i+k/2
  N_endpoints(i): window 내에 5' 또는 3' 말단이 있는 fragment 수

계산 단위: position별 (bp 단위 배열)
병렬: 염색체별 독립 ProcessPoolExecutor

정규화 (normalizeWPSwigs.py):
  adjusted_WPS = (WPS - windowMedian) / Coverage * 100
  windowMedian: 1kb 슬라이딩 윈도우 median

출력:
  {prefix}.{mode}_{frag}.WPS.bw      : raw WPS bigWig
  {prefix}.{mode}_{frag}.WPS.norm.bw : adjusted WPS bigWig
  {prefix}.{mode}_{frag}.Coverage.bw : coverage bigWig

사용법:
  python wps_compute.py \\
      --bam sample.bam \\
      --out-prefix ./results/SID001/wps/SID001 \\
      --mode L --frag long \\
      --jobs 8

  # 특정 염색체만
  python wps_compute.py \\
      --bam sample.bam \\
      --out-prefix ./results/SID001/wps/SID001 \\
      --mode L --frag long \\
      --chroms chr1 chr2 chr21 \\
      --jobs 4
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

# WPS 파라미터 (Snyder et al. 2016)
WPS_PARAMS = {
    "L": {"frag_min": 120, "frag_max": 180, "window": 120},
    "S": {"frag_min":  35, "frag_max":  80, "window":  16},
}
SHORT_MAX   = 150
MAX_FRAG    = 1000   # 비정상 template_length 필터

try:
    import pyBigWig as _pbw
    _HAS_BW = True
except ImportError:
    _HAS_BW = False
    log.warning("pyBigWig 미설치 — bigWig 저장 불가. pip install pyBigWig")


# ─────────────────────────────────────────────────────────────────────
# 염색체 단위 WPS 계산 (모듈 레벨 — pickle 가능)
# ─────────────────────────────────────────────────────────────────────
def _calc_chrom(
    chrom:       str,
    chrom_len:   int,
    bam_path:    str,
    mode:        str,
    frag_filter: str,   # "short" | "long" | "all"
    min_mapq:    int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    염색체 전체에 대해 bp 단위 WPS 배열과 Coverage 배열을 반환합니다.

    Returns
    -------
    wps      : int32 [chrom_len]   WPS(i) = spanning(i) - endpoints(i)
    coverage : int32 [chrom_len]   fragment coverage depth
    """
    prm    = WPS_PARAMS[mode]
    fmin   = prm["frag_min"]
    fmax   = prm["frag_max"]
    half_k = prm["window"] // 2

    # position별 배열
    spanning  = np.zeros(chrom_len, dtype=np.int32)
    endpoints = np.zeros(chrom_len, dtype=np.int32)
    coverage  = np.zeros(chrom_len, dtype=np.int32)

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

            # fragment 좌표 결정
            if read.is_paired and read.template_length != 0:
                tlen = abs(read.template_length)
                if tlen > MAX_FRAG or tlen < (read.query_length or 0):
                    fs = read.reference_start
                    fe = read.reference_end
                else:
                    fs = read.reference_start
                    fe = fs + tlen
            else:
                fs = read.reference_start
                fe = read.reference_end

            flen = fe - fs
            if flen <= 0 or flen > MAX_FRAG:
                continue

            # fragment 필터
            is_short = (flen <= SHORT_MAX)
            if frag_filter == "short" and not is_short:
                continue
            if frag_filter == "long"  and is_short:
                continue

            # ── coverage: [fs, fe) 범위 ──────────────────────────
            lo = max(fs, 0)
            hi = min(fe, chrom_len)
            if lo < hi:
                coverage[lo:hi] += 1

            # ── WPS: fragment 길이 필터 적용 ─────────────────────
            if not (fmin <= flen <= fmax):
                continue

            # spanning: position i 에서 [i-k/2, i+k/2] 를 완전히 덮음
            # → frag_start < i - half_k  AND  frag_end > i + half_k
            # → i 범위: (fs+half_k, fe-half_k)  exclusive
            sp_lo = max(fs + half_k + 1, 0)
            sp_hi = min(fe - half_k,     chrom_len)
            if sp_lo < sp_hi:
                spanning[sp_lo:sp_hi] += 1

            # endpoints: 5'(fs) 와 3'(fe-1) 말단
            # 말단 ep 에 대해 window [ep-half_k, ep+half_k] 내 position 에 기여
            for ep in (fs, fe - 1):
                w_lo = max(ep - half_k,     0)
                w_hi = min(ep + half_k + 1, chrom_len)
                if w_lo < w_hi:
                    endpoints[w_lo:w_hi] += 1

    wps = (spanning - endpoints).astype(np.int32)
    return wps, coverage


# ─────────────────────────────────────────────────────────────────────
# 정규화: (WPS - windowMedian) / Coverage * 100
# ─────────────────────────────────────────────────────────────────────
def _normalize(
    wps:      np.ndarray,
    coverage: np.ndarray,
    win_size: int = 1000,
) -> np.ndarray:
    """
    원저자 normalizeWPSwigs.py 방식.
    windowMedian: 각 position i 주변 win_size bp median.
    Coverage = 0 → 0.0.
    """
    n      = len(wps)
    wps_f  = wps.astype(np.float64)
    half_w = win_size // 2

    # 비중첩 블록 median으로 근사 (속도 최적화)
    wmed = np.zeros(n, dtype=np.float64)
    for s in range(0, n, win_size):
        e   = min(s + win_size, n)
        med = float(np.median(wps_f[s:e]))
        wmed[s:e] = med

    adjusted = wps_f - wmed
    cov      = coverage.astype(np.float64)

    norm = np.zeros(n, dtype=np.float32)
    mask = cov > 0
    norm[mask] = (adjusted[mask] / cov[mask] * 100.0).astype(np.float32)
    return norm


# ─────────────────────────────────────────────────────────────────────
# bigWig 저장
# ─────────────────────────────────────────────────────────────────────
def _write_bw(
    path:        str,
    chrom_arrs:  dict[str, np.ndarray],
    chrom_sizes: dict[str, int],
    chrom_order: list[str],
) -> None:
    if not _HAS_BW:
        log.warning("pyBigWig 없음 — %s 생략", path)
        return

    present = [c for c in chrom_order
               if c in chrom_arrs and len(chrom_arrs[c]) > 0]
    if not present:
        log.warning("bigWig 저장 생략: 데이터 없음 (%s)", path)
        return

    bw = _pbw.open(path, "w")
    bw.addHeader([(c, chrom_sizes[c]) for c in present])

    for chrom in present:
        arr = np.asarray(chrom_arrs[chrom], dtype=np.float32)
        arr = np.where(np.isfinite(arr), arr, 0.0)

        # 0이 아닌 연속 구간만 저장 (파일 크기 절감)
        nonzero = arr != 0.0
        if not nonzero.any():
            continue

        edges  = np.diff(nonzero.astype(np.int8), prepend=0, append=0)
        starts = np.where(edges ==  1)[0]
        ends   = np.where(edges == -1)[0]

        for s, e in zip(starts.tolist(), ends.tolist()):
            seg = arr[s:e].tolist()
            if not seg:
                continue
            try:
                bw.addEntries(chrom, int(s), values=seg, span=1, step=1)
            except Exception as ex:
                log.debug("bigWig addEntries 실패 [%s:%d-%d]: %s", chrom, s, e, ex)

    bw.close()
    log.info("저장: %s", path)


# ─────────────────────────────────────────────────────────────────────
# 모듈 레벨 worker (ProcessPoolExecutor pickle 호환)
# ─────────────────────────────────────────────────────────────────────
def _chrom_worker(args: tuple):
    chrom, chrom_len, bam_path, mode, frag_filter, min_mapq = args
    return chrom, _calc_chrom(chrom, chrom_len, bam_path,
                              mode, frag_filter, min_mapq)


# ─────────────────────────────────────────────────────────────────────
# 공개 인터페이스
# ─────────────────────────────────────────────────────────────────────
def run(
    bam_path:    str,
    out_prefix:  str,
    mode:        str             = "L",
    frag_filter: str             = "long",
    chroms:      Optional[list[str]] = None,
    min_mapq:    int             = 20,
    win_size:    int             = 1000,
    n_jobs:      int             = 4,
    save_raw_bw: bool            = False,   # raw WPS bigWig (용량 큼)
    # 하위 호환
    marker_bed:  Optional[str]   = None,
    extend:      int             = 1000,
    save_bw:     bool            = True,
    **kwargs,
) -> dict[str, str]:
    """
    BAM → bp 단위 WPS (염색체별 병렬) → bigWig + npy

    Parameters
    ----------
    mode        : "L" (뉴클레오솜) | "S" (TF footprint)
    frag_filter : "short" (≤150bp) | "long" (≥151bp) | "all"
    chroms      : 처리할 염색체 목록 (None=전체, chrM 제외)
    win_size    : adjusted WPS 의 windowMedian 블록 크기 (bp)
    save_raw_bw : True 이면 raw WPS bigWig 도 저장 (기본 False)
    """
    os.makedirs(os.path.dirname(os.path.abspath(out_prefix)) or ".", exist_ok=True)

    # 염색체 크기 로드
    chrom_sizes: dict[str, int] = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for sq in bam.header["SQ"]:
            chrom_sizes[sq["SN"]] = sq["LN"]

    if chroms is None:
        chroms = [c for c in chrom_sizes if c not in ("chrM", "MT", "M")]
    else:
        chroms = [c for c in chroms if c in chrom_sizes]

    log.info("WPS 계산 시작: mode=%s  frag=%s  chroms=%d  jobs=%d",
             mode, frag_filter, len(chroms), n_jobs)

    # 염색체별 병렬 계산
    wps_dict:  dict[str, np.ndarray] = {}
    norm_dict: dict[str, np.ndarray] = {}
    cov_dict:  dict[str, np.ndarray] = {}

    task_args = [
        (chrom, chrom_sizes[chrom], bam_path, mode, frag_filter, min_mapq)
        for chrom in chroms
    ]

    with ProcessPoolExecutor(max_workers=n_jobs) as ex:
        futures = {ex.submit(_chrom_worker, a): a[0] for a in task_args}
        for fut in as_completed(futures):
            chrom = futures[fut]
            try:
                chrom, (wps, cov) = fut.result()
                wps_dict[chrom]  = wps
                norm_dict[chrom] = _normalize(wps, cov, win_size=win_size)
                cov_dict[chrom]  = cov
                log.info("  ✓ %s  wps=[%d~%d]  norm=[%.2f~%.2f]",
                         chrom,
                         int(wps.min()), int(wps.max()),
                         float(norm_dict[chrom].min()),
                         float(norm_dict[chrom].max()))
            except Exception as exc:
                log.error("  ✗ %s: %s", chrom, exc)

    if not wps_dict:
        log.error("WPS 계산 결과 없음")
        return {}

    suffix    = f"{mode}_{frag_filter}"
    norm_bw   = f"{out_prefix}.{suffix}.WPS.norm.bw"
    cov_bw    = f"{out_prefix}.{suffix}.Coverage.bw"
    raw_bw    = f"{out_prefix}.{suffix}.WPS.raw.bw"
    npy_path  = f"{out_prefix}.{suffix}.npy"

    chrom_order = [c for c in
                   [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
                   + sorted(wps_dict.keys())
                   if c in wps_dict]

    # adjusted WPS bigWig (항상 저장)
    _write_bw(norm_bw, {c: norm_dict[c] for c in chrom_order},
              chrom_sizes, chrom_order)

    # coverage bigWig
    _write_bw(cov_bw, {c: cov_dict[c].astype(np.float32) for c in chrom_order},
              chrom_sizes, chrom_order)

    # raw WPS bigWig (선택)
    if save_raw_bw:
        _write_bw(raw_bw,
                  {c: wps_dict[c].astype(np.float32) for c in chrom_order},
                  chrom_sizes, chrom_order)

    # npy 저장 (marker profile 재사용용)
    np.save(npy_path, {
        "wps_norm": {c: norm_dict[c] for c in chrom_order},
        "coverage": {c: cov_dict[c]  for c in chrom_order},
        "mode":     mode,
        "frag":     frag_filter,
    }, allow_pickle=True)
    log.info("npy 저장: %s", npy_path)

    log.info("완료: %s.*", out_prefix)
    return {
        "wps_norm_bw": norm_bw,
        "coverage_bw": cov_bw,
        "npy": npy_path,
        "mode": mode,
        "frag_filter": frag_filter,
        "chroms": chrom_order,
    }


# ─────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────
def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="wps_compute.py",
        description="BAM → bp 단위 WPS bigWig (염색체별 병렬)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--bam",        required=True)
    p.add_argument("--out-prefix", required=True, dest="out_prefix",
                   help="출력 prefix (디렉토리/파일명)")
    p.add_argument("--mode",   choices=["L", "S"], default="L",
                   help="L=뉴클레오솜(k=120, frag 120-180bp)  S=TF(k=16, frag 35-80bp)")
    p.add_argument("--frag",   choices=["short", "long", "all"],
                   default="long", dest="frag_filter",
                   help="short(≤150bp) | long(≥151bp) | all")
    p.add_argument("--chroms", nargs="*", default=None,
                   help="처리할 염색체 (미지정=전체)")
    p.add_argument("--min-mapq",  type=int, default=20, dest="min_mapq")
    p.add_argument("--win-size",  type=int, default=1000, dest="win_size",
                   help="adjusted WPS windowMedian 블록 크기 (bp)")
    p.add_argument("--jobs",      type=int, default=4)
    p.add_argument("--save-raw-bw", action="store_true", dest="save_raw_bw",
                   help="raw WPS bigWig 도 저장 (용량 큼)")
    p.add_argument("--log-level", default="INFO", dest="log_level",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR"])
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
        win_size    = args.win_size,
        n_jobs      = args.jobs,
        save_raw_bw = args.save_raw_bw,
    )


if __name__ == "__main__":
    main()