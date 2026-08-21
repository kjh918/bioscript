"""
cnv/binning.py

BAM → bin-level fragment count TSV  (chromosome 단위 multiprocessing)

- 리드 중앙(mid = start + tlen//2)이 bin 안에 들어오는 경우만 count
- duplicate / read2 / unmapped 제외
- short(<cutoff bp) / long(>=cutoff bp) / total 분류
- mappability bigWig hard mask (bin 평균 < threshold 제외)
- GC / mappability / count 세 단계 모두 chromosome 단위 병렬 처리

출력 컬럼:
  chrom, bin_start, bin_end, gc, mappability, total, short, long
"""
from __future__ import annotations

import logging
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Optional

import numpy as np
import pysam
import pyBigWig

from .utils import HG38_CHROM_SIZES, STANDARD_CHROMS, write_tsv

logger = logging.getLogger(__name__)

# ── 상수 ──────────────────────────────────────────────────────────────
_MAX_FRAG = 1000

FIELDNAMES = ["chrom", "bin_start", "bin_end", "gc", "mappability",
              "total", "short", "long"]


# ── chromosome 단위 worker (module-level, pickle 호환) ────────────────
def _chrom_worker(args: tuple) -> list[dict]:
    """
    단일 chromosome 처리 → bin dict 리스트 반환 (masked bin 제외)

    args = (chrom, chrom_len, bam_path, bw_path, fasta_path,
            bin_size, short_long_cutoff, mapq_min, mappability_threshold)
    """
    (chrom, chrom_len, bam_path, bw_path, fasta_path,
     bin_size, short_long_cutoff, mapq_min, mappability_threshold) = args

    n_bins = (chrom_len + bin_size - 1) // bin_size

    # numpy array로 bin 데이터 관리
    bin_starts      = np.arange(n_bins, dtype=np.int32) * bin_size
    bin_ends        = np.minimum(bin_starts + bin_size, chrom_len).astype(np.int32)
    gc_arr          = np.zeros(n_bins, dtype=np.float32)
    mappability_arr = np.zeros(n_bins, dtype=np.float32)
    masked_arr      = np.ones(n_bins,  dtype=bool)       # 기본 전체 mask
    total_arr       = np.zeros(n_bins, dtype=np.int32)
    short_arr       = np.zeros(n_bins, dtype=np.int32)
    long_arr        = np.zeros(n_bins, dtype=np.int32)

    # ── step 1: mappability ───────────────────────────────────────────
    try:
        bw = pyBigWig.open(bw_path)
        bw_chroms = bw.chroms()
        if chrom in bw_chroms:
            for i in range(n_bins):
                vals = bw.stats(chrom, int(bin_starts[i]), int(bin_ends[i]), type="mean")
                m = vals[0] if vals and vals[0] is not None else 0.0
                mappability_arr[i] = m
                masked_arr[i]      = m < mappability_threshold
        else:
            masked_arr[:] = True   # chrom 없으면 전체 mask
        bw.close()
    except Exception as e:
        logger.warning("[%s] mappability 오류: %s", chrom, e)
        masked_arr[:] = True

    # ── step 2: GC content ────────────────────────────────────────────
    try:
        fa = pysam.FastaFile(fasta_path)
        if chrom in fa.references:
            for i in range(n_bins):
                if masked_arr[i]:
                    continue   # masked bin은 GC 계산 스킵
                seq = fa.fetch(chrom, int(bin_starts[i]), int(bin_ends[i])).upper()
                if seq:
                    gc_arr[i] = sum(1 for c in seq if c in ("G", "C")) / len(seq)
        fa.close()
    except Exception as e:
        logger.warning("[%s] GC 계산 오류: %s", chrom, e)

    # ── step 3: read count ────────────────────────────────────────────
    try:
        bam = pysam.AlignmentFile(bam_path, "rb", threads=1)
        for read in bam.fetch(chrom):
            # 기본 필터
            if (read.is_unmapped or read.mate_is_unmapped or
                    read.is_duplicate or not read.is_paired or
                    not read.is_proper_pair or read.is_read2):
                continue
            if read.mapping_quality < mapq_min:
                continue

            tlen = abs(read.template_length)
            if tlen == 0 or tlen > _MAX_FRAG:
                continue

            mid     = read.reference_start + tlen // 2
            bin_idx = mid // bin_size
            if bin_idx >= n_bins:
                continue
            if masked_arr[bin_idx]:
                continue

            total_arr[bin_idx] += 1
            if tlen < short_long_cutoff:
                short_arr[bin_idx] += 1
            else:
                long_arr[bin_idx] += 1
        bam.close()
    except Exception as e:
        logger.warning("[%s] BAM count 오류: %s", chrom, e)

    # ── 결과 조립 (masked bin 제외) ───────────────────────────────────
    rows: list[dict] = []
    for i in range(n_bins):
        if masked_arr[i]:
            continue
        rows.append({
            "chrom":       chrom,
            "bin_start":   int(bin_starts[i]),
            "bin_end":     int(bin_ends[i]),
            "gc":          f"{gc_arr[i]:.4f}",
            "mappability": f"{mappability_arr[i]:.4f}",
            "total":       int(total_arr[i]),
            "short":       int(short_arr[i]),
            "long":        int(long_arr[i]),
        })
    return rows


# ── 공개 인터페이스 ────────────────────────────────────────────────────
def run_binning(
    bam_path:              str,
    bw_path:               str,
    fasta_path:            str,
    output_path:           str,
    bin_size:              int   = 100_000,
    short_long_cutoff:     int   = 150,
    mapq_min:              int   = 30,
    mappability_threshold: float = 0.9,
    chroms:                Optional[list[str]] = None,
    n_jobs:                int   = 4,
) -> None:
    """
    BAM → bin count TSV (chromosome 단위 multiprocessing)

    Parameters
    ----------
    bam_path              : 입력 BAM (sorted, indexed)
    bw_path               : mappability bigWig
    fasta_path            : hg38 reference FASTA
    output_path           : 출력 TSV 경로
    bin_size              : bin 크기 (default 100kb)
    short_long_cutoff     : short/long 분류 기준 (bp)
    mapq_min              : 최소 MAPQ
    mappability_threshold : bin 평균 mappability 하한
    chroms                : 처리할 chrom 목록 (None이면 전체 standard)
    n_jobs                : 병렬 프로세스 수
    """
    target = chroms or STANDARD_CHROMS

    # BAM 헤더에서 chrom 크기 확인 (헤더에 없는 chrom 제외)
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        bam_chroms = {sq["SN"]: sq["LN"] for sq in bam.header.get("SQ", [])}

    target = [c for c in target if c in bam_chroms]
    if not target:
        raise ValueError("BAM 헤더에 처리 가능한 chrom 없음")

    logger.info("binning 시작: %d chroms, bin=%dbp, jobs=%d",
                len(target), bin_size, n_jobs)

    task_args = [
        (chrom,
         bam_chroms.get(chrom, HG38_CHROM_SIZES.get(chrom, 0)),
         bam_path, bw_path, fasta_path,
         bin_size, short_long_cutoff, mapq_min, mappability_threshold)
        for chrom in target
    ]

    # chrom 순서 유지를 위해 결과를 dict로 수집
    chrom_rows: dict[str, list[dict]] = {}

    with ProcessPoolExecutor(max_workers=n_jobs) as ex:
        futures = {ex.submit(_chrom_worker, a): a[0] for a in task_args}
        for fut in as_completed(futures):
            chrom = futures[fut]
            try:
                rows = fut.result()
                chrom_rows[chrom] = rows
                logger.info("  ✓ %s: %d bins", chrom, len(rows))
            except Exception as exc:
                logger.exception("  ✗ %s: %s", chrom, exc)
                chrom_rows[chrom] = []

    # chrom 순서대로 합치기
    all_rows: list[dict] = []
    for chrom in target:
        all_rows.extend(chrom_rows.get(chrom, []))

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    write_tsv(all_rows, output_path, FIELDNAMES)
    logger.info("binning 완료 → %s (%d bins)", output_path, len(all_rows))