"""
bin_extractor.py
================
BAM 파일을 염색체 단위로 1회 순차 스캔하여
bin 단위 원시 지표를 수집합니다.

수집 지표 (bin 당):
  short / long raw_count     — midpoint 귀속, fragment 당 1회만 집계
  short / long breadth_ratio — coverage 마스크 기반
  gc                         — FASTA 참조 GC 비율
  mappability                — bigWig 평균

WPS 는 여기서 계산하지 않습니다.
bp 단위 WPS 는 wps_compute.py 에서 별도 수행합니다 (논문 방식).
"""

from __future__ import annotations

import bisect
import logging
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Optional

import numpy as np
import pandas as pd
import pysam

from nipt_fragmentomics.core.constants import DEFAULT_BIN_SIZE, MIN_MAPQ
from nipt_fragmentomics.core.schema import FragmentScore

log = logging.getLogger(__name__)

try:
    import pyBigWig as _pbw
    _HAS_BW = True
except ImportError:
    _HAS_BW = False


# ─────────────────────────────────────────────────────────────────────
# bin 단위 누적기
# ─────────────────────────────────────────────────────────────────────
class _BinAccumulator:
    """
    단일 bin 의 count / coverage 지표를 누적합니다.
    WPS 는 수집하지 않습니다.
    """
    __slots__ = (
        "chrom", "start", "end", "bin_len",
        "short_count", "long_count",
        "short_cov", "long_cov",
    )

    def __init__(self, chrom: str, start: int, end: int):
        self.chrom    = chrom
        self.start    = start
        self.end      = end
        self.bin_len  = end - start

        self.short_count = 0
        self.long_count  = 0
        self.short_cov   = np.zeros(self.bin_len, dtype=bool)
        self.long_cov    = np.zeros(self.bin_len, dtype=bool)

    def add_fragment(self, fs: FragmentScore) -> None:
        """FragmentScore 를 count / coverage 에 누적."""
        if fs.is_short:
            self.short_count += 1
        else:
            self.long_count += 1

        lo = max(fs.frag_start, self.start) - self.start
        hi = min(fs.frag_end,   self.end)   - self.start
        if lo < hi:
            if fs.is_short:
                self.short_cov[lo:hi] = True
            else:
                self.long_cov[lo:hi]  = True

    def to_row(self, gc: float, mappability: float) -> dict:
        bl = float(self.bin_len)
        return {
            "chrom":         self.chrom,
            "start":         self.start,
            "end":           self.end,
            "short_count":   self.short_count,
            "long_count":    self.long_count,
            "short_breadth": float(self.short_cov.sum() / bl),
            "long_breadth":  float(self.long_cov.sum()  / bl),
            "gc":            gc,
            "mappability":   mappability,
        }


# ─────────────────────────────────────────────────────────────────────
# 어노테이션 헬퍼
# ─────────────────────────────────────────────────────────────────────
def _gc_of_region(fasta: Optional[pysam.FastaFile],
                  chrom: str, start: int, end: int) -> float:
    if fasta is None:
        return float("nan")
    try:
        seq = fasta.fetch(chrom, start, end).upper()
        return (seq.count("G") + seq.count("C")) / len(seq) if seq else float("nan")
    except Exception:
        return float("nan")


def _mappability_of_region(bw, chrom: str, start: int, end: int) -> float:
    if bw is None:
        return float("nan")
    try:
        val = bw.stats(chrom, start, end, type="mean")[0]
        return float(val) if val is not None else float("nan")
    except Exception:
        return float("nan")


# ─────────────────────────────────────────────────────────────────────
# BED 로더
# ─────────────────────────────────────────────────────────────────────
def _read_bed(bed_path: str) -> pd.DataFrame:
    """
    bin.bed.gz (annotate_bin_metadata + apply_final_filters 출력) 를 읽어
    스캔 대상 bin 목록을 반환합니다.

    - .bed.gz bgzip 자동 감지
    - header 행 자동 감지
    - is_filtered 컬럼 있으면 True 인 bin 제외 (성염색체 보호)
    - gc / mappability 컬럼 있으면 보존 (재계산 생략)
    """
    import gzip

    opener = gzip.open if bed_path.endswith(".gz") else open
    with opener(bed_path, "rt") as f:
        first = f.readline().rstrip("\n")

    fields = first.split("\t")
    try:
        int(fields[1])
        has_header = False
    except (ValueError, IndexError):
        has_header = True

    bed = pd.read_csv(
        bed_path, sep="\t",
        header=0 if has_header else None,
        comment="#",
        compression="gzip" if bed_path.endswith(".gz") else None,
    )
    if not has_header:
        cols = list(bed.columns)
        cols[0], cols[1], cols[2] = "chrom", "start", "end"
        bed.columns = cols
    else:
        bed = bed.rename(columns={
            bed.columns[0]: "chrom",
            bed.columns[1]: "start",
            bed.columns[2]: "end",
        })

    before = len(bed)
    bed["start"] = pd.to_numeric(bed["start"], errors="coerce")
    bed["end"]   = pd.to_numeric(bed["end"],   errors="coerce")
    bed = bed.dropna(subset=["start", "end"]).copy()
    bed["start"] = bed["start"].astype(int)
    bed["end"]   = bed["end"].astype(int)

    if len(bed) < before:
        log.warning("BED 파싱: %d 행 제거", before - len(bed))
    if bed.empty:
        raise ValueError(f"BED 파일에서 유효한 행을 읽지 못했습니다: {bed_path}")

    if "is_filtered" in bed.columns:
        sex_mask  = bed["chrom"].isin({"chrX", "chrY", "X", "Y"})
        keep_mask = (~bed["is_filtered"].astype(bool)) | sex_mask
        n_excl    = (~keep_mask).sum()
        bed       = bed[keep_mask].copy()
        log.info("is_filtered 적용: %d bins 제외", n_excl)

    log.info("BED 로드: %d bins  [%s]", len(bed), bed_path)
    return bed


# ─────────────────────────────────────────────────────────────────────
# 염색체 단위 스캔 (worker)
# ─────────────────────────────────────────────────────────────────────
def _scan_chrom(
    chrom:      str,
    bins:       list[tuple[int, int]],
    bam_path:   str,
    fasta_path: Optional[str],
    bw_path:    Optional[str],
    min_mapq:   int,
    gc_precomp:  Optional[list[float]] = None,
    map_precomp: Optional[list[float]] = None,
) -> list[dict]:
    """
    chrom 의 bins 를 1회 순차 스캔.
    각 fragment 는 midpoint 가 속한 bin 에만 귀속 (bisect O(log n)).
    gc_precomp / map_precomp 가 있으면 FASTA / bw 재계산 생략.
    """
    if not bins:
        return []

    bin_starts   = [s for s, _ in bins]
    accumulators = [_BinAccumulator(chrom, s, e) for s, e in bins]

    if gc_precomp is not None:
        gc_vals = gc_precomp
    else:
        fasta_h = pysam.FastaFile(fasta_path) if fasta_path else None
        gc_vals = [_gc_of_region(fasta_h, chrom, s, e) for s, e in bins]
        if fasta_h: fasta_h.close()

    if map_precomp is not None:
        map_vals = map_precomp
    else:
        bw_h = (_pbw.open(bw_path) if (_HAS_BW and bw_path) else None)
        map_vals = [_mappability_of_region(bw_h, chrom, s, e) for s, e in bins]
        if bw_h: bw_h.close()

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

            fs = FragmentScore.from_read(read)
            if fs is None:
                continue

            idx = bisect.bisect_right(bin_starts, fs.midpoint) - 1
            if idx < 0 or idx >= len(bins):
                continue
            s, e = bins[idx]
            if not (s <= fs.midpoint < e):
                continue

            accumulators[idx].add_fragment(fs)

    return [acc.to_row(gc_vals[i], map_vals[i])
            for i, acc in enumerate(accumulators)]


# ─────────────────────────────────────────────────────────────────────
# 공개 인터페이스
# ─────────────────────────────────────────────────────────────────────
def run(
    bam_path:   str,
    out_path:   str,
    bed_path:   Optional[str] = None,
    fasta_path: Optional[str] = None,
    bw_path:    Optional[str] = None,
    bin_size:   int           = DEFAULT_BIN_SIZE,
    min_mapq:   int           = MIN_MAPQ,
    n_jobs:     int           = 4,
    # 하위 호환용 (무시됨)
    vcf_path:   Optional[str] = None,
    min_baseq:  int           = 20,
) -> pd.DataFrame:
    """BAM → bins_raw.parquet (count / breadth / gc / mappability)"""
    if vcf_path is not None:
        log.warning("vcf_path 는 더 이상 사용되지 않습니다. 무시합니다.")

    # bin 목록 구성
    if bed_path and os.path.exists(bed_path):
        bed = _read_bed(bed_path)
    else:
        log.info("BED 없음 — BAM 헤더에서 자동 bin 생성 (bin_size=%d)", bin_size)
        rows = []
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for sq in bam.header["SQ"]:
                chrom, length = sq["SN"], sq["LN"]
                for s in range(0, length, bin_size):
                    rows.append({"chrom": chrom, "start": s,
                                 "end": min(s + bin_size, length)})
        bed = pd.DataFrame(rows)

    bed = bed.sort_values(["chrom", "start"]).reset_index(drop=True)

    # BED 에 gc / mappability 있으면 재계산 생략
    has_gc  = "gc"          in bed.columns
    has_map = "mappability" in bed.columns
    if has_gc:  log.info("BED 에 gc 컬럼 존재 → FASTA 재계산 생략")
    if has_map: log.info("BED 에 mappability 컬럼 존재 → bigWig 재계산 생략")

    chrom_data: dict[str, dict] = {}
    for chrom, g in bed.groupby("chrom"):
        chrom_data[chrom] = {
            "bins": list(zip(g["start"], g["end"])),
            "gc":   g["gc"].tolist()          if has_gc  else None,
            "map":  g["mappability"].tolist() if has_map else None,
        }

    log.info("총 %d bins, %d chroms", len(bed), len(chrom_data))

    all_rows: list[dict] = []
    with ProcessPoolExecutor(max_workers=n_jobs) as ex:
        futures = {
            ex.submit(
                _scan_chrom, chrom,
                chrom_data[chrom]["bins"],
                bam_path,
                None if has_gc  else fasta_path,
                None if has_map else bw_path,
                min_mapq,
                chrom_data[chrom]["gc"],
                chrom_data[chrom]["map"],
            ): chrom
            for chrom in chrom_data
        }
        for fut in as_completed(futures):
            chrom = futures[fut]
            try:
                rows = fut.result()
                all_rows.extend(rows)
                log.info("  ✓ %s (%d bins)", chrom, len(rows))
            except Exception as exc:
                log.error("  ✗ %s: %s", chrom, exc)

    df = pd.DataFrame(all_rows)

    for c in ["short_count", "long_count"]:
        if c in df: df[c] = df[c].astype("int32")
    for c in ["short_breadth", "long_breadth", "gc", "mappability"]:
        if c in df: df[c] = df[c].astype("float32")

    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    df.to_parquet(out_path, engine="pyarrow", index=False)
    log.info("bins_raw 저장: %s (%d rows)", out_path, len(df))
    return df
