"""
fixtures.py
===========
테스트용 합성 데이터 생성 유틸.

생성 항목
---------
- make_synthetic_bam()  : 제어된 short/long fragment 를 담은 BAM 파일
- make_bin_bed()        : bin BED 파일 (단일 염색체)
- make_bins_raw_df()    : bins_raw parquet 대신 사용할 DataFrame
- make_bins_corrected_df(): bins_corrected DataFrame (GC/mappability 포함)
- CHROM / BIN_SIZE / CHROM_LEN : 공통 상수

모든 생성 함수는 tmp_path (pytest fixture) 를 받아 파일 경로를 반환합니다.
"""

from __future__ import annotations

import os
import struct
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import pysam

# ── 공통 상수 ────────────────────────────────────────────────────────
CHROM     = "chr1"
CHROM_LEN = 1_000_000   # 1 Mb (테스트용 축소)
BIN_SIZE  = 100_000     # 100 kb


# ─────────────────────────────────────────────────────────────────────
# BAM 생성
# ─────────────────────────────────────────────────────────────────────
def make_synthetic_bam(
    tmp_path: Path,
    n_short:  int = 200,
    n_long:   int = 200,
    short_len: int = 120,
    long_len:  int = 180,
    chrom:    str = CHROM,
    chrom_len: int = CHROM_LEN,
    seed:     int = 42,
) -> str:
    """
    제어된 short / long fragment 를 담은 sorted + indexed BAM 을 생성합니다.

    - short fragment: length = short_len (≤ 150 bp)
    - long  fragment: length = long_len  (> 150 bp)
    - 각 fragment 는 paired-end (read1 만 실질적으로 생성, template_length 설정)
    - fragment start position 은 균등 분포 (chrom_len 내)
    """
    rng      = np.random.default_rng(seed)
    bam_path = str(tmp_path / "test.bam")
    tmp_sam  = str(tmp_path / "test.sam")

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": chrom, "LN": chrom_len}],
    }

    with pysam.AlignmentFile(tmp_sam, "w", header=header) as out:
        qname_idx = 0

        def _write_read(start: int, flen: int, is_short_frag: bool) -> None:
            nonlocal qname_idx
            a                       = pysam.AlignedSegment()
            a.query_name            = f"read_{qname_idx:06d}"
            a.flag                  = 67          # paired, read1, mapped
            a.reference_id          = 0
            a.reference_start       = start
            a.mapping_quality       = 60
            a.cigar                 = [(0, min(flen, 150))]   # MATCH
            a.next_reference_id     = 0
            a.next_reference_start  = start + flen - 1
            a.template_length       = flen
            a.query_sequence        = "A" * min(flen, 150)
            a.query_qualities       = pysam.qualitystring_to_array("I" * min(flen, 150))
            out.write(a)
            qname_idx += 1

        # short fragments
        for _ in range(n_short):
            start = int(rng.integers(0, chrom_len - short_len))
            _write_read(start, short_len, is_short_frag=True)

        # long fragments
        for _ in range(n_long):
            start = int(rng.integers(0, chrom_len - long_len))
            _write_read(start, long_len, is_short_frag=False)

    # SAM → 정렬 BAM → 인덱스
    sorted_bam = str(tmp_path / "test_sorted.bam")
    pysam.sort("-o", sorted_bam, tmp_sam)
    pysam.index(sorted_bam)
    return sorted_bam


# ─────────────────────────────────────────────────────────────────────
# BED 파일 생성
# ─────────────────────────────────────────────────────────────────────
def make_bin_bed(
    tmp_path: Path,
    chrom:    str = CHROM,
    chrom_len: int = CHROM_LEN,
    bin_size: int = BIN_SIZE,
) -> str:
    """chrom 단일 염색체에 대한 bin BED 파일을 생성합니다."""
    bed_path = str(tmp_path / "bins.bed")
    rows = []
    for s in range(0, chrom_len, bin_size):
        e = min(s + bin_size, chrom_len)
        rows.append(f"{chrom}\t{s}\t{e}")
    Path(bed_path).write_text("\n".join(rows) + "\n")
    return bed_path


# ─────────────────────────────────────────────────────────────────────
# bins_raw DataFrame (parquet 없이 직접 생성)
# ─────────────────────────────────────────────────────────────────────
def make_bins_raw_df(
    n_bins:     int   = 10,
    chrom:      str   = CHROM,
    bin_size:   int   = BIN_SIZE,
    short_base: int   = 80,
    long_base:  int   = 80,
    gc_mean:    float = 0.45,
    seed:       int   = 0,
) -> pd.DataFrame:
    """
    n_bins 개의 정상 bin 을 담은 bins_raw DataFrame 을 생성합니다.
    CNV 테스트용으로 특정 bin 에 gain/loss 를 주입할 수 있습니다.
    """
    rng = np.random.default_rng(seed)
    starts = np.arange(n_bins) * bin_size

    df = pd.DataFrame({
        "chrom":            chrom,
        "start":            starts,
        "end":              starts + bin_size,
        "short_count":      rng.integers(short_base - 20, short_base + 20, n_bins).astype(np.int32),
        "long_count":       rng.integers(long_base  - 20, long_base  + 20, n_bins).astype(np.int32),
        "short_breadth":    rng.uniform(0.8, 1.0, n_bins).astype(np.float32),
        "long_breadth":     rng.uniform(0.8, 1.0, n_bins).astype(np.float32),
        "short_wps_L":      rng.normal(0, 1, n_bins).astype(np.float32),
        "short_wps_S":      rng.normal(0, 1, n_bins).astype(np.float32),
        "long_wps_L":       rng.normal(0, 1, n_bins).astype(np.float32),
        "long_wps_S":       rng.normal(0, 1, n_bins).astype(np.float32),
        "short_vaf":        rng.uniform(0.3, 0.7, n_bins).astype(np.float32),
        "long_vaf":         rng.uniform(0.3, 0.7, n_bins).astype(np.float32),
        "short_alt_depth":  rng.integers(5, 15, n_bins).astype(np.int32),
        "short_total_depth":rng.integers(20, 40, n_bins).astype(np.int32),
        "long_alt_depth":   rng.integers(5, 15, n_bins).astype(np.int32),
        "long_total_depth": rng.integers(20, 40, n_bins).astype(np.int32),
        "gc":               rng.normal(gc_mean, 0.05, n_bins).clip(0.1, 0.9).astype(np.float32),
        "mappability":      np.ones(n_bins, dtype=np.float32),
    })
    return df


def make_bins_corrected_df(
    n_bins:     int   = 10,
    chrom:      str   = CHROM,
    bin_size:   int   = BIN_SIZE,
    seed:       int   = 0,
    inject_gain_idx: Optional[int] = None,   # gain 주입할 bin index
    inject_loss_idx: Optional[int] = None,   # loss 주입할 bin index
) -> pd.DataFrame:
    """
    bins_corrected 형식의 DataFrame 을 생성합니다.
    inject_gain_idx / inject_loss_idx 로 CNV 테스트용 이상 bin 주입.
    """
    raw = make_bins_raw_df(n_bins=n_bins, chrom=chrom,
                           bin_size=bin_size, seed=seed)
    rng = np.random.default_rng(seed)

    log2_s = rng.normal(0, 0.15, n_bins).astype(np.float32)
    log2_l = rng.normal(0, 0.15, n_bins).astype(np.float32)

    if inject_gain_idx is not None:
        log2_s[inject_gain_idx] =  1.2   # trisomy 수준 gain
    if inject_loss_idx is not None:
        log2_s[inject_loss_idx] = -1.2   # loss

    raw["mappability_pass"]        = True
    raw["log2_corrected_short"]    = log2_s
    raw["log2_corrected_long"]     = log2_l
    raw["gc_correction_short"]     = np.zeros(n_bins, dtype=np.float32)
    raw["gc_correction_long"]      = np.zeros(n_bins, dtype=np.float32)
    raw["gc_fit_short"]            = np.full(n_bins, np.nan, dtype=np.float32)
    raw["gc_fit_long"]             = np.full(n_bins, np.nan, dtype=np.float32)
    raw["log2_ratio"]              = (log2_s - log2_l).astype(np.float32)
    for f in ("short","long"):
        for t in ("L","S"):
            raw[f"{f}_wps_{t}_corrected"] = raw[f"{f}_wps_{t}"].values

    return raw
