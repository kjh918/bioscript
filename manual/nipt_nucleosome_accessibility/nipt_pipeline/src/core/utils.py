# BAM and BED parsing utilities


import os
import sys
import argparse
import pysam
from pathlib import Path


def open_alignment(path: str) -> pysam.AlignmentFile:
    lower = path.lower()
    if lower.endswith(".bam"):
        return pysam.AlignmentFile(path, "rb")
    elif lower.endswith(".cram"):
        return pysam.AlignmentFile(path, "rc")
    else:
        raise ValueError("지원하지 않는 확장자입니다. (.bam / .cram)")


def read_bed_regions(bed_path: str):
    regions = []
    with open(bed_path, "r") as bed:
        for line in bed:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            if start >= end:
                raise ValueError(f"잘못된 구간(start >= end): {chrom}:{start}-{end}")
            regions.append((chrom, start, end))
    return regions


def count_total_pass_fragments(bam: pysam.AlignmentFile, min_mapq: int) -> int:
    """
    전체 BAM에서 QC 통과 paired-end fragment 수 계산.
    read1만 count해서 fragment 단위로 계산.
    """
    total = 0

    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or not read.is_paired or not read.is_proper_pair:
            continue
        if read.mate_is_unmapped:
            continue
        if read.mapping_quality < min_mapq:
            continue
        if read.is_duplicate or read.is_secondary or read.is_supplementary:
            continue
        if abs(read.template_length) == 0:
            continue
        if not read.is_read1:
            continue

        total += 1

    return total