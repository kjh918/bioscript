#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Calculate Nucleosome Accessibility Score using explicitly managed Fragment objects.
"""

import os
import sys
import argparse
import pysam
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent)) 
from src.core.fragment import Fragment

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

# ==========================================
# [Core / Utils]
# ==========================================
def parse_args():
    parser = argparse.ArgumentParser(description="Fragment 객체 기반 Nucleosome Accessibility Score 계산")
    parser.add_argument("-b", "--bam", required=True, help="정렬된 BAM/CRAM 파일 (index 필수)")
    parser.add_argument("-r", "--bed", required=True, help="입력 BED 파일 (0-based half-open)")
    parser.add_argument("-o", "--out", required=True, help="출력 TSV 경로")
    parser.add_argument("--short-max", type=int, required=True, help="NFR 최대 길이")
    parser.add_argument("--mono-min", type=int, required=True, help="Mono-nucleosome 최소 길이")
    parser.add_argument("--mono-max", type=int, required=True, help="Mono-nucleosome 최대 길이")
    parser.add_argument("--min-mapq", type=int, required=True, help="최소 MAPQ")
    
    args = parser.parse_args()

    if args.short_max <= 0 or args.mono_min <= 0 or args.mono_max <= 0:
        raise ValueError("Fragment 길이 기준은 0보다 커야 합니다.")
    if args.short_max >= args.mono_min:
        raise ValueError(f"short-max({args.short_max})는 mono-min({args.mono_min})보다 작아야 합니다.")
    if args.mono_min >= args.mono_max:
        raise ValueError(f"mono-min({args.mono_min})는 mono-max({args.mono_max})보다 작아야 합니다.")
    if not os.path.exists(args.bam):
        raise FileNotFoundError(f"BAM 파일 누락: {args.bam}")
    if not os.path.exists(args.bed):
        raise FileNotFoundError(f"BED 파일 누락: {args.bed}")

    return args


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


# ==========================================
# [Tasks / Pipeline]
# ==========================================
def extract_fragments_in_region(bam: pysam.AlignmentFile, chrom: str, start: int, end: int, min_mapq: int) -> list:
    """
    특정 구간에서 Read들을 수집하여 2쌍이 온전한 Fragment 객체 리스트를 반환합니다.
    """
    fragments = []
    # [MODIFIED] Re-introduced cache to track read pairs on the fly to fulfill the 2-read requirement
    cache = {}
    seen_qnames = set()

    try:
        iterator = bam.fetch(chrom, start, end)
    except ValueError as e:
        raise ValueError(f"BAM fetch 실패 ({chrom}:{start}-{end}): {str(e)}")

    for read in iterator:
        if read.is_unmapped or not read.is_paired or not read.is_proper_pair:
            continue
        if read.mapping_quality < min_mapq:
            continue
        if read.is_duplicate or read.is_secondary or read.is_supplementary:
            continue

        qname = read.query_name
        if qname in seen_qnames:
            continue

        # [MODIFIED] If the mate is already in the cache, pop it and construct the two-read Fragment
        if qname in cache:
            mate = cache.pop(qname)
            frag = Fragment(read, mate)
            
            if frag.overlaps(start, end):
                fragments.append(frag)
                seen_qnames.add(qname)
        else:
            cache[qname] = read

    # [MODIFIED] Process leftovers. If a read is left in the cache, its mate lies outside the fetched BED region.
    # We must explicitly fetch it using bam.mate() to construct the Fragment object.
    for qname, read in cache.items():
        try:
            mate = bam.mate(read)
            frag = Fragment(read, mate)
            
            if frag.overlaps(start, end):
                fragments.append(frag)
        except ValueError:
            # pysam cannot find the mate (e.g., unmapped, edge cases)
            continue

    return fragments


def calculate_accessibility(fragments: list, args) -> tuple:
    """Fragment 객체 리스트를 순회하며 Score 및 카운트 계산"""
    short_count = 0
    mono_count = 0
    
    # [MODIFIED] The absolute total of all fragments overlapping the BED region
    actual_region_total = len(fragments)

    for frag in fragments:
        if frag.is_short_fragment(args.short_max):
            short_count += 1
        elif frag.is_mono_nucleosome(args.mono_min, args.mono_max):
            mono_count += 1

    # Accessibility score is traditionally calculated as short / (short + mono)
    # This remains focused purely on the signal ratio, ignoring intermediate/long fragments
    target_total = short_count + mono_count
    score = (short_count / target_total) if target_total > 0 else 0.0

    # [MODIFIED] Returning actual_region_total to represent true fragment density in the region
    return short_count, mono_count, actual_region_total, score

def main():
    args = parse_args()

    header = [
        "chrom", "start", "end",
        "short_fragments",
        "mono_fragments",
        "region_total_fragments",
        "total_pass_fragments",
        "cpm_short",
        "cpm_mono",
        "cpm_total",
        "accessibility_score",
    ]

    regions = read_bed_regions(args.bed)
    print(f"Loaded {len(regions)} regions from BED.")
    
    with open_alignment(args.bam) as bam:
        total_pass_fragments = count_total_pass_fragments(bam, args.min_mapq)
    print(f"Total QC-passed fragments in alignment: {total_pass_fragments}")
    
    if total_pass_fragments == 0:
        raise ValueError("total_pass_fragments가 0입니다. BAM 또는 filtering 조건을 확인하세요.")

    with open_alignment(args.bam) as bam, open(args.out, "w") as out:
        print("\t".join(header), file=out)

        for chrom, start, end in regions:
            fragments = extract_fragments_in_region(
                bam, chrom, start, end, args.min_mapq
            )

            # [MODIFIED] Unpack the newly returned actual_region_total
            short_count, mono_count, region_total, score = calculate_accessibility(fragments, args)

            cpm_short = short_count / total_pass_fragments * 1_000_000
            cpm_mono = mono_count / total_pass_fragments * 1_000_000
            # [MODIFIED] CPM total now reflects the true overlapping fragment coverage
            cpm_total = region_total / total_pass_fragments * 1_000_000

            row = [
                chrom,
                str(start),
                str(end),
                str(short_count),
                str(mono_count),
                str(region_total),
                str(total_pass_fragments),
                f"{cpm_short:.6f}",
                f"{cpm_mono:.6f}",
                f"{cpm_total:.6f}",
                f"{score:.4f}",
            ]
            print("\t".join(row), file=out)

if __name__ == "__main__":
    main()