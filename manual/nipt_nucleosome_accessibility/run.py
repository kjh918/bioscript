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

    # 파라미터 논리 검증 (Default 없이 엄격한 에러 처리)
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
    cache = {}
    fragments = []

    try:
        iterator = bam.fetch(chrom, start, end)
    except ValueError as e:
        raise ValueError(f"BAM fetch 실패 ({chrom}:{start}-{end}): {str(e)}")

    for read in iterator:
        # 기본 품질 및 Pair 조건 필터링
        if read.is_unmapped or not read.is_paired or not read.is_proper_pair:
            continue
        if read.mapping_quality < min_mapq:
            continue
        if read.is_duplicate or read.is_secondary or read.is_supplementary:
            continue
        if abs(read.template_length) == 0:
            continue

        qname = read.query_name
        
        # [MODIFIED] Cache를 활용하여 동일 구간 내의 Mate를 찾아 Fragment 생성
        if qname in cache:
            mate = cache.pop(qname)
            fragments.append(Fragment(read, mate))
        else:
            cache[qname] = read

    # [MODIFIED] Cache에 남은 리드 = Mate가 BED 구간 밖에 있는 경우 (bam.mate로 명시적 검색)
    for qname, read in cache.items():
        try:
            mate = bam.mate(read)
            fragments.append(Fragment(read, mate))
        except ValueError:
            # pysam이 mate를 찾을 수 없는 경우 (unmapped, 다른 contig 등) 무시
            continue

    return fragments


def calculate_accessibility(fragments: list, args) -> tuple:
    """Fragment 객체 리스트를 순회하며 Score 계산"""
    short_count = 0
    mono_count = 0

    for frag in fragments:
        # 로직을 도메인 객체(Fragment)로 위임
        if frag.is_short_fragment(args.short_max):
            short_count += 1
        elif frag.is_mono_nucleosome(args.mono_min, args.mono_max):
            mono_count += 1

    total = short_count + mono_count
    score = (short_count / total) if total > 0 else 0.0

    return short_count, mono_count, score


def main():
    args = parse_args()
    header = ["chrom", "start", "end", "short_fragments", "mono_fragments", "accessibility_score"]
    regions = read_bed_regions(args.bed)
    
    with open_alignment(args.bam) as bam, open(args.out, "w") as out:
        print("\t".join(header), file=out)
        
        for chrom, start, end in regions:
            # 1. Pipeline: Read -> Fragment 생성
            fragments = extract_fragments_in_region(bam, chrom, start, end, args.min_mapq)
            
            # 2. Pipeline: Fragment -> Score 계산
            short_count, mono_count, score = calculate_accessibility(fragments, args)
            
            row = [chrom, str(start), str(end), str(short_count), str(mono_count), f"{score:.4f}"]
            print("\t".join(row), file=out)

if __name__ == "__main__":
    main()