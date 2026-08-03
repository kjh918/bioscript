# Accessibility calculation logic
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# [MODIFIED] New file created specifically for CNV calculations to isolate functionality
import math
import sys 
import pysam
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent)) 
from src.core.fragment import Fragment


def calculate_cnv_metrics(fragments: list, region_total: int, baseline_ref=None) -> tuple:
    """
    Calculates Read Depth (RD) and standardized Z-scores for CNV detection.
    
    Args:
        fragments: List of Fragment objects in the current bin/region.
        region_total: Total valid fragments in the region.
        baseline_ref: Baseline dictionary or model for normalization (e.g., GC correction).
        
    Returns:
        tuple: (read_depth, z_score)
    """
    
    # 1. Calculate Raw Read Depth (RD)
    # For paired-end CNV, fragment count represents the depth.
    read_depth = region_total
    
    # 2. Apply Normalization (Placeholder for GC-bias correction)
    # normalized_rd = apply_gc_correction(read_depth, gc_content)
    normalized_rd = read_depth # Placeholder
    
    # 3. Calculate Z-score based on reference
    z_score = 0.0
    if baseline_ref:
        # [MODIFIED] Placeholder for standard Z-score math: (x - mean) / std_dev
        # mean = baseline_ref.get_mean()
        # std_dev = baseline_ref.get_std()
        # z_score = (normalized_rd - mean) / std_dev
        pass
        
    return read_depth, z_score


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
