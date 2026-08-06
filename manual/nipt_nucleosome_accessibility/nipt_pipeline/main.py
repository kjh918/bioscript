# Main CLI entry point
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent)) 
from src.core.utils import open_alignment, read_bed_regions, count_total_pass_fragments
from src.accessibility.calculator import extract_fragments_in_region, calculate_accessibility
# [MODIFIED] Imported the new CNV calculator module
#from src.cnv.calculator import calculate_cnv_metrics 

def parse_args():
    parser = argparse.ArgumentParser(description="Multi-feature NIPT Pipeline: Accessibility & CNV")
    parser.add_argument("-b", "--bam", required=True, help="Sorted BAM/CRAM file (index required)")
    parser.add_argument("-r", "--bed", required=True, help="Input BED file for binning/regions")
    parser.add_argument("-o", "--out", required=True, help="Output TSV path")
    
    # Accessibility params
    parser.add_argument("--short-max", type=int, default=150, help="Max length for NFR")
    parser.add_argument("--mono-min", type=int, default=151, help="Min length for Mono-nucleosome")
    parser.add_argument("--mono-max", type=int, default=220, help="Max length for Mono-nucleosome")
    parser.add_argument("--min-mapq", type=int, default=20, help="Minimum MAPQ")
    
    # [MODIFIED] Added arguments specifically for CNV normalization and scoring
    parser.add_argument("--enable-cnv", action="store_true", help="Enable CNV calculation module")
    parser.add_argument("--ref-baseline", type=str, default=None, help="Optional: Reference baseline for Z-score")
    
    args = parser.parse_args()
    # (Validation logic omitted for brevity, identical to your original script)
    return args

def main():
    args = parse_args()
    regions = read_bed_regions(args.bed)
    
    with open_alignment(args.bam) as bam:
        total_pass_fragments = count_total_pass_fragments(bam, args.min_mapq)

    # [MODIFIED] Updated header to include CNV specific metrics (Read Depth, Z-score)
    header = [
        "chrom", "start", "end",
        "short_fragments", "mono_fragments", "region_total_fragments",
        "accessibility_score", "cpm_total",
        "cnv_read_depth", "cnv_z_score" 
    ]

    with open_alignment(args.bam) as bam, open(args.out, "w") as out:
        print("\t".join(header), file=out)

        for chrom, start, end in regions:
            # 1. Accessibility Feature Extraction
            fragments = extract_fragments_in_region(bam, chrom, start, end, args.min_mapq)
            short_count, mono_count, region_total, acc_score = calculate_accessibility(fragments, args)
            cpm_total = (region_total / total_pass_fragments) * 1_000_000 if total_pass_fragments else 0

            # 2. CNV Feature Extraction
            cnv_rd = 0
            cnv_z_score = 0.0
            
            # [MODIFIED] Trigger CNV module if requested
            if args.enable_cnv:
                # We can reuse the extracted fragments or query the BAM specifically for depth
                cnv_rd, cnv_z_score = calculate_cnv_metrics(
                    fragments=fragments, 
                    region_total=region_total, 
                    baseline_ref=args.ref_baseline
                )

            row = [
                chrom, str(start), str(end),
                str(short_count), str(mono_count), str(region_total),
                f"{acc_score:.4f}", f"{cpm_total:.6f}",
                # [MODIFIED] Append CNV outputs to the row
                str(cnv_rd), f"{cnv_z_score:.4f}"
            ]
            print("\t".join(row), file=out)

if __name__ == "__main__":
    main()