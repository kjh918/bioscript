#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = samtools
# VERSION = 1.10
# THREADS = 8
# PROFILE = bam_sort_index

"""
Tool: samtools (1.10)
Profile: bam_sort_index
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="samtools Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='Sequence identifier used for file naming (Default: )')
    parser.add_argument('--BamDir', required=True, default='', help='Directory containing the primary (unsorted) BAM file (Default: )')
    parser.add_argument('--sorted_bam', required=False, default='[BamDir]/[SeqID].[InputSuffix].bam', help='Coordinate sorted BAM file (Default: [BamDir]/[SeqID].[InputSuffix].bam)')
    parser.add_argument('--sorted_bai', required=False, default='[BamDir]/[SeqID].[OutputSuffix].bam.bai', help='BAM index file (BAI) (Default: [BamDir]/[SeqID].[OutputSuffix].bam.bai)')
    parser.add_argument('--InputSuffix', required=False, default='primary', help='Suffix of input BAM (e.g., primary, initial) (Default: primary)')
    parser.add_argument('--OutputSuffix', required=False, default='sorted', help='Suffix for output BAM (e.g., sorted, coord_sorted) (Default: sorted)')
    parser.add_argument('--samtools_bin', required=False, default='samtools', help='Path to samtools executable or binary name (Default: samtools)')
    parser.add_argument('--Threads', required=False, default='8', help='Number of threads for sorting and indexing (Default: 8)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    BamDir = args.BamDir
    sorted_bam = args.sorted_bam
    sorted_bai = args.sorted_bai
    InputSuffix = args.InputSuffix
    OutputSuffix = args.OutputSuffix
    samtools_bin = args.samtools_bin
    Threads = args.Threads

    # --- [Output Paths] ---
    if not sorted_bam:
    sorted_bam = f"{BamDir}/{SeqID}.{InputSuffix}.bam"
    if not sorted_bai:
    sorted_bai = f"{BamDir}/{SeqID}.{OutputSuffix}.bam.bai"

    # --- [Command Execution] ---
    cmd = f"{samtools_bin} sort -@ {Threads} -o {sorted_bam} {BamDir}/{SeqID}.bam && {samtools_bin} index -b -@ {Threads} {sorted_bam}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    if BamDir:
        _tgt = os.path.dirname(BamDir) if os.path.splitext(BamDir)[1] else BamDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if sorted_bam:
        _tgt = os.path.dirname(sorted_bam) if os.path.splitext(sorted_bam)[1] else sorted_bam
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if sorted_bai:
        _tgt = os.path.dirname(sorted_bai) if os.path.splitext(sorted_bai)[1] else sorted_bai
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()