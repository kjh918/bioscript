#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = samtools
# VERSION = 1.10
# THREADS = 1
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
    parser.add_argument('--InputSuffix', required=False, default='primary', help='Suffix for input BAM (e.g., sorted, coord_sorted) (Default: primary)')
    parser.add_argument('--OutputSuffix', required=False, default='sorted', help='Suffix for output BAM (e.g., sorted, coord_sorted) (Default: sorted)')
    parser.add_argument('--samtools_bin', required=False, default='samtools', help='Path to samtools executable or binary name (Default: samtools)')
    parser.add_argument('--Threads', required=False, default='8', help='Number of threads for sorting and indexing (Default: 8)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    BamDir = args.BamDir
    InputSuffix = args.InputSuffix
    OutputSuffix = args.OutputSuffix
    samtools_bin = args.samtools_bin
    Threads = args.Threads

    # --- [Output Paths] ---
    sorted_bam = f"{BamDir}/{SeqID}.{OutputSuffix}.bam"
    sorted_bai = f"{BamDir}/{SeqID}.{OutputSuffix}.bam.bai"

    # --- [Command Execution] ---
    cmd = f"{samtools_bin} sort -@ {Threads} -o {sorted_bam} {BamDir}/{SeqID}.{InputSuffix}.bam && {samtools_bin} index -b -@ {Threads} {sorted_bam}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    os.makedirs(os.path.dirname(BamDir) if '.' in os.path.basename(BamDir) else BamDir, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()