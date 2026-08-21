#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = samtools_flagstats
# VERSION = 1.10
# THREADS = 1
# PROFILE = summary_flagstats

"""
Tool: samtools_flagstats (1.10)
Profile: summary_flagstats
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="samtools_flagstats Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='Sequence identifier (Default: )')
    parser.add_argument('--BamDir', required=True, default='', help='Directory containing BAM files (Default: )')
    parser.add_argument('--flag_stats_txt', required=False, default='[QcDir]/[SeqID][InputSuffix].[OutputSuffix].txt', help='Flag Stats from BAM file (Default: [QcDir]/[SeqID][InputSuffix].[OutputSuffix].txt)')
    parser.add_argument('--InputSuffix', required=False, default='.sorted', help='Suffix of the input BAM (e.g., primary, sorted, recal) (Default: .sorted)')
    parser.add_argument('--OutputSuffix', required=False, default='filtered', help='Suffix for the filtered output (Default: filtered)')
    parser.add_argument('--samtools_bin', required=False, default='samtools', help='No description (Default: samtools)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    BamDir = args.BamDir
    flag_stats_txt = args.flag_stats_txt
    InputSuffix = args.InputSuffix
    OutputSuffix = args.OutputSuffix
    samtools_bin = args.samtools_bin

    # --- [Output Paths] ---
    if not flag_stats_txt:
    flag_stats_txt = f"[QcDir]/{SeqID}{InputSuffix}.{OutputSuffix}.txt"

    # --- [Command Execution] ---
    cmd = f"{samtools_bin} flagstats {BamDir}/[SeqId]{InputSuffix}.bam > {flag_stats_txt}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    if flag_stats_txt:
        _tgt = os.path.dirname(flag_stats_txt) if os.path.splitext(flag_stats_txt)[1] else flag_stats_txt
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if BamDir:
        _tgt = os.path.dirname(BamDir) if os.path.splitext(BamDir)[1] else BamDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()