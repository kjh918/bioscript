#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = samtools_filter_index
# VERSION = 1.10
# THREADS = 1
# PROFILE = flexible_bam_filter

"""
Tool: samtools_filter_index (1.10)
Profile: flexible_bam_filter
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="samtools_filter_index Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='Sequence identifier (Default: )')
    parser.add_argument('--BamDir', required=True, default='', help='Directory containing BAM files (Default: )')
    parser.add_argument('--InputSuffix', required=False, default='recal', help='Suffix of the input BAM (e.g., primary, sorted, recal) (Default: recal)')
    parser.add_argument('--OutputSuffix', required=False, default='filtered', help='Suffix for the filtered output (Default: filtered)')
    parser.add_argument('--samtools_bin', required=False, default='samtools', help='No description (Default: samtools)')
    parser.add_argument('--Threads', required=False, default='8', help='No description (Default: 8)')
    parser.add_argument('--min_mapq', required=False, default='20', help='No description (Default: 20)')
    parser.add_argument('--include_flag', required=False, default='0x2', help='No descriptison (Default: 0x2)')
    parser.add_argument('--exclude_flag', required=False, default='0x104', help='No description (Default: 0x104)')
    parser.add_argument('--expr_nm', required=False, default='[NM] < 12', help='No description (Default: [NM] < 12)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    BamDir = args.BamDir
    InputSuffix = args.InputSuffix
    OutputSuffix = args.OutputSuffix
    samtools_bin = args.samtools_bin
    Threads = args.Threads
    min_mapq = args.min_mapq
    include_flag = args.include_flag
    exclude_flag = args.exclude_flag
    expr_nm = args.expr_nm

    # --- [Output Paths] ---
    filtered_bam = f"{BamDir}/{SeqID}.{InputSuffix}.{OutputSuffix}.bam"
    filtered_bai = f"{BamDir}/{SeqID}.{InputSuffix}.{OutputSuffix}.bam.bai"
    analysis_ready_bam = f"{BamDir}/{SeqID}.analysisReady.bam"
    analysis_ready_bai = f"{BamDir}/{SeqID}.analysisReady.bam.bai"

    # --- [Command Execution] ---
    cmd = f"{samtools_bin} view -b -h -q {min_mapq} -f {include_flag} -F {exclude_flag} -e '{expr_nm}' --threads {Threads} {BamDir}/{SeqID}.{InputSuffix}.bam > {filtered_bam} && {samtools_bin} index --threads {Threads} {filtered_bam} && ln -Tsf {filtered_bam} {analysis_ready_bam} && ln -Tsf {filtered_bai} {analysis_ready_bai}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    os.makedirs(os.path.dirname(BamDir) if '.' in os.path.basename(BamDir) else BamDir, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()