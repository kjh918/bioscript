#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = samtools_filter_index
# VERSION = 1.10
# THREADS = 8
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
    parser.add_argument('--filtered_bam', required=False, default='[BamDir]/[SeqID].[InputSuffix].[OutputSuffix].bam', help='Filtered BAM file (Default: [BamDir]/[SeqID].[InputSuffix].[OutputSuffix].bam)')
    parser.add_argument('--filtered_bai', required=False, default='[BamDir]/[SeqID].[InputSuffix].[OutputSuffix].bam.bai', help='No description (Default: [BamDir]/[SeqID].[InputSuffix].[OutputSuffix].bam.bai)')
    parser.add_argument('--analysis_ready_bam', required=False, default='[BamDir]/[SeqID].analysisReady.bam', help='No description (Default: [BamDir]/[SeqID].analysisReady.bam)')
    parser.add_argument('--analysis_ready_bai', required=False, default='[BamDir]/[SeqID].analysisReady.bam.bai', help='No description (Default: [BamDir]/[SeqID].analysisReady.bam.bai)')
    parser.add_argument('--InputSuffix', required=False, default='recal', help='Suffix of the input BAM (e.g., primary, sorted, recal) (Default: recal)')
    parser.add_argument('--OutputSuffix', required=False, default='filtered', help='Suffix for the filtered output (Default: filtered)')
    parser.add_argument('--samtools_bin', required=False, default='samtools', help='No description (Default: samtools)')
    parser.add_argument('--ln_bin', required=False, default='ln', help='No description (Default: ln)')
    parser.add_argument('--Threads', required=False, default='8', help='No description (Default: 8)')
    parser.add_argument('--min_mapq', required=False, default='20', help='No description (Default: 20)')
    parser.add_argument('--include_flag', required=False, default='0x2', help='No description (Default: 0x2)')
    parser.add_argument('--exclude_flag', required=False, default='0x104', help='No description (Default: 0x104)')
    parser.add_argument('--expr_nm', required=False, default='[NM] < 12', help='No description (Default: [NM] < 12)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    BamDir = args.BamDir
    filtered_bam = args.filtered_bam
    filtered_bai = args.filtered_bai
    analysis_ready_bam = args.analysis_ready_bam
    analysis_ready_bai = args.analysis_ready_bai
    InputSuffix = args.InputSuffix
    OutputSuffix = args.OutputSuffix
    samtools_bin = args.samtools_bin
    ln_bin = args.ln_bin
    Threads = args.Threads
    min_mapq = args.min_mapq
    include_flag = args.include_flag
    exclude_flag = args.exclude_flag
    expr_nm = args.expr_nm

    # --- [Output Paths] ---
    if not filtered_bam:
    filtered_bam = f"{BamDir}/{SeqID}.{InputSuffix}.{OutputSuffix}.bam"
    if not filtered_bai:
    filtered_bai = f"{BamDir}/{SeqID}.{InputSuffix}.{OutputSuffix}.bam.bai"
    if not analysis_ready_bam:
    analysis_ready_bam = f"{BamDir}/{SeqID}.analysisReady.bam"
    if not analysis_ready_bai:
    analysis_ready_bai = f"{BamDir}/{SeqID}.analysisReady.bam.bai"

    # --- [Command Execution] ---
    cmd = f"{samtools_bin} view -b -h -q {min_mapq} -f {include_flag} -F {exclude_flag} -e '{expr_nm}' --threads {Threads} {BamDir}/{SeqID}.{InputSuffix}.bam > {filtered_bam} && {samtools_bin} index --threads {Threads} {filtered_bam} && {ln_bin} -Tsf {filtered_bam} {analysis_ready_bam} && {ln_bin} -Tsf {filtered_bai} {analysis_ready_bai}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    if filtered_bam:
        _tgt = os.path.dirname(filtered_bam) if os.path.splitext(filtered_bam)[1] else filtered_bam
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if analysis_ready_bam:
        _tgt = os.path.dirname(analysis_ready_bam) if os.path.splitext(analysis_ready_bam)[1] else analysis_ready_bam
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if BamDir:
        _tgt = os.path.dirname(BamDir) if os.path.splitext(BamDir)[1] else BamDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if analysis_ready_bai:
        _tgt = os.path.dirname(analysis_ready_bai) if os.path.splitext(analysis_ready_bai)[1] else analysis_ready_bai
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if filtered_bai:
        _tgt = os.path.dirname(filtered_bai) if os.path.splitext(filtered_bai)[1] else filtered_bai
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()