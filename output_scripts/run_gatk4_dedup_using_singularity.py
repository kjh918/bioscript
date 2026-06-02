#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = gatk4
# VERSION = 4.4.0.0
# THREADS = 14
# PROFILE = dedup_using_singularity

"""
Tool: gatk4 (4.4.0.0)
Profile: dedup_using_singularity
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="gatk4 Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='Sequence identifier used for file naming (Default: )')
    parser.add_argument('--BamDir', required=True, default='', help='Directory containing the input BAM file (Default: )')
    parser.add_argument('--qcResDir', required=True, default='', help='Directory for duplicate metrics output (Default: )')
    parser.add_argument('--dedup_bam', required=False, default='[BamDir]/[SeqID].[InputSuffix].[OutputSuffix].bam', help='BAM file with marked/removed duplicates (Default: [BamDir]/[SeqID].[InputSuffix].[OutputSuffix].bam)')
    parser.add_argument('--dedup_bai', required=False, default='[BamDir]/[SeqID].[InputSuffix].[OutputSuffix].bam.bai', help='Index file for the deduplicated BAM (Default: [BamDir]/[SeqID].[InputSuffix].[OutputSuffix].bam.bai)')
    parser.add_argument('--metrics_txt', required=False, default='[qcResDir]/[SeqID].[InputSuffix].[OutputSuffix].metrics.txt', help='Text file containing duplication metrics (Default: [qcResDir]/[SeqID].[InputSuffix].[OutputSuffix].metrics.txt)')
    parser.add_argument('--InputSuffix', required=False, default='sorted', help='Suffix of the input BAM (e.g., sorted, primary) (Default: sorted)')
    parser.add_argument('--OutputSuffix', required=False, default='dedup', help='Suffix for the output file (e.g., dedup, md) (Default: dedup)')
    parser.add_argument('--singularity_bin', required=False, default='singularity', help='Path to singularity executable (Default: singularity)')
    parser.add_argument('--gatk_bin', required=False, default='gatk', help='GATK binary path inside container (Default: gatk)')
    parser.add_argument('--sif', required=False, default='/storage/images/gatk-4.4.0.0.sif', help='GATK singularity image (.sif) (Default: /storage/images/gatk-4.4.0.0.sif)')
    parser.add_argument('--bind', required=False, default='/storage,/data', help='Mount paths for singularity (Default: /storage,/data)')
    parser.add_argument('--Threads', required=False, default='14', help='Number of threads for ParallelGC (Default: 14)')
    parser.add_argument('--xmx_mb', required=False, default='16384', help='Maximum Java heap memory (MB) (Default: 16384)')
    parser.add_argument('--remove_seq_dups', required=False, default='true', help='Whether to remove sequencing duplicates (Default: true)')
    parser.add_argument('--create_index', required=False, default='true', help='Whether to create a BAM index for the output (Default: true)')
    parser.add_argument('--optical_duplicate_pixel_distance', required=False, default='2500', help='Pixel distance for optical duplicate detection (Default: 2500)')
    parser.add_argument('--other_md_args', required=False, default='', help='Other additional GATK MarkDuplicates arguments (Default: )')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    BamDir = args.BamDir
    qcResDir = args.qcResDir
    dedup_bam = args.dedup_bam
    dedup_bai = args.dedup_bai
    metrics_txt = args.metrics_txt
    InputSuffix = args.InputSuffix
    OutputSuffix = args.OutputSuffix
    singularity_bin = args.singularity_bin
    gatk_bin = args.gatk_bin
    sif = args.sif
    bind = args.bind
    Threads = args.Threads
    xmx_mb = args.xmx_mb
    remove_seq_dups = args.remove_seq_dups
    create_index = args.create_index
    optical_duplicate_pixel_distance = args.optical_duplicate_pixel_distance
    other_md_args = args.other_md_args

    # --- [Output Paths] ---
    if not dedup_bam:
    dedup_bam = f"{BamDir}/{SeqID}.{InputSuffix}.{OutputSuffix}.bam"
    if not dedup_bai:
    dedup_bai = f"{BamDir}/{SeqID}.{InputSuffix}.{OutputSuffix}.bam.bai"
    if not metrics_txt:
    metrics_txt = f"{qcResDir}/{SeqID}.{InputSuffix}.{OutputSuffix}.metrics.txt"

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind} {sif} {gatk_bin} MarkDuplicates --java-options '-XX:ParallelGCThreads={Threads} -Xmx{xmx_mb}m' --INPUT {BamDir}/{SeqID}.{InputSuffix}.bam --OUTPUT {dedup_bam} --METRICS_FILE {metrics_txt} --CREATE_INDEX {create_index} --REMOVE_SEQUENCING_DUPLICATES {remove_seq_dups} --OPTICAL_DUPLICATE_PIXEL_DISTANCE {optical_duplicate_pixel_distance} {other_md_args}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    if dedup_bai:
        _tgt = os.path.dirname(dedup_bai) if os.path.splitext(dedup_bai)[1] else dedup_bai
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if dedup_bam:
        _tgt = os.path.dirname(dedup_bam) if os.path.splitext(dedup_bam)[1] else dedup_bam
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if BamDir:
        _tgt = os.path.dirname(BamDir) if os.path.splitext(BamDir)[1] else BamDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if metrics_txt:
        _tgt = os.path.dirname(metrics_txt) if os.path.splitext(metrics_txt)[1] else metrics_txt
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if qcResDir:
        _tgt = os.path.dirname(qcResDir) if os.path.splitext(qcResDir)[1] else qcResDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()