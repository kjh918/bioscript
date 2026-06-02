#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = mosdepth
# VERSION = 0.3.6
# THREADS = 8
# PROFILE = singularity_coverage

"""
Tool: mosdepth (0.3.6)
Profile: singularity_coverage
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="mosdepth Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='Sequence identifier used for file tracking (Default: )')
    parser.add_argument('--BamDir', required=True, default='', help='Directory containing the BAM file to be analyzed (Default: )')
    parser.add_argument('--qcResDir', required=True, default='', help='Output directory for coverage results (Default: )')
    parser.add_argument('--prefix', required=False, default='[qcResDir]/[SeqID].[InputSuffix]', help='Output prefix for all mosdepth result files (Default: [qcResDir]/[SeqID].[InputSuffix])')
    parser.add_argument('--regions_bed_gz', required=False, default='[prefix].regions.bed.gz', help='Compressed BED file containing region depth (Default: [prefix].regions.bed.gz)')
    parser.add_argument('--regions_bed_gz_csi', required=False, default='[prefix].regions.bed.gz.csi', help='Index file for the regions BED (Default: [prefix].regions.bed.gz.csi)')
    parser.add_argument('--summary_txt', required=False, default='[prefix].mosdepth.summary.txt', help='Summary of mean depths per chromosome/region (Default: [prefix].mosdepth.summary.txt)')
    parser.add_argument('--global_dist_txt', required=False, default='[prefix].global.dist.txt', help='Global distribution of coverage (Default: [prefix].global.dist.txt)')
    parser.add_argument('--InputSuffix', required=False, default='analysisReady', help='Suffix of the input BAM (e.g., analysisReady, recal, sorted) (Default: analysisReady)')
    parser.add_argument('--singularity_bin', required=False, default='singularity', help='Path to singularity executable (Default: singularity)')
    parser.add_argument('--mosdepth_bin', required=False, default='/opt/mosdepth', help='Path to mosdepth binary inside container (Default: /opt/mosdepth)')
    parser.add_argument('--sif', required=False, default='/storage/images/mosdepth-0.3.6.sif', help='mosdepth singularity image file (Default: /storage/images/mosdepth-0.3.6.sif)')
    parser.add_argument('--bind', required=False, default='/storage,/data', help='Mount paths for singularity (Default: /storage,/data)')
    parser.add_argument('--Threads', required=False, default='8', help='Number of threads for decompression and processing (Default: 8)')
    parser.add_argument('--bin', required=False, default='100000', help='Window size for quantized output (used with --by) (Default: 100000)')
    parser.add_argument('--mapq', required=False, default='20', help='Minimum mapping quality for a read to be counted (Default: 20)')
    parser.add_argument('--mosdepth_args', required=False, default='--no-per-base --fast-mode', help='Additional mosdepth arguments (Default: --no-per-base --fast-mode)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    BamDir = args.BamDir
    qcResDir = args.qcResDir
    prefix = args.prefix
    regions_bed_gz = args.regions_bed_gz
    regions_bed_gz_csi = args.regions_bed_gz_csi
    summary_txt = args.summary_txt
    global_dist_txt = args.global_dist_txt
    InputSuffix = args.InputSuffix
    singularity_bin = args.singularity_bin
    mosdepth_bin = args.mosdepth_bin
    sif = args.sif
    bind = args.bind
    Threads = args.Threads
    bin = args.bin
    mapq = args.mapq
    mosdepth_args = args.mosdepth_args

    # --- [Output Paths] ---
    if not prefix:
    prefix = f"{qcResDir}/{SeqID}.{InputSuffix}"
    if not regions_bed_gz:
    regions_bed_gz = f"{prefix}.regions.bed.gz"
    if not regions_bed_gz_csi:
    regions_bed_gz_csi = f"{prefix}.regions.bed.gz.csi"
    if not summary_txt:
    summary_txt = f"{prefix}.mosdepth.summary.txt"
    if not global_dist_txt:
    global_dist_txt = f"{prefix}.global.dist.txt"

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind} {sif} {mosdepth_bin} --threads {Threads} {mosdepth_args} --by {bin} --mapq {mapq} {prefix} {BamDir}/{SeqID}.{InputSuffix}.bam"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    if global_dist_txt:
        _tgt = os.path.dirname(global_dist_txt) if os.path.splitext(global_dist_txt)[1] else global_dist_txt
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if regions_bed_gz:
        _tgt = os.path.dirname(regions_bed_gz) if os.path.splitext(regions_bed_gz)[1] else regions_bed_gz
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if qcResDir:
        _tgt = os.path.dirname(qcResDir) if os.path.splitext(qcResDir)[1] else qcResDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if summary_txt:
        _tgt = os.path.dirname(summary_txt) if os.path.splitext(summary_txt)[1] else summary_txt
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if regions_bed_gz_csi:
        _tgt = os.path.dirname(regions_bed_gz_csi) if os.path.splitext(regions_bed_gz_csi)[1] else regions_bed_gz_csi
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if BamDir:
        _tgt = os.path.dirname(BamDir) if os.path.splitext(BamDir)[1] else BamDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if prefix:
        _tgt = os.path.dirname(prefix) if os.path.splitext(prefix)[1] else prefix
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()