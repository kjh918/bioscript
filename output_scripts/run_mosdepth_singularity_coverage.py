#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = mosdepth
# VERSION = 0.3.6
# THREADS = 1
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
    prefix = f"{qcResDir}/{SeqID}.{InputSuffix}"
    regions_bed_gz = f"{prefix}.regions.bed.gz"
    regions_bed_gz_csi = f"{prefix}.regions.bed.gz.csi"
    summary_txt = f"{prefix}.mosdepth.summary.txt"
    global_dist_txt = f"{prefix}.global.dist.txt"

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind} {sif} {mosdepth_bin} --threads {Threads} {mosdepth_args} --by {bin} --mapq {mapq} {prefix} {BamDir}/{SeqID}.{InputSuffix}.bam"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    os.makedirs(os.path.dirname(BamDir) if '.' in os.path.basename(BamDir) else BamDir, exist_ok=True)
    os.makedirs(os.path.dirname(qcResDir) if '.' in os.path.basename(qcResDir) else qcResDir, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()