#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = mosdepth
# VERSION = 0.3.6
# THREADS = 1
# PROFILE = target_coverage_calculation

"""
Tool: mosdepth (0.3.6)
Profile: target_coverage_calculation
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="mosdepth Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='Sequence identifier used for file tracking (Default: )')
    parser.add_argument('--BamDir', required=True, default='', help='Directory containing the BAM file to be analyzed (Default: )')
    parser.add_argument('--qcResDir', required=True, default='', help='Directory for quality control results (Default: )')
    parser.add_argument('--TargetBed', required=True, default='', help='Target intervals (BED) for coverage calculation (Default: )')
    parser.add_argument('--InputSuffix', required=False, default='analysisReady', help='No description (Default: analysisReady)')
    parser.add_argument('--singularity_bin', required=False, default='singularity', help='No description (Default: singularity)')
    parser.add_argument('--mosdepth_sif', required=False, default='/storage/images/mosdepth-0.3.6.sif', help='No description (Default: /storage/images/mosdepth-0.3.6.sif)')
    parser.add_argument('--mosdepth_bin', required=False, default='/opt/mosdepth', help='No description (Default: /opt/mosdepth)')
    parser.add_argument('--bind', required=False, default='/storage,/data', help='No description (Default: /storage,/data)')
    parser.add_argument('--Threads', required=False, default='4', help='Number of threads for decompression (Default: 4)')
    parser.add_argument('--MapQ', required=False, default='20', help='Mapping quality threshold (default: 20) (Default: 20)')
    parser.add_argument('--extra_args', required=False, default='', help='Additional mosdepth flags (e.g., --no-per-base) (Default: )')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    BamDir = args.BamDir
    qcResDir = args.qcResDir
    TargetBed = args.TargetBed
    InputSuffix = args.InputSuffix
    singularity_bin = args.singularity_bin
    mosdepth_sif = args.mosdepth_sif
    mosdepth_bin = args.mosdepth_bin
    bind = args.bind
    Threads = args.Threads
    MapQ = args.MapQ
    extra_args = args.extra_args

    # --- [Output Paths] ---
    global_dist = f"{qcResDir}/{SeqID}.mosdepth.global.dist.txt"
    region_dist = f"{qcResDir}/{SeqID}.mosdepth.region.dist.txt"
    summary_txt = f"{qcResDir}/{SeqID}.mosdepth.summary.txt"
    regions_bed = f"{qcResDir}/{SeqID}.regions.bed.gz"

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind} {mosdepth_sif} {mosdepth_bin} --threads {Threads} --by {TargetBed} --mapq {MapQ} {extra_args} {qcResDir}/{SeqID} {BamDir}/{SeqID}.{InputSuffix}.bam"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    os.makedirs(os.path.dirname(BamDir) if '.' in os.path.basename(BamDir) else BamDir, exist_ok=True)
    os.makedirs(os.path.dirname(qcResDir) if '.' in os.path.basename(qcResDir) else qcResDir, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()