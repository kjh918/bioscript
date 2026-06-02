#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = mosdepth
# VERSION = 0.3.6
# THREADS = 4
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
    parser.add_argument('--global_dist', required=False, default='[qcResDir]/[SeqID].mosdepth.global.dist.txt', help='No description (Default: [qcResDir]/[SeqID].mosdepth.global.dist.txt)')
    parser.add_argument('--region_dist', required=False, default='[qcResDir]/[SeqID].mosdepth.region.dist.txt', help='No description (Default: [qcResDir]/[SeqID].mosdepth.region.dist.txt)')
    parser.add_argument('--summary_txt', required=False, default='[qcResDir]/[SeqID].mosdepth.summary.txt', help='No description (Default: [qcResDir]/[SeqID].mosdepth.summary.txt)')
    parser.add_argument('--regions_bed', required=False, default='[qcResDir]/[SeqID].regions.bed.gz', help='No description (Default: [qcResDir]/[SeqID].regions.bed.gz)')
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
    global_dist = args.global_dist
    region_dist = args.region_dist
    summary_txt = args.summary_txt
    regions_bed = args.regions_bed
    InputSuffix = args.InputSuffix
    singularity_bin = args.singularity_bin
    mosdepth_sif = args.mosdepth_sif
    mosdepth_bin = args.mosdepth_bin
    bind = args.bind
    Threads = args.Threads
    MapQ = args.MapQ
    extra_args = args.extra_args

    # --- [Output Paths] ---
    if not global_dist:
    global_dist = f"{qcResDir}/{SeqID}.mosdepth.global.dist.txt"
    if not region_dist:
    region_dist = f"{qcResDir}/{SeqID}.mosdepth.region.dist.txt"
    if not summary_txt:
    summary_txt = f"{qcResDir}/{SeqID}.mosdepth.summary.txt"
    if not regions_bed:
    regions_bed = f"{qcResDir}/{SeqID}.regions.bed.gz"

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind} {mosdepth_sif} {mosdepth_bin} --threads {Threads} --by {TargetBed} --mapq {MapQ} {extra_args} {qcResDir}/{SeqID} {BamDir}/{SeqID}.{InputSuffix}.bam"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    if qcResDir:
        _tgt = os.path.dirname(qcResDir) if os.path.splitext(qcResDir)[1] else qcResDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if region_dist:
        _tgt = os.path.dirname(region_dist) if os.path.splitext(region_dist)[1] else region_dist
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if summary_txt:
        _tgt = os.path.dirname(summary_txt) if os.path.splitext(summary_txt)[1] else summary_txt
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if BamDir:
        _tgt = os.path.dirname(BamDir) if os.path.splitext(BamDir)[1] else BamDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if regions_bed:
        _tgt = os.path.dirname(regions_bed) if os.path.splitext(regions_bed)[1] else regions_bed
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if global_dist:
        _tgt = os.path.dirname(global_dist) if os.path.splitext(global_dist)[1] else global_dist
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()