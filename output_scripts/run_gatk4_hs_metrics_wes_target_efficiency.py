#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = gatk4_hs_metrics
# VERSION = 4.4.0.0
# THREADS = 14
# PROFILE = wes_target_efficiency

"""
Tool: gatk4_hs_metrics (4.4.0.0)
Profile: wes_target_efficiency
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="gatk4_hs_metrics Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='No description (Default: )')
    parser.add_argument('--BamDir', required=True, default='', help='No description (Default: )')
    parser.add_argument('--qcResDir', required=True, default='', help='No description (Default: )')
    parser.add_argument('--ReferenceFasta', required=True, default='', help='No description (Default: )')
    parser.add_argument('--BaitIntervals', required=True, default='', help='Target capture kit bait intervals (Picard-style) (Default: )')
    parser.add_argument('--TargetIntervals', required=True, default='', help='Target capture kit target intervals (Picard-style) (Default: )')
    parser.add_argument('--InputSuffix', required=True, default='analysisReady', help='Target BAM suffix (Default: analysisReady)')
    parser.add_argument('--hs_metrics', required=False, default='[qcResDir]/[SeqID].[InputSuffix].HS.metrics.txt', help='Hybrid Selection (HS) metrics report (Default: [qcResDir]/[SeqID].[InputSuffix].HS.metrics.txt)')
    parser.add_argument('--singularity_bin', required=False, default='singularity', help='No description (Default: singularity)')
    parser.add_argument('--gatk4_sif', required=False, default='/storage/images/gatk-4.4.0.0.sif', help='No description (Default: /storage/images/gatk-4.4.0.0.sif)')
    parser.add_argument('--bind', required=False, default='/storage,/data', help='No description (Default: /storage,/data)')
    parser.add_argument('--xmx_mb', required=False, default='16384', help='No description (Default: 16384)')
    parser.add_argument('--Threads', required=False, default='14', help='No description (Default: 14)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    BamDir = args.BamDir
    qcResDir = args.qcResDir
    ReferenceFasta = args.ReferenceFasta
    BaitIntervals = args.BaitIntervals
    TargetIntervals = args.TargetIntervals
    InputSuffix = args.InputSuffix
    hs_metrics = args.hs_metrics
    singularity_bin = args.singularity_bin
    gatk4_sif = args.gatk4_sif
    bind = args.bind
    xmx_mb = args.xmx_mb
    Threads = args.Threads

    # --- [Output Paths] ---
    if not hs_metrics:
    hs_metrics = f"{qcResDir}/{SeqID}.{InputSuffix}.HS.metrics.txt"

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind} {gatk4_sif} gatk CollectHsMetrics --java-options '-XX:ParallelGCThreads={Threads} -Xmx{xmx_mb}m' --INPUT {BamDir}/{SeqID}.{InputSuffix}.bam --OUTPUT {hs_metrics} --REFERENCE_SEQUENCE {ReferenceFasta} --BAIT_INTERVALS {BaitIntervals} --TARGET_INTERVALS {TargetIntervals}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    if hs_metrics:
        _tgt = os.path.dirname(hs_metrics) if os.path.splitext(hs_metrics)[1] else hs_metrics
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if BamDir:
        _tgt = os.path.dirname(BamDir) if os.path.splitext(BamDir)[1] else BamDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if qcResDir:
        _tgt = os.path.dirname(qcResDir) if os.path.splitext(qcResDir)[1] else qcResDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()