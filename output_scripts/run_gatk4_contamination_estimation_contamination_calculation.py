#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = gatk4_contamination_estimation
# VERSION = 4.4.0.0
# THREADS = 14
# PROFILE = contamination_calculation

"""
Tool: gatk4_contamination_estimation (4.4.0.0)
Profile: contamination_calculation
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="gatk4_contamination_estimation Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='No description (Default: )')
    parser.add_argument('--BamDir', required=True, default='', help='No description (Default: )')
    parser.add_argument('--qcResDir', required=True, default='', help='결과 테이블(.table) 저장 경로 (Default: )')
    parser.add_argument('--ReferenceFasta', required=True, default='', help='No description (Default: )')
    parser.add_argument('--TargetInterval', required=True, default='', help='분석 대상 영역 (Interval_list) (Default: )')
    parser.add_argument('--VcfGnomad', required=True, default='', help='gnomAD germline resource VCF (Required for pileup) (Default: )')
    parser.add_argument('--pileup_table', required=False, default='[qcResDir]/[SeqID].targeted_sequencing.table', help='Pileup summary statistics (Default: [qcResDir]/[SeqID].targeted_sequencing.table)')
    parser.add_argument('--contamination_table', required=False, default='[qcResDir]/[SeqID].contamination.table', help='Final contamination estimation (Default: [qcResDir]/[SeqID].contamination.table)')
    parser.add_argument('--InputSuffix', required=False, default='analysisReady', help='No description (Default: analysisReady)')
    parser.add_argument('--singularity_bin', required=False, default='singularity', help='No description (Default: singularity)')
    parser.add_argument('--gatk4_sif', required=False, default='/storage/images/gatk-4.4.0.0.sif', help='No description (Default: /storage/images/gatk-4.4.0.0.sif)')
    parser.add_argument('--bind', required=False, default='/storage,/data', help='No description (Default: /storage,/data)')
    parser.add_argument('--Threads', required=False, default='14', help='No description (Default: 14)')
    parser.add_argument('--pileup_xmx', required=False, default='8192', help='Xmx for GetPileupSummaries (8G) (Default: 8192)')
    parser.add_argument('--contam_xmx', required=False, default='16384', help='Xmx for CalculateContamination (16G) (Default: 16384)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    BamDir = args.BamDir
    qcResDir = args.qcResDir
    ReferenceFasta = args.ReferenceFasta
    TargetInterval = args.TargetInterval
    VcfGnomad = args.VcfGnomad
    pileup_table = args.pileup_table
    contamination_table = args.contamination_table
    InputSuffix = args.InputSuffix
    singularity_bin = args.singularity_bin
    gatk4_sif = args.gatk4_sif
    bind = args.bind
    Threads = args.Threads
    pileup_xmx = args.pileup_xmx
    contam_xmx = args.contam_xmx

    # --- [Output Paths] ---
    if not pileup_table:
    pileup_table = f"{qcResDir}/{SeqID}.targeted_sequencing.table"
    if not contamination_table:
    contamination_table = f"{qcResDir}/{SeqID}.contamination.table"

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind} {gatk4_sif} gatk GetPileupSummaries --java-options '-XX:ParallelGCThreads={Threads} -Xmx{pileup_xmx}m' -I {BamDir}/{SeqID}.{InputSuffix}.bam -V {VcfGnomad} -L {TargetInterval} -R {ReferenceFasta} -O {pileup_table} && {singularity_bin} exec -B {bind} {gatk4_sif} gatk CalculateContamination --java-options '-XX:ParallelGCThreads={Threads} -Xmx{contam_xmx}m' -I {pileup_table} -O {contamination_table}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    if BamDir:
        _tgt = os.path.dirname(BamDir) if os.path.splitext(BamDir)[1] else BamDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if pileup_table:
        _tgt = os.path.dirname(pileup_table) if os.path.splitext(pileup_table)[1] else pileup_table
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if contamination_table:
        _tgt = os.path.dirname(contamination_table) if os.path.splitext(contamination_table)[1] else contamination_table
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if qcResDir:
        _tgt = os.path.dirname(qcResDir) if os.path.splitext(qcResDir)[1] else qcResDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()