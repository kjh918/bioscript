#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = gatk4_filter_mutect_calls
# VERSION = 4.4.0.0
# THREADS = 1
# PROFILE = somatic_variant_filtering

"""
Tool: gatk4_filter_mutect_calls (4.4.0.0)
Profile: somatic_variant_filtering
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="gatk4_filter_mutect_calls Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='No description (Default: )')
    parser.add_argument('--vcfDir', required=True, default='', help='입력 및 출력 VCF가 위치한 경로 (Default: )')
    parser.add_argument('--qcResDir', required=True, default='', help='Contamination 및 OB-prior 파일 경로 (Default: )')
    parser.add_argument('--ReferenceFasta', required=True, default='', help='No description (Default: )')
    parser.add_argument('--TargetInterval', required=True, default='', help='No description (Default: )')
    parser.add_argument('--ContaminationTable', required=True, default='', help='[SeqID].contamination.table (Default: )')
    parser.add_argument('--Suffix', required=False, default='', help='File name suffix (e.g., \'keep.germline.\') (Default: )')
    parser.add_argument('--singularity_bin', required=False, default='singularity', help='No description (Default: singularity)')
    parser.add_argument('--gatk4_sif', required=False, default='/storage/images/gatk-4.4.0.0.sif', help='No description (Default: /storage/images/gatk-4.4.0.0.sif)')
    parser.add_argument('--bind', required=False, default='/storage,/data', help='No description (Default: /storage,/data)')
    parser.add_argument('--xmx_mb', required=False, default='16384', help='16GB memory from original script (Default: 16384)')
    parser.add_argument('--Threads', required=False, default='14', help='No description (Default: 14)')
    parser.add_argument('--extra_args', required=False, default='', help='No description (Default: )')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    vcfDir = args.vcfDir
    qcResDir = args.qcResDir
    ReferenceFasta = args.ReferenceFasta
    TargetInterval = args.TargetInterval
    ContaminationTable = args.ContaminationTable
    Suffix = args.Suffix
    singularity_bin = args.singularity_bin
    gatk4_sif = args.gatk4_sif
    bind = args.bind
    xmx_mb = args.xmx_mb
    Threads = args.Threads
    extra_args = args.extra_args

    # --- [Output Paths] ---
    filtered_vcf = f"{vcfDir}/{SeqID}.mutect2.{Suffix}filtered.vcf"

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind} {gatk4_sif} gatk FilterMutectCalls --java-options '-XX:ParallelGCThreads={Threads} -Xmx{xmx_mb}m' -V {vcfDir}/{SeqID}.mutect2.{Suffix}vcf -L {TargetInterval} --reference {ReferenceFasta} --contamination-table {ContaminationTable} --ob-priors {qcResDir}/{SeqID}.{Suffix}read-orientation-model.tar.gz -O {filtered_vcf} {extra_args}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    os.makedirs(os.path.dirname(vcfDir) if '.' in os.path.basename(vcfDir) else vcfDir, exist_ok=True)
    os.makedirs(os.path.dirname(qcResDir) if '.' in os.path.basename(qcResDir) else qcResDir, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()