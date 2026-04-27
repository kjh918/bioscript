#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = gatk4_mutect2
# VERSION = 4.4.0.0
# THREADS = 1
# PROFILE = tumor_only_with_f1r2

"""
Tool: gatk4_mutect2 (4.4.0.0)
Profile: tumor_only_with_f1r2
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="gatk4_mutect2 Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='No description (Default: )')
    parser.add_argument('--BamDir', required=True, default='', help='No description (Default: )')
    parser.add_argument('--vcfDir', required=True, default='', help='Raw VCF 결과 저장 경로 (Default: )')
    parser.add_argument('--qcResDir', required=True, default='', help='f1r2.tar.gz 및 중간 파일 저장 경로 (Default: )')
    parser.add_argument('--ReferenceFasta', required=True, default='', help='No description (Default: )')
    parser.add_argument('--TargetInterval', required=True, default='', help='Target BED 또는 interval_list (Default: )')
    parser.add_argument('--VcfGnomad', required=True, default='', help='gnomAD germline resource VCF (Default: )')
    parser.add_argument('--VcfPon', required=True, default='', help='Panel of Normals VCF (Default: )')
    parser.add_argument('--InputSuffix', required=False, default='analysisReady', help='No description (Default: analysisReady)')
    parser.add_argument('--singularity_bin', required=False, default='singularity', help='No description (Default: singularity)')
    parser.add_argument('--gatk4_sif', required=False, default='/storage/images/gatk-4.4.0.0.sif', help='No description (Default: /storage/images/gatk-4.4.0.0.sif)')
    parser.add_argument('--bind', required=False, default='/storage,/data', help='No description (Default: /storage,/data)')
    parser.add_argument('--xmx_mb', required=False, default='32768', help='Java Max Heap Size (MB) (Default: 32768)')
    parser.add_argument('--Threads', required=False, default='14', help='Parallel GC 및 OpenMP 스레드 수 (Default: 14)')
    parser.add_argument('--pairHMM', required=False, default='AVX_LOGLESS_CACHING_OMP', help='No description (Default: AVX_LOGLESS_CACHING_OMP)')
    parser.add_argument('--f1r2_max_depth', required=False, default='2500', help='No description (Default: 2500)')
    parser.add_argument('--min_base_q', required=False, default='20', help='No description (Default: 20)')
    parser.add_argument('--extra_args', required=False, default='', help='No description (Default: )')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    BamDir = args.BamDir
    vcfDir = args.vcfDir
    qcResDir = args.qcResDir
    ReferenceFasta = args.ReferenceFasta
    TargetInterval = args.TargetInterval
    VcfGnomad = args.VcfGnomad
    VcfPon = args.VcfPon
    InputSuffix = args.InputSuffix
    singularity_bin = args.singularity_bin
    gatk4_sif = args.gatk4_sif
    bind = args.bind
    xmx_mb = args.xmx_mb
    Threads = args.Threads
    pairHMM = args.pairHMM
    f1r2_max_depth = args.f1r2_max_depth
    min_base_q = args.min_base_q
    extra_args = args.extra_args

    # --- [Output Paths] ---
    raw_vcf = f"{vcfDir}/{SeqID}.mutect2.vcf"
    f1r2_tar_gz = f"{qcResDir}/{SeqID}.f1r2.tar.gz"
    mutect_stats = f"{vcfDir}/{SeqID}.mutect2.vcf.stats"

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind} {gatk4_sif} gatk Mutect2 --java-options '-XX:+UseParallelGC -XX:ParallelGCThreads={Threads} -Xmx{xmx_mb}m' -pairHMM {pairHMM} --reference {ReferenceFasta} --intervals {TargetInterval} --input {BamDir}/{SeqID}.{InputSuffix}.bam --germline-resource {VcfGnomad} --panel-of-normals {VcfPon} --f1r2-tar-gz {f1r2_tar_gz} --f1r2-max-depth {f1r2_max_depth} --min-base-quality-score {min_base_q} --base-quality-score-threshold {min_base_q} --output {raw_vcf} {extra_args}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    os.makedirs(os.path.dirname(BamDir) if '.' in os.path.basename(BamDir) else BamDir, exist_ok=True)
    os.makedirs(os.path.dirname(vcfDir) if '.' in os.path.basename(vcfDir) else vcfDir, exist_ok=True)
    os.makedirs(os.path.dirname(qcResDir) if '.' in os.path.basename(qcResDir) else qcResDir, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()