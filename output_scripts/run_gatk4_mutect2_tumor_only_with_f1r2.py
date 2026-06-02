#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = gatk4_mutect2
# VERSION = 4.4.0.0
# THREADS = 14
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
    parser.add_argument('--raw_vcf', required=False, default='[vcfDir]/[SeqID].mutect2.vcf', help='Filtered 전의 원본 VCF (Default: [vcfDir]/[SeqID].mutect2.vcf)')
    parser.add_argument('--f1r2_tar_gz', required=False, default='[qcResDir]/[SeqID].f1r2.tar.gz', help='Read orientation bias 통계 파일 (Default: [qcResDir]/[SeqID].f1r2.tar.gz)')
    parser.add_argument('--mutect_stats', required=False, default='[vcfDir]/[SeqID].mutect2.vcf.stats', help='필터링에 필요한 통계 메타데이터 (Default: [vcfDir]/[SeqID].mutect2.vcf.stats)')
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
    raw_vcf = args.raw_vcf
    f1r2_tar_gz = args.f1r2_tar_gz
    mutect_stats = args.mutect_stats
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
    if not raw_vcf:
    raw_vcf = f"{vcfDir}/{SeqID}.mutect2.vcf"
    if not f1r2_tar_gz:
    f1r2_tar_gz = f"{qcResDir}/{SeqID}.f1r2.tar.gz"
    if not mutect_stats:
    mutect_stats = f"{vcfDir}/{SeqID}.mutect2.vcf.stats"

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind} {gatk4_sif} gatk Mutect2 --java-options '-XX:+UseParallelGC -XX:ParallelGCThreads={Threads} -Xmx{xmx_mb}m' -pairHMM {pairHMM} --reference {ReferenceFasta} --intervals {TargetInterval} --input {BamDir}/{SeqID}.{InputSuffix}.bam --germline-resource {VcfGnomad} --panel-of-normals {VcfPon} --f1r2-tar-gz {f1r2_tar_gz} --f1r2-max-depth {f1r2_max_depth} --min-base-quality-score {min_base_q} --base-quality-score-threshold {min_base_q} --output {raw_vcf} {extra_args}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    if mutect_stats:
        _tgt = os.path.dirname(mutect_stats) if os.path.splitext(mutect_stats)[1] else mutect_stats
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if qcResDir:
        _tgt = os.path.dirname(qcResDir) if os.path.splitext(qcResDir)[1] else qcResDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if f1r2_tar_gz:
        _tgt = os.path.dirname(f1r2_tar_gz) if os.path.splitext(f1r2_tar_gz)[1] else f1r2_tar_gz
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if vcfDir:
        _tgt = os.path.dirname(vcfDir) if os.path.splitext(vcfDir)[1] else vcfDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if raw_vcf:
        _tgt = os.path.dirname(raw_vcf) if os.path.splitext(raw_vcf)[1] else raw_vcf
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if BamDir:
        _tgt = os.path.dirname(BamDir) if os.path.splitext(BamDir)[1] else BamDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()