#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = sob_detector
# VERSION = 1.0.4
# THREADS = 1
# PROFILE = somatic_artifact_cleanup

"""
Tool: sob_detector (1.0.4)
Profile: somatic_artifact_cleanup
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="sob_detector Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='No description (Default: )')
    parser.add_argument('--vcfDir', required=True, default='', help='필터링할 VCF가 위치한 경로 (Default: )')
    parser.add_argument('--BamDir', required=True, default='', help='분석용 BAM(analysisReady) 경로 (Default: )')
    parser.add_argument('--InputVcfSuffix', required=False, default='filtered', help='No description (Default: filtered)')
    parser.add_argument('--InputBamSuffix', required=False, default='analysisReady', help='No description (Default: analysisReady)')
    parser.add_argument('--singularity_bin', required=False, default='singularity', help='No description (Default: singularity)')
    parser.add_argument('--sob_sif', required=False, default='/storage/images/sobdetector.sif', help='No description (Default: /storage/images/sobdetector.sif)')
    parser.add_argument('--sob_bin', required=False, default='SOBDetector', help='No description (Default: SOBDetector)')
    parser.add_argument('--bind', required=False, default='/storage,/data', help='No description (Default: /storage,/data)')
    parser.add_argument('--minBaseQ', required=False, default='20', help='No description (Default: 20)')
    parser.add_argument('--minMapQ', required=False, default='20', help='No description (Default: 20)')
    parser.add_argument('--only_passed', required=False, default='false', help='PASS 필터된 변이만 처리할지 여부 (Default: false)')
    parser.add_argument('--extra_args', required=False, default='', help='No description (Default: )')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    vcfDir = args.vcfDir
    BamDir = args.BamDir
    InputVcfSuffix = args.InputVcfSuffix
    InputBamSuffix = args.InputBamSuffix
    singularity_bin = args.singularity_bin
    sob_sif = args.sob_sif
    sob_bin = args.sob_bin
    bind = args.bind
    minBaseQ = args.minBaseQ
    minMapQ = args.minMapQ
    only_passed = args.only_passed
    extra_args = args.extra_args

    # --- [Output Paths] ---
    bias_filtered_vcf = f"{vcfDir}/{SeqID}.mutect2.bias.filtered.vcf"

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind} {sob_sif} {sob_bin}  --input-type VCF --input-variants {vcfDir}/{SeqID}.mutect2.{InputVcfSuffix}.vcf --input-bam {BamDir}/{SeqID}.{InputBamSuffix}.bam --output-variants {bias_filtered_vcf} --minBaseQuality {minBaseQ} --minMappingQuality {minMapQ} --only-passed {only_passed} {extra_args}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    os.makedirs(os.path.dirname(vcfDir) if '.' in os.path.basename(vcfDir) else vcfDir, exist_ok=True)
    os.makedirs(os.path.dirname(BamDir) if '.' in os.path.basename(BamDir) else BamDir, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()