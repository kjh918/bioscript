#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = gatk3
# VERSION = 3.8
# THREADS = 1
# PROFILE = indel_realign_using_singularity

"""
Tool: gatk3 (3.8)
Profile: indel_realign_using_singularity
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="gatk3 Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='Sequence identifier used for file naming (Default: )')
    parser.add_argument('--BamDir', required=True, default='', help='Directory containing the input BAM and target output BAM (Default: )')
    parser.add_argument('--qcResDir', required=True, default='', help='Directory for the target intervals file output (Default: )')
    parser.add_argument('--ReferenceFasta', required=True, default='', help='Path to the reference genome FASTA file (Default: )')
    parser.add_argument('--KnownIndel1', required=True, default='', help='Path to known Indels VCF (e.g., Mills and 1000G) (Default: )')
    parser.add_argument('--KnownIndel2', required=True, default='', help='Path to additional known Indels VCF (Default: )')
    parser.add_argument('--InputSuffix', required=False, default='dedup', help='Suffix of input BAM (e.g., dedup, sorted, primary) (Default: dedup)')
    parser.add_argument('--OutputSuffix', required=False, default='realign', help='Suffix for output BAM (e.g., realign, ir) (Default: realign)')
    parser.add_argument('--singularity_bin', required=False, default='singularity', help='Path to singularity executable (Default: singularity)')
    parser.add_argument('--java_bin', required=False, default='java', help='Java executable (inside container) (Default: java)')
    parser.add_argument('--gatk_jar', required=False, default='/usr/GenomeAnalysisTK.jar', help='GATK 3.8 jar path inside container (Default: /usr/GenomeAnalysisTK.jar)')
    parser.add_argument('--sif', required=False, default='/storage/images/gatk-3.8-1.sif', help='GATK 3.8 singularity image file (Default: /storage/images/gatk-3.8-1.sif)')
    parser.add_argument('--bind', required=False, default='/storage,/data', help='Mount paths for singularity (Default: /storage,/data)')
    parser.add_argument('--Threads', required=False, default='8', help='Number of threads for RealignerTargetCreator (-nt) (Default: 8)')
    parser.add_argument('--xmx_mb', required=False, default='16384', help='Maximum Java heap memory (MB) (Default: 16384)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    BamDir = args.BamDir
    qcResDir = args.qcResDir
    ReferenceFasta = args.ReferenceFasta
    KnownIndel1 = args.KnownIndel1
    KnownIndel2 = args.KnownIndel2
    InputSuffix = args.InputSuffix
    OutputSuffix = args.OutputSuffix
    singularity_bin = args.singularity_bin
    java_bin = args.java_bin
    gatk_jar = args.gatk_jar
    sif = args.sif
    bind = args.bind
    Threads = args.Threads
    xmx_mb = args.xmx_mb

    # --- [Output Paths] ---
    target_intervals = f"{qcResDir}/{SeqID}.{OutputSuffix}.realign_target.intervals"
    realigned_bam = f"{BamDir}/{SeqID}.{OutputSuffix}.bam"

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind} {sif} {java_bin} -Xmx{xmx_mb}m -jar {gatk_jar} -T RealignerTargetCreator -R {ReferenceFasta} -I {BamDir}/{SeqID}.{InputSuffix}.bam -o {target_intervals} -known {KnownIndel1} -known {KnownIndel2} -nt {Threads} && {singularity_bin} exec -B {bind} {sif} {java_bin} -Xmx{xmx_mb}m -jar {gatk_jar} -T IndelRealigner -R {ReferenceFasta} -targetIntervals {target_intervals} -known {KnownIndel1} -known {KnownIndel2} -I {BamDir}/{SeqID}.{InputSuffix}.bam -o {realigned_bam}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    os.makedirs(os.path.dirname(BamDir) if '.' in os.path.basename(BamDir) else BamDir, exist_ok=True)
    os.makedirs(os.path.dirname(qcResDir) if '.' in os.path.basename(qcResDir) else qcResDir, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()