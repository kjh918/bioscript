#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = gatk4
# VERSION = 4.4.0.0
# THREADS = 1
# PROFILE = bqsr_using_singularity

"""
Tool: gatk4 (4.4.0.0)
Profile: bqsr_using_singularity
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="gatk4 Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='Sequence identifier used for file naming (Default: )')
    parser.add_argument('--BamDir', required=True, default='', help='Directory containing the input BAM and target output BAM (Default: )')
    parser.add_argument('--qcResDir', required=True, default='', help='Directory for the recalibration table output (Default: )')
    parser.add_argument('--ReferenceFasta', required=True, default='', help='Path to the reference genome FASTA file (Default: )')
    parser.add_argument('--KnownSnp', required=True, default='', help='Path to known SNPs VCF (e.g., dbSNP) (Default: )')
    parser.add_argument('--KnownIndel1', required=True, default='', help='Path to known Indels VCF (e.g., Mills and 1000G) (Default: )')
    parser.add_argument('--KnownIndel2', required=True, default='', help='Path to additional known Indels VCF (Default: )')
    parser.add_argument('--InputSuffix', required=False, default='realign', help='Suffix of input BAM (e.g., dedup, sorted, primary) (Default: dedup)')
    parser.add_argument('--OutputSuffix', required=False, default='recal', help='Suffix for output BAM (e.g., recal, bqsr) (Default: recal)')
    parser.add_argument('--singularity_bin', required=False, default='singularity', help='Path to singularity executable (Default: singularity)')
    parser.add_argument('--gatk_bin', required=False, default='gatk', help='GATK wrapper script or binary inside container (Default: gatk)')
    parser.add_argument('--sif', required=False, default='/storage/images/gatk-4.4.0.0.sif', help='GATK singularity image file (Default: /storage/images/gatk-4.4.0.0.sif)')
    parser.add_argument('--bind', required=False, default='/storage,/data', help='Mount paths for singularity (Default: /storage,/data)')
    parser.add_argument('--Threads', required=False, default='14', help='Threads for ParallelGC and GATK engine (Default: 14)')
    parser.add_argument('--xmx_mb', required=False, default='16384', help='Maximum Java heap memory (MB) (Default: 16384)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    BamDir = args.BamDir
    qcResDir = args.qcResDir
    ReferenceFasta = args.ReferenceFasta
    KnownSnp = args.KnownSnp
    KnownIndel1 = args.KnownIndel1
    KnownIndel2 = args.KnownIndel2
    InputSuffix = args.InputSuffix
    OutputSuffix = args.OutputSuffix
    singularity_bin = args.singularity_bin
    gatk_bin = args.gatk_bin
    sif = args.sif
    bind = args.bind
    Threads = args.Threads
    xmx_mb = args.xmx_mb

    # --- [Output Paths] ---
    recal_table = f"{qcResDir}/{SeqID}.{InputSuffix}.recal_table.txt"
    recal_bam = f"{BamDir}/{SeqID}.{OutputSuffix}.bam"
    recal_bai = f"{BamDir}/{SeqID}.{OutputSuffix}.bam.bai"

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind} {sif} {gatk_bin} BaseRecalibrator --java-options '-XX:ParallelGCThreads={Threads} -Xmx{xmx_mb}m' --input {BamDir}/{SeqID}.{InputSuffix}.bam --reference {ReferenceFasta} --output {recal_table} --known-sites {KnownSnp} --known-sites {KnownIndel1} --known-sites {KnownIndel2} && {singularity_bin} exec -B {bind} {sif} {gatk_bin} ApplyBQSR --java-options '-XX:ParallelGCThreads={Threads} -Xmx{xmx_mb}m' --input {BamDir}/{SeqID}.{InputSuffix}.bam --bqsr-recal-file {recal_table} --output {recal_bam} --create-output-bam-index true"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    os.makedirs(os.path.dirname(BamDir) if '.' in os.path.basename(BamDir) else BamDir, exist_ok=True)
    os.makedirs(os.path.dirname(qcResDir) if '.' in os.path.basename(qcResDir) else qcResDir, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()