#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = gatk4
# VERSION = 4.4.0.0
# THREADS = 14
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
    parser.add_argument('--KnownSnp1', required=True, default='', help='Path to known SNPs VCF (e.g., dbSNP) (Default: )')
    parser.add_argument('--KnownSnp2', required=True, default='', help='Path to known SNPs VCF (e.g., dbSNP) (Default: )')
    parser.add_argument('--KnownIndel1', required=True, default='', help='Path to known Indels VCF (e.g., Mills and 1000G) (Default: )')
    parser.add_argument('--KnownIndel2', required=True, default='', help='Path to additional known Indels VCF (Default: )')
    parser.add_argument('--recal_table', required=False, default='[qcResDir]/[SeqID].[InputSuffix].recal_table.txt', help='Recalibration report table used by ApplyBQSR (Default: [qcResDir]/[SeqID].[InputSuffix].recal_table.txt)')
    parser.add_argument('--recal_bam', required=False, default='[BamDir]/[SeqID].[OutputSuffix].bam', help='Final recalibrated BAM file (Default: [BamDir]/[SeqID].[OutputSuffix].bam)')
    parser.add_argument('--recal_bai', required=False, default='[BamDir]/[SeqID].[OutputSuffix].bam.bai', help='Index file for the recalibrated BAM (Default: [BamDir]/[SeqID].[OutputSuffix].bam.bai)')
    parser.add_argument('--InputSuffix', required=False, default='merged.dup.marked.realign', help='Suffix of input BAM (e.g., dedup, sorted, primary) (Default: merged.dup.marked.realign)')
    parser.add_argument('--OutputSuffix', required=False, default='merged.dup.marked.realign.recal', help='Suffix for output BAM (e.g., recal, bqsr) (Default: merged.dup.marked.realign.recal)')
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
    KnownSnp1 = args.KnownSnp1
    KnownSnp2 = args.KnownSnp2
    KnownIndel1 = args.KnownIndel1
    KnownIndel2 = args.KnownIndel2
    recal_table = args.recal_table
    recal_bam = args.recal_bam
    recal_bai = args.recal_bai
    InputSuffix = args.InputSuffix
    OutputSuffix = args.OutputSuffix
    singularity_bin = args.singularity_bin
    gatk_bin = args.gatk_bin
    sif = args.sif
    bind = args.bind
    Threads = args.Threads
    xmx_mb = args.xmx_mb

    # --- [Output Paths] ---
    if not recal_table:
    recal_table = f"{qcResDir}/{SeqID}.{InputSuffix}.recal_table.txt"
    if not recal_bam:
    recal_bam = f"{BamDir}/{SeqID}.{OutputSuffix}.bam"
    if not recal_bai:
    recal_bai = f"{BamDir}/{SeqID}.{OutputSuffix}.bam.bai"

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind} {sif} {gatk_bin} BaseRecalibrator --java-options '-XX:ParallelGCThreads={Threads} -Xmx{xmx_mb}m' --input {BamDir}/{SeqID}.{InputSuffix}.bam --reference {ReferenceFasta} --output {recal_table} --known-sites {KnownSnp1} --known-sites {KnownSnp2} --known-sites {KnownIndel1} --known-sites {KnownIndel2} && {singularity_bin} exec -B {bind} {sif} {gatk_bin} ApplyBQSR --java-options '-XX:ParallelGCThreads={Threads} -Xmx{xmx_mb}m' --input {BamDir}/{SeqID}.{InputSuffix}.bam --bqsr-recal-file {recal_table} --output {recal_bam} --create-output-bam-index true"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    if recal_bam:
        _tgt = os.path.dirname(recal_bam) if os.path.splitext(recal_bam)[1] else recal_bam
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if qcResDir:
        _tgt = os.path.dirname(qcResDir) if os.path.splitext(qcResDir)[1] else qcResDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if recal_table:
        _tgt = os.path.dirname(recal_table) if os.path.splitext(recal_table)[1] else recal_table
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if BamDir:
        _tgt = os.path.dirname(BamDir) if os.path.splitext(BamDir)[1] else BamDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if recal_bai:
        _tgt = os.path.dirname(recal_bai) if os.path.splitext(recal_bai)[1] else recal_bai
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()