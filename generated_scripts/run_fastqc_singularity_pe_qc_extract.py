#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = fastqc
# VERSION = 0.12.1
# THREADS = 8
# PROFILE = singularity_pe_qc_extract

"""
Tool: fastqc (0.12.1)
Profile: singularity_pe_qc_extract
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="fastqc Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='Sequence identifier used to locate FASTQ files (Default: )')
    parser.add_argument('--RawFastqDir', required=True, default='', help='Directory containing raw or trimmed FASTQ files (Default: )')
    parser.add_argument('--qcResDir', required=True, default='', help='Directory where FastQC reports will be saved (Default: )')
    parser.add_argument('--singularity_bin', required=False, default='singularity', help='Path to singularity executable (Default: singularity)')
    parser.add_argument('--fastqc_bin', required=False, default='fastqc', help='FastQC binary name inside container (Default: fastqc)')
    parser.add_argument('--sif', required=False, default='/storage/images/fastqc-0.12.1.sif', help='Path to FastQC SIF image (Default: /storage/images/fastqc-0.12.1.sif)')
    parser.add_argument('--bind', required=False, default='/storage,/data', help='Mount paths for singularity exec (Default: /storage,/data)')
    parser.add_argument('--threads', required=False, default='8', help='Number of threads (FastQC processes files in parallel) (Default: 8)')
    parser.add_argument('--fastqc_args', required=False, default='--extract', help='Additional FastQC flags (e.g., --extract, --nogroup) (Default: --extract)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    RawFastqDir = args.RawFastqDir
    qcResDir = args.qcResDir
    singularity_bin = args.singularity_bin
    fastqc_bin = args.fastqc_bin
    sif = args.sif
    bind = args.bind
    threads = args.threads
    fastqc_args = args.fastqc_args

    # --- [Output Paths] ---
    r1_html = f"{qcResDir}/{SeqID}_R1_fastqc.html"
    r2_html = f"{qcResDir}/{SeqID}_R2_fastqc.html"
    r1_zip = f"{qcResDir}/{SeqID}_R1_fastqc.zip"
    r2_zip = f"{qcResDir}/{SeqID}_R2_fastqc.zip"
    r1_dir = f"{qcResDir}/{SeqID}_R1_fastqc"
    r2_dir = f"{qcResDir}/{SeqID}_R2_fastqc"

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind} {sif} {fastqc_bin} {fastqc_args} --threads {threads} --outdir {qcResDir} {RawFastqDir}/{SeqID}_R1.fastq.gz {RawFastqDir}/{SeqID}_R2.fastq.gz"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    os.makedirs(os.path.dirname(RawFastqDir) if '.' in os.path.basename(RawFastqDir) else RawFastqDir, exist_ok=True)
    os.makedirs(os.path.dirname(qcResDir) if '.' in os.path.basename(qcResDir) else qcResDir, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()