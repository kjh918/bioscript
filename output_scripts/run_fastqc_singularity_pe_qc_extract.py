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
    parser.add_argument('--r1_html', required=False, default='[qcResDir]/[SeqID]_R1_fastqc.html', help='Quality report for Read 1 in HTML (Default: [qcResDir]/[SeqID]_R1_fastqc.html)')
    parser.add_argument('--r2_html', required=False, default='[qcResDir]/[SeqID]_R2_fastqc.html', help='Quality report for Read 2 in HTML (Default: [qcResDir]/[SeqID]_R2_fastqc.html)')
    parser.add_argument('--r1_zip', required=False, default='[qcResDir]/[SeqID]_R1_fastqc.zip', help='Data zip file for Read 1 (Default: [qcResDir]/[SeqID]_R1_fastqc.zip)')
    parser.add_argument('--r2_zip', required=False, default='[qcResDir]/[SeqID]_R2_fastqc.zip', help='Data zip file for Read 2 (Default: [qcResDir]/[SeqID]_R2_fastqc.zip)')
    parser.add_argument('--r1_dir', required=False, default='[qcResDir]/[SeqID]_R1_fastqc', help='Extracted report directory for Read 1 (Default: [qcResDir]/[SeqID]_R1_fastqc)')
    parser.add_argument('--r2_dir', required=False, default='[qcResDir]/[SeqID]_R2_fastqc', help='Extracted report directory for Read 2 (Default: [qcResDir]/[SeqID]_R2_fastqc)')
    parser.add_argument('--singularity_bin', required=False, default='singularity', help='Path to singularity executable (Default: singularity)')
    parser.add_argument('--fastqc_bin', required=False, default='fastqc', help='FastQC binary name inside container (Default: fastqc)')
    parser.add_argument('--sif', required=False, default='/storage/images/fastqc-0.12.1.sif', help='Path to FastQC SIF image (Default: /storage/images/fastqc-0.12.1.sif)')
    parser.add_argument('--bind', required=False, default='/storage,/data', help='Mount paths for singularity exec (Default: /storage,/data)')
    parser.add_argument('--Threads', required=False, default='8', help='Number of threads (FastQC processes files in parallel) (Default: 8)')
    parser.add_argument('--fastqc_args', required=False, default='--extract', help='Additional FastQC flags (e.g., --extract, --nogroup) (Default: --extract)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    RawFastqDir = args.RawFastqDir
    qcResDir = args.qcResDir
    r1_html = args.r1_html
    r2_html = args.r2_html
    r1_zip = args.r1_zip
    r2_zip = args.r2_zip
    r1_dir = args.r1_dir
    r2_dir = args.r2_dir
    singularity_bin = args.singularity_bin
    fastqc_bin = args.fastqc_bin
    sif = args.sif
    bind = args.bind
    Threads = args.Threads
    fastqc_args = args.fastqc_args

    # --- [Output Paths] ---
    if not r1_html:
    r1_html = f"{qcResDir}/{SeqID}_R1_fastqc.html"
    if not r2_html:
    r2_html = f"{qcResDir}/{SeqID}_R2_fastqc.html"
    if not r1_zip:
    r1_zip = f"{qcResDir}/{SeqID}_R1_fastqc.zip"
    if not r2_zip:
    r2_zip = f"{qcResDir}/{SeqID}_R2_fastqc.zip"
    if not r1_dir:
    r1_dir = f"{qcResDir}/{SeqID}_R1_fastqc"
    if not r2_dir:
    r2_dir = f"{qcResDir}/{SeqID}_R2_fastqc"

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind} {sif} {fastqc_bin} {fastqc_args} --threads {Threads} --outdir {qcResDir} {RawFastqDir}/{SeqID}_R1.fastq.gz {RawFastqDir}/{SeqID}_R2.fastq.gz"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    if r2_zip:
        _tgt = os.path.dirname(r2_zip) if os.path.splitext(r2_zip)[1] else r2_zip
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if r2_html:
        _tgt = os.path.dirname(r2_html) if os.path.splitext(r2_html)[1] else r2_html
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if r1_zip:
        _tgt = os.path.dirname(r1_zip) if os.path.splitext(r1_zip)[1] else r1_zip
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if r2_dir:
        _tgt = os.path.dirname(r2_dir) if os.path.splitext(r2_dir)[1] else r2_dir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if r1_dir:
        _tgt = os.path.dirname(r1_dir) if os.path.splitext(r1_dir)[1] else r1_dir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if RawFastqDir:
        _tgt = os.path.dirname(RawFastqDir) if os.path.splitext(RawFastqDir)[1] else RawFastqDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if qcResDir:
        _tgt = os.path.dirname(qcResDir) if os.path.splitext(qcResDir)[1] else qcResDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if r1_html:
        _tgt = os.path.dirname(r1_html) if os.path.splitext(r1_html)[1] else r1_html
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()