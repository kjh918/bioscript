#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = fastp
# VERSION = 0.23.4
# THREADS = 1
# PROFILE = pe_trim_using_singularity

"""
Tool: fastp (0.23.4)
Profile: pe_trim_using_singularity
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="fastp Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='Sequence identifier used for file naming (Default: )')
    parser.add_argument('--RawFastqDir', required=True, default='', help='Directory containing the raw input FASTQ files. Suffix = [SeqID]_R1.fastq.gz and [SeqID]_R2.fastq.gz (Default: )')
    parser.add_argument('--TrimFastqDir', required=True, default='', help='Directory where trimmed FASTQ files will be stored. Trimmed Read = Suffix = [SeqID].trimmed_R1.fastq.gz and [SeqID].trimmed_R2.fastq.gz (Default: )')
    parser.add_argument('--qcResDir', required=True, default='', help='Directory where fastp JSON and HTML reports will be stored (Default: )')
    parser.add_argument('--singularity_bin', required=False, default='singularity', help='Path to singularity executable (Default: singularity)')
    parser.add_argument('--fastp_bin', required=False, default='fastp', help='fastp binary name or path inside the container (Default: fastp)')
    parser.add_argument('--sif', required=False, default='/storage/images/fastp-0.23.4.sif', help='Path to fastp Singularity image file (Default: /storage/images/fastp-0.23.4.sif)')
    parser.add_argument('--bind', required=False, default='/storage,/data', help='Mount paths for Singularity execution (Default: /storage,/data)')
    parser.add_argument('--Threads', required=False, default='8', help='Number of threads to be used by fastp (Default: 8)')
    parser.add_argument('--length_required', required=False, default='100', help='Reads shorter than this length will be discarded (Default: 100)')
    parser.add_argument('--average_qual', required=False, default='10', help='Reads with average quality score lower than this will be discarded (Default: 10)')
    parser.add_argument('--qualified_quality_phred', required=False, default='15', help='The quality threshold that a base is considered as qualified (Default: 15)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    RawFastqDir = args.RawFastqDir
    TrimFastqDir = args.TrimFastqDir
    qcResDir = args.qcResDir
    singularity_bin = args.singularity_bin
    fastp_bin = args.fastp_bin
    sif = args.sif
    bind = args.bind
    Threads = args.Threads
    length_required = args.length_required
    average_qual = args.average_qual
    qualified_quality_phred = args.qualified_quality_phred

    # --- [Output Paths] ---
    out_read1 = f"{TrimFastqDir}/{SeqID}.trimmed_R1.fastq.gz"
    out_read2 = f"{TrimFastqDir}/{SeqID}.trimmed_R2.fastq.gz"
    json = f"{qcResDir}/{SeqID}.fastp.json"
    html = f"{qcResDir}/{SeqID}.fastp.html"

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind} {sif} {fastp_bin} --thread {Threads} --in1 {RawFastqDir}/{SeqID}_R1.fastq.gz --in2 {RawFastqDir}/{SeqID}_R2.fastq.gz --out1 {out_read1} --out2 {out_read2} --json {json} --html {html} --trim_poly_g --detect_adapter_for_pe --length_required {length_required} --average_qual {average_qual} --qualified_quality_phred {qualified_quality_phred}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    os.makedirs(os.path.dirname(RawFastqDir) if '.' in os.path.basename(RawFastqDir) else RawFastqDir, exist_ok=True)
    os.makedirs(os.path.dirname(TrimFastqDir) if '.' in os.path.basename(TrimFastqDir) else TrimFastqDir, exist_ok=True)
    os.makedirs(os.path.dirname(qcResDir) if '.' in os.path.basename(qcResDir) else qcResDir, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()