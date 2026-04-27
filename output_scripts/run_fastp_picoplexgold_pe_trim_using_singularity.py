#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = fastp
# VERSION = 0.23.4
# THREADS = 1
# PROFILE = PicoplexGold_pe_trim_using_singularity

"""
Tool: fastp (0.23.4)
Profile: PicoplexGold_pe_trim_using_singularity
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="fastp Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='Sequence identifier used for file naming and sample tracking (Default: )')
    parser.add_argument('--RawFastqDir', required=True, default='', help='Directory containing the raw input FASTQ files (Default: )')
    parser.add_argument('--TrimFastqDir', required=True, default='', help='Target directory for the processed (trimmed) FASTQ files (Default: )')
    parser.add_argument('--qcResDir', required=True, default='', help='Directory for quality control reports (JSON/HTML) (Default: )')
    parser.add_argument('--singularity_bin', required=False, default='singularity', help='Path to the singularity executable (Default: singularity)')
    parser.add_argument('--fastp_bin', required=False, default='fastp', help='fastp executable name or path inside container (Default: fastp)')
    parser.add_argument('--sif', required=False, default='/storage/images/fastp-0.23.4.sif', help='Path to fastp Singularity image (.sif) (Default: /storage/images/fastp-0.23.4.sif)')
    parser.add_argument('--bind', required=False, default='/storage,/data', help='Directories to mount for Singularity (Default: /storage,/data)')
    parser.add_argument('--Threads', required=False, default='8', help='Number of threads for parallel processing (Default: 8)')
    parser.add_argument('--length_required', required=False, default='100', help='Minimum read length filter (Default: 100)')
    parser.add_argument('--average_qual', required=False, default='10', help='Minimum average quality filter (Default: 10)')
    parser.add_argument('--qualified_quality_phred', required=False, default='15', help='Phred quality score threshold (base quality) (Default: 15)')
    parser.add_argument('--trim_front1', required=False, default='14', help='Number of bases to trim from the start of Read 1 (Default: 14)')
    parser.add_argument('--trim_front2', required=False, default='14', help='Number of bases to trim from the start of Read 2 (Default: 14)')
    parser.add_argument('--adapter_sequence', required=False, default='AGATCGGAAGAGCACACGTCTGAACTCCAGTCA', help='Adapter sequence for Read 1 (Default: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA)')
    parser.add_argument('--adapter_sequence_r2', required=False, default='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT', help='Adapter sequence for Read 2 (Default: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT)')
    
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
    trim_front1 = args.trim_front1
    trim_front2 = args.trim_front2
    adapter_sequence = args.adapter_sequence
    adapter_sequence_r2 = args.adapter_sequence_r2

    # --- [Output Paths] ---
    out_read1 = f"{TrimFastqDir}/{SeqID}.trimmed_R1.fastq.gz"
    out_read2 = f"{TrimFastqDir}/{SeqID}.trimmed_R2.fastq.gz"
    json = f"{qcResDir}/{SeqID}.fastp.json"
    html = f"{qcResDir}/{SeqID}.fastp.html"

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind} {sif} {fastp_bin} --thread {Threads} --in1 {RawFastqDir}/{SeqID}_R1.fastq.gz --in2 {RawFastqDir}/{SeqID}_R2.fastq.gz --out1 {out_read1} --out2 {out_read2} --json {json} --html {html} --trim_poly_g --detect_adapter_for_pe --adapter_sequence {adapter_sequence} --adapter_sequence_r2 {adapter_sequence_r2} --length_required {length_required} --average_qual {average_qual} --qualified_quality_phred {qualified_quality_phred} --trim_front1 {trim_front1} --trim_front2 {trim_front2}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    os.makedirs(os.path.dirname(RawFastqDir) if '.' in os.path.basename(RawFastqDir) else RawFastqDir, exist_ok=True)
    os.makedirs(os.path.dirname(TrimFastqDir) if '.' in os.path.basename(TrimFastqDir) else TrimFastqDir, exist_ok=True)
    os.makedirs(os.path.dirname(qcResDir) if '.' in os.path.basename(qcResDir) else qcResDir, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()