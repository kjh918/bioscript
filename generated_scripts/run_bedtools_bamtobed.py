#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = bedtools
# VERSION = 2.27.1
# THREADS = 1
# PROFILE = bamtobed

"""
Tool: bedtools (2.27.1)
Profile: bamtobed
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="bedtools Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='Sequence identifier used for file naming (Default: )')
    parser.add_argument('--BamDir', required=True, default='', help='Directory containing the input BAM file (Default: )')
    parser.add_argument('--BamToBedDir', required=True, default='', help='Target directory for final delivery of BED files (Default: )')
    parser.add_argument('--InputSuffix', required=False, default='analysisReady', help='Suffix of the input BAM (e.g., analysisReady, recal, sorted) (Default: analysisReady)')
    parser.add_argument('--bedtools_bin', required=False, default='bedtools', help='Path to bedtools binary (Default: bedtools)')
    parser.add_argument('--bgzip_bin', required=False, default='bgzip', help='Path to bgzip binary (Default: bgzip)')
    parser.add_argument('--ln_bin', required=False, default='ln', help='Path to ln binary (Default: ln)')
    parser.add_argument('--cp_bin', required=False, default='cp', help='Path to cp binary (Default: cp)')
    parser.add_argument('--Threads', required=False, default='8', help='Number of threads for bgzip compression (Default: 8)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    BamDir = args.BamDir
    BamToBedDir = args.BamToBedDir
    InputSuffix = args.InputSuffix
    bedtools_bin = args.bedtools_bin
    bgzip_bin = args.bgzip_bin
    ln_bin = args.ln_bin
    cp_bin = args.cp_bin
    Threads = args.Threads

    # --- [Output Paths] ---
    bed_gz = f"{BamDir}/{SeqID}.{InputSuffix}.bed.gz"
    final_bed = f"{BamToBedDir}/{SeqID}.bed.gz"

    # --- [Command Execution] ---
    cmd = f"{bedtools_bin} bamtobed -i {BamDir}/{SeqID}.{InputSuffix}.bam > {BamDir}/{SeqID}.{InputSuffix}.bed && {bgzip_bin} -f -@ {Threads} {BamDir}/{SeqID}.{InputSuffix}.bed && {ln_bin} -Tsf {bed_gz} {BamDir}/{SeqID}.bed.gz && {cp_bin} {bed_gz} {final_bed}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    os.makedirs(os.path.dirname(BamDir) if '.' in os.path.basename(BamDir) else BamDir, exist_ok=True)
    os.makedirs(os.path.dirname(BamToBedDir) if '.' in os.path.basename(BamToBedDir) else BamToBedDir, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()