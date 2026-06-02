#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = bedtools
# VERSION = 2.27.1
# THREADS = 8
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
    parser.add_argument('--bed_gz', required=False, default='[BamDir]/[SeqID].[InputSuffix].bed.gz', help='Compressed BED file in the working directory (Default: [BamDir]/[SeqID].[InputSuffix].bed.gz)')
    parser.add_argument('--final_bed', required=False, default='[BamToBedDir]/[SeqID].bed.gz', help='Final BED file delivered to the target directory (Default: [BamToBedDir]/[SeqID].bed.gz)')
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
    bed_gz = args.bed_gz
    final_bed = args.final_bed
    InputSuffix = args.InputSuffix
    bedtools_bin = args.bedtools_bin
    bgzip_bin = args.bgzip_bin
    ln_bin = args.ln_bin
    cp_bin = args.cp_bin
    Threads = args.Threads

    # --- [Output Paths] ---
    if not bed_gz:
    bed_gz = f"{BamDir}/{SeqID}.{InputSuffix}.bed.gz"
    if not final_bed:
    final_bed = f"{BamToBedDir}/{SeqID}.bed.gz"

    # --- [Command Execution] ---
    cmd = f"{bedtools_bin} bamtobed -i {BamDir}/{SeqID}.{InputSuffix}.bam > {BamDir}/{SeqID}.{InputSuffix}.bed && {bgzip_bin} -f -@ {Threads} {BamDir}/{SeqID}.{InputSuffix}.bed && {ln_bin} -Tsf {bed_gz} {BamDir}/{SeqID}.bed.gz && {cp_bin} {bed_gz} {final_bed}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    if bed_gz:
        _tgt = os.path.dirname(bed_gz) if os.path.splitext(bed_gz)[1] else bed_gz
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if final_bed:
        _tgt = os.path.dirname(final_bed) if os.path.splitext(final_bed)[1] else final_bed
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if BamToBedDir:
        _tgt = os.path.dirname(BamToBedDir) if os.path.splitext(BamToBedDir)[1] else BamToBedDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if BamDir:
        _tgt = os.path.dirname(BamDir) if os.path.splitext(BamDir)[1] else BamDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()