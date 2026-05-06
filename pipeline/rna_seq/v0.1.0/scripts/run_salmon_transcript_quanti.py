#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = salmon
# VERSION = 1.10.0
# THREADS = 1
# PROFILE = transcript_quanti

"""
Tool: salmon (1.10.0)
Profile: transcript_quanti
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="salmon Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='No description (Default: )')
    parser.add_argument('--BamDir', required=True, default='', help='Directory containing input BAM files (Default: )')
    parser.add_argument('--OutputDir', required=True, default='', help='Path to output directory (Default: )')
    parser.add_argument('--SalmonIndex', required=True, default='', help='Path to Salmon index directory (Default: )')
    parser.add_argument('--Threads', required=False, default='12', help='No description (Default: 12)')
    parser.add_argument('--libType', required=False, default='A', help='Library type (A: Auto-detect) (Default: A)')
    parser.add_argument('--InputSuffix', required=False, default='.Aligned.toTranscriptome.out.bam', help='No description (Default: .Aligned.toTranscriptome.out.bam)')
    parser.add_argument('--extra_args', required=False, default='', help='No description (Default: )')
    parser.add_argument('--samlon_bin', required=False, default='/storage/apps/salmon-1.10.0/bin/salmon', help='No description (Default: /storage/apps/salmon-1.10.0/bin/salmon)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    BamDir = args.BamDir
    OutputDir = args.OutputDir
    SalmonIndex = args.SalmonIndex
    Threads = args.Threads
    libType = args.libType
    InputSuffix = args.InputSuffix
    extra_args = args.extra_args
    samlon_bin = args.samlon_bin

    # --- [Output Paths] ---
    quant_sf = f"{OutputDir}/{SeqID}_quant/quant.sf"
    quant_log = f"{OutputDir}/{SeqID}_quant/logs/salmon_quant.log"
    lib_format = f"{OutputDir}/{SeqID}_quant/lib_format_counts.json"

    # --- [Command Execution] ---
    cmd = f"{samlon_bin} quant -t {SalmonIndex} --threads {Threads} --libType {libType} -a {BamDir}/{SeqID}{InputSuffix} -o {OutputDir} {extra_args}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    os.makedirs(os.path.dirname(BamDir) if '.' in os.path.basename(BamDir) else BamDir, exist_ok=True)
    os.makedirs(os.path.dirname(OutputDir) if '.' in os.path.basename(OutputDir) else OutputDir, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()