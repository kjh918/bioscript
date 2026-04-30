#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = salmon_quant
# VERSION = 1.10.1
# THREADS = 1
# PROFILE = transcript_quantification

"""
Tool: salmon_quant (1.10.1)
Profile: transcript_quantification
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="salmon_quant Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='No description (Default: )')
    parser.add_argument('--FastqDir', required=True, default='', help='No description (Default: )')
    parser.add_argument('--OutputDir', required=True, default='', help='정량화 결과(quant 폴더)가 저장될 경로 (Default: )')
    parser.add_argument('--SalmonIndex', required=True, default='', help='Path to Salmon index directory (Default: )')
    parser.add_argument('--Threads', required=False, default='12', help='No description (Default: 12)')
    parser.add_argument('--libType', required=False, default='A', help='Library type (A: Auto-detect) (Default: A)')
    parser.add_argument('--InputSuffix', required=False, default='fastq.gz', help='No description (Default: fastq.gz)')
    parser.add_argument('--validateMappings', required=False, default='--validateMappings', help='No description (Default: --validateMappings)')
    parser.add_argument('--extra_args', required=False, default='', help='No description (Default: )')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    FastqDir = args.FastqDir
    OutputDir = args.OutputDir
    SalmonIndex = args.SalmonIndex
    Threads = args.Threads
    libType = args.libType
    InputSuffix = args.InputSuffix
    validateMappings = args.validateMappings
    extra_args = args.extra_args

    # --- [Output Paths] ---
    quant_sf = f"{OutputDir}/{SeqID}_quant/quant.sf"
    quant_log = f"{OutputDir}/{SeqID}_quant/logs/salmon_quant.log"
    lib_format = f"{OutputDir}/{SeqID}_quant/lib_format_counts.json"

    # --- [Command Execution] ---
    cmd = f"salmon quant --threads {Threads} --index {SalmonIndex} --libType {libType} -1 {FastqDir}/{SeqID}_R1.{InputSuffix} -2 {FastqDir}/{SeqID}_R2.{InputSuffix} {validateMappings} -o {OutputDir}/{SeqID}_quant {extra_args}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    os.makedirs(os.path.dirname(FastqDir) if '.' in os.path.basename(FastqDir) else FastqDir, exist_ok=True)
    os.makedirs(os.path.dirname(OutputDir) if '.' in os.path.basename(OutputDir) else OutputDir, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()