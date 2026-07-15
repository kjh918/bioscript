#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = cbnipt_cnv_call_cli
# VERSION = 1.0.0
# THREADS = 1
# PROFILE = cnv_calling

"""
Tool: cbnipt_cnv_call_cli (1.0.0)
Profile: cnv_calling
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="cbnipt_cnv_call_cli Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--ReferenceFasta', required=True, default='', help='Absolute path to the hg38 reference genome FASTA file. (Default: )')
    parser.add_argument('--MappabilityBW', required=False, default='/storage/home/jhkim/Apps/Python-3.11.13/python', help='Absolute path to the target Python interpreter environment. (Default: /storage/home/jhkim/Apps/Python-3.11.13/python)')
    parser.add_argument('--OutBinFile', required=True, default='', help='Absolute path to the specific run output directory. (Default: )')
    parser.add_argument('--CliScript', required=False, default='/storage/home/jhkim/scripts/bioscript/manual/cbnipt_manual_cnv_call/scripts/cli.py', help='Absolute path to the target cli.py script tool. (Default: /storage/home/jhkim/scripts/bioscript/manual/cbnipt_manual_cnv_call/scripts/cli.py)')
    parser.add_argument('--BinSize', type=int, default=100000, help='Bin 크기 (100kb)')
    parser.add_argument('--MinMappability', type=float, default=0.9, help='최소 Mappability 점수')
    parser.add_argument('--MinGC', type=float, default=0.3, help='최소 GC 함량 필터')
    parser.add_argument('--MaxGC', type=float, default=0.7, help='최대 GC 함량 필터')
    parser.add_argument('--Threads', required=False, default='1', help='The parsed identifier extracted from the BAM filename. (Default: 4)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    BamPath = args.BamPath
    ReferenceFasta = args.ReferenceFasta
    OutDir = args.OutDir
    PythonBin = args.PythonBin
    CliScript = args.CliScript
    Threads = args.Threads


    # --- [Command Execution] ---
    cmd = f"{PythonBin} {CliScript} make-bins --BamPath {BamPath} --ReferenceFasta {ReferenceFasta} --OutDir {OutDir} --Threads {Threads}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    os.makedirs(os.path.dirname(OutDir) if '.' in os.path.basename(OutDir) else OutDir, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()