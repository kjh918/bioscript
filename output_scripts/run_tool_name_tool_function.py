#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = TOOL_NAME
# VERSION = TASK_VERSION
# THREADS = 8
# PROFILE = TOOL_FUNCTION

"""
Tool: TOOL_NAME (TASK_VERSION)
Profile: TOOL_FUNCTION
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="TOOL_NAME Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='Unique sample identifier used to generate the output file names (Default: )')
    parser.add_argument('--InputFileDir', required=True, default='', help='Path of Input sample data. (Default: )')
    parser.add_argument('--qcResDir', required=True, default='', help='Path of QC Results. (Default: )')
    parser.add_argument('--OutputDir', required=True, default='', help='Path of Output files. (Default: )')
    parser.add_argument('--OutputSummary', required=True, default='', help='Path of Output files. (Default: )')
    parser.add_argument('--InputSuffix', required=False, default='', help='File name input-suffix (e.g., \'.bam\') (Default: )')
    parser.add_argument('--OutputSuffix', required=False, default='', help='File name output-suffix (e.g., \'.summary.txt\') (Default: )')
    parser.add_argument('--singularity_bin', required=False, default='singularity', help='No description (Default: singularity)')
    parser.add_argument('--bind', required=False, default='/storage,/data', help='No description (Default: /storage,/data)')
    parser.add_argument('--xmx_mb', required=False, default='16384', help='16GB memory from original script (Default: 16384)')
    parser.add_argument('--Threads', required=False, default='8', help='Threads (default : 8) (Default: 8)')
    parser.add_argument('--extra_args', required=False, default='', help='Extra arguments (Default: )')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    InputFileDir = args.InputFileDir
    qcResDir = args.qcResDir
    OutputDir = args.OutputDir
    OutputSummary = args.OutputSummary
    InputSuffix = args.InputSuffix
    OutputSuffix = args.OutputSuffix
    singularity_bin = args.singularity_bin
    bind = args.bind
    xmx_mb = args.xmx_mb
    Threads = args.Threads
    extra_args = args.extra_args

    # --- [Output Paths] ---

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind}  {extra_args}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    if InputFileDir:
        _tgt = os.path.dirname(InputFileDir) if os.path.splitext(InputFileDir)[1] else InputFileDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if OutputDir:
        _tgt = os.path.dirname(OutputDir) if os.path.splitext(OutputDir)[1] else OutputDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if OutputSummary:
        _tgt = os.path.dirname(OutputSummary) if os.path.splitext(OutputSummary)[1] else OutputSummary
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if qcResDir:
        _tgt = os.path.dirname(qcResDir) if os.path.splitext(qcResDir)[1] else qcResDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()