#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = TOOL_NAME
# VERSION = VERSION
# THREADS = 8
# PROFILE = TASK_PROFILE

"""
Tool: TOOL_NAME (VERSION)
Profile: TASK_PROFILE
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="TOOL_NAME Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='분석 고유 ID (Default: )')
    parser.add_argument('--InputDir', required=True, default='', help='입력 파일 디렉토리 (Default: )')
    parser.add_argument('--OutDir', required=True, default='', help='결과 저장 디렉토리 (Default: )')
    parser.add_argument('--output_file', required=False, default='[OutDir]/[SeqID].out', help='주요 결과 파일 (Default: [OutDir]/[SeqID].out)')
    parser.add_argument('--tool_bin', required=False, default='TOOL_COMMAND', help='Tool path (Default: TOOL_COMMAND)')
    parser.add_argument('--InputSuffix', required=False, default='.UseAnalysis.bam', help='Input bam suffix (Default: .UseAnalysis.bam)')
    parser.add_argument('--OutputSuffix', required=False, default='.txt', help='Output txt suffix (Default: .txt)')
    parser.add_argument('--Threads', required=False, default='8', help='사용할 스레드 수 (Default: 8)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    InputDir = args.InputDir
    OutDir = args.OutDir
    output_file = args.output_file
    tool_bin = args.tool_bin
    InputSuffix = args.InputSuffix
    OutputSuffix = args.OutputSuffix
    Threads = args.Threads

    # --- [Output Paths] ---
    if not output_file:
    output_file = f"{OutDir}/{SeqID}.out"

    # --- [Command Execution] ---
    cmd = f"{tool_bin} INPUT_OPTIONS OUTPUT_OPTIONS"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    if output_file:
        _tgt = os.path.dirname(output_file) if os.path.splitext(output_file)[1] else output_file
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if OutDir:
        _tgt = os.path.dirname(OutDir) if os.path.splitext(OutDir)[1] else OutDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if InputDir:
        _tgt = os.path.dirname(InputDir) if os.path.splitext(InputDir)[1] else InputDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()