#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = ribodetector
# VERSION = 0.3.3
# THREADS = 1
# PROFILE = rrna_depletion

"""
Tool: ribodetector (0.3.3)
Profile: rrna_depletion
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="ribodetector Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='No description (Default: )')
    parser.add_argument('--FastqDir', required=True, default='', help='No description (Default: )')
    parser.add_argument('--OutputDir', required=True, default='', help='결과 FASTQ가 저장될 경로 (Default: )')
    parser.add_argument('--Threads', required=False, default='10', help='No description (Default: 10)')
    parser.add_argument('--ReadLen', required=False, default='100', help='시퀀싱 리드 길이 (Default: 100)')
    parser.add_argument('--ChunkSize', required=False, default='256', help='모델 로딩 및 처리 청크 크기 (Default: 256)')
    parser.add_argument('--ExcludeMode', required=False, default='rrna', help='제거할 타입 (rrna) (Default: rrna)')
    parser.add_argument('--InputSuffix', required=False, default='fastq.gz', help='No description (Default: fastq.gz)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    FastqDir = args.FastqDir
    OutputDir = args.OutputDir
    Threads = args.Threads
    ReadLen = args.ReadLen
    ChunkSize = args.ChunkSize
    ExcludeMode = args.ExcludeMode
    InputSuffix = args.InputSuffix

    # --- [Output Paths] ---
    non_rrna_r1 = f"{OutputDir}/{SeqID}.nonrrna.1.fq"
    non_rrna_r2 = f"{OutputDir}/{SeqID}.nonrrna.2.fq"

    # --- [Command Execution] ---
    cmd = f"ribodetector_cpu -t {Threads} -l {ReadLen} -i {FastqDir}/{SeqID}_R1.{InputSuffix} {FastqDir}/{SeqID}_R2.{InputSuffix} -e {ExcludeMode} --chunk_size {ChunkSize} -o {non_rrna_r1} {non_rrna_r2}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    os.makedirs(os.path.dirname(FastqDir) if '.' in os.path.basename(FastqDir) else FastqDir, exist_ok=True)
    os.makedirs(os.path.dirname(OutputDir) if '.' in os.path.basename(OutputDir) else OutputDir, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()