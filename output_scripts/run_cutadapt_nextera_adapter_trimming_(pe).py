#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = Cutadapt
# VERSION = 4.4
# THREADS = 4
# PROFILE = Nextera Adapter Trimming (PE)

"""
Tool: Cutadapt (4.4)
Profile: Nextera Adapter Trimming (PE)
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="Cutadapt Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='분석 고유 ID (Sample ID) (Default: )')
    parser.add_argument('--RawFastqDir', required=True, default='', help='원본 FASTQ 파일이 위치한 디렉토리 경로 (Default: )')
    parser.add_argument('--OutDir', required=True, default='', help='Trimmed FASTQ 파일이 저장될 디렉토리 경로 (Default: )')
    parser.add_argument('--r1_clean', required=False, default='[OutDir]/[SeqID]_R1.trimmed.fastq.gz', help='Trimmed Read 1 FASTQ (Default: [OutDir]/[SeqID]_R1.trimmed.fastq.gz)')
    parser.add_argument('--r2_clean', required=False, default='[OutDir]/[SeqID]_R2.trimmed.fastq.gz', help='Trimmed Read 2 FASTQ (Default: [OutDir]/[SeqID]_R2.trimmed.fastq.gz)')
    parser.add_argument('--trim_log', required=False, default='[OutDir]/[SeqID]_cutadapt.log', help='Cutadapt 처리 통계 로그 파일 (Default: [OutDir]/[SeqID]_cutadapt.log)')
    parser.add_argument('--InputSuffix', required=False, default='.fastq.gz', help='입력 FASTQ 파일 확장자 (예: .fq.gz, .fastq.gz) (Default: .fastq.gz)')
    parser.add_argument('--AdapterR1', required=False, default='CTGTCTCTTATACACATCT', help='Read 1 3\' 어댑터 서열 (-a) (Default: CTGTCTCTTATACACATCT)')
    parser.add_argument('--AdapterR2', required=False, default='CTGTCTCTTATACACATCT', help='Read 2 3\' 어댑터 서열 (-A) (Default: CTGTCTCTTATACACATCT)')
    parser.add_argument('--QualityCutoff', required=False, default='20', help='3\' 말단 Quality Trimming 임계값 (-q) (Default: 20)')
    parser.add_argument('--MinLength', required=False, default='30', help='Trimming 후 남길 최소 Read 길이 (-m) (Default: 30)')
    parser.add_argument('--Threads', required=False, default='4', help='병렬 처리에 사용할 스레드 수 (-j) (Default: 4)')
    parser.add_argument('--cutadapt_bin', required=False, default='/storage/home/jhkim/Apps/Python-3.11.13/python -m cutadapt', help='Cutadapt 실행 명령어 또는 절대 경로 (Default: /storage/home/jhkim/Apps/Python-3.11.13/python -m cutadapt)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    RawFastqDir = args.RawFastqDir
    OutDir = args.OutDir
    r1_clean = args.r1_clean
    r2_clean = args.r2_clean
    trim_log = args.trim_log
    InputSuffix = args.InputSuffix
    AdapterR1 = args.AdapterR1
    AdapterR2 = args.AdapterR2
    QualityCutoff = args.QualityCutoff
    MinLength = args.MinLength
    Threads = args.Threads
    cutadapt_bin = args.cutadapt_bin

    # --- [Output Paths] ---
    if not r1_clean:
    r1_clean = f"{OutDir}/{SeqID}_R1.trimmed.fastq.gz"
    if not r2_clean:
    r2_clean = f"{OutDir}/{SeqID}_R2.trimmed.fastq.gz"
    if not trim_log:
    trim_log = f"{OutDir}/{SeqID}_cutadapt.log"

    # --- [Command Execution] ---
    cmd = f"{cutadapt_bin} -j {Threads} -a {AdapterR1} -A {AdapterR2} -q {QualityCutoff} -m {MinLength} -o {r1_clean} -p {r2_clean} {RawFastqDir}/{SeqID}_R1{InputSuffix} {RawFastqDir}/{SeqID}_R2{InputSuffix} > {trim_log}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    if trim_log:
        _tgt = os.path.dirname(trim_log) if os.path.splitext(trim_log)[1] else trim_log
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if RawFastqDir:
        _tgt = os.path.dirname(RawFastqDir) if os.path.splitext(RawFastqDir)[1] else RawFastqDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if OutDir:
        _tgt = os.path.dirname(OutDir) if os.path.splitext(OutDir)[1] else OutDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if r2_clean:
        _tgt = os.path.dirname(r2_clean) if os.path.splitext(r2_clean)[1] else r2_clean
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if r1_clean:
        _tgt = os.path.dirname(r1_clean) if os.path.splitext(r1_clean)[1] else r1_clean
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()