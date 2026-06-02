#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = sortmerna
# VERSION = 4.3.7
# THREADS = 8
# PROFILE = rrna_filtering

"""
Tool: sortmerna (4.3.7)
Profile: rrna_filtering
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="sortmerna Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='No description (Default: )')
    parser.add_argument('--FastqDir', required=True, default='', help='No description (Default: )')
    parser.add_argument('--qcResDir', required=True, default='', help='Log 및 결과 파일 저장 경로 (Default: )')
    parser.add_argument('--RefArgs', required=True, default='', help='List of --ref references (Default: )')
    parser.add_argument('--IndexDir', required=False, default='', help='Pre-built index directory (--idx-dir) (Default: )')
    parser.add_argument('--non_rrna_r1', required=False, default='[qcResDir]/[SeqID]_1.non_rRNA.fastq.gz', help='No description (Default: [qcResDir]/[SeqID]_1.non_rRNA.fastq.gz)')
    parser.add_argument('--non_rrna_r2', required=False, default='[qcResDir]/[SeqID]_2.non_rRNA.fastq.gz', help='No description (Default: [qcResDir]/[SeqID]_2.non_rRNA.fastq.gz)')
    parser.add_argument('--sortmerna_log', required=False, default='[qcResDir]/[SeqID].sortmerna.log', help='No description (Default: [qcResDir]/[SeqID].sortmerna.log)')
    parser.add_argument('--Threads', required=False, default='8', help='No description (Default: 8)')
    parser.add_argument('--InputSuffix', required=False, default='fastq.gz', help='No description (Default: fastq.gz)')
    parser.add_argument('--paired_cmd', required=False, default='--paired_in --out2', help='No description (Default: --paired_in --out2)')
    parser.add_argument('--extra_args', required=False, default='', help='No description (Default: )')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    FastqDir = args.FastqDir
    qcResDir = args.qcResDir
    RefArgs = args.RefArgs
    IndexDir = args.IndexDir
    non_rrna_r1 = args.non_rrna_r1
    non_rrna_r2 = args.non_rrna_r2
    sortmerna_log = args.sortmerna_log
    Threads = args.Threads
    InputSuffix = args.InputSuffix
    paired_cmd = args.paired_cmd
    extra_args = args.extra_args

    # --- [Output Paths] ---
    if not non_rrna_r1:
    non_rrna_r1 = f"{qcResDir}/{SeqID}_1.non_rRNA.fastq.gz"
    if not non_rrna_r2:
    non_rrna_r2 = f"{qcResDir}/{SeqID}_2.non_rRNA.fastq.gz"
    if not sortmerna_log:
    sortmerna_log = f"{qcResDir}/{SeqID}.sortmerna.log"

    # --- [Command Execution] ---
    cmd = f"sortmerna {RefArgs} --reads {FastqDir}/{SeqID}_1.{InputSuffix} --reads {FastqDir}/{SeqID}_2.{InputSuffix} --threads {Threads} --workdir . --aligned {qcResDir}/rRNA_reads --fastx --other {qcResDir}/non_rRNA_reads {paired_cmd} {extra_args} && mv {qcResDir}/non_rRNA_reads_fwd.f*q.gz {non_rrna_r1} && mv {qcResDir}/non_rRNA_reads_rev.f*q.gz {non_rrna_r2} && mv {qcResDir}/rRNA_reads.log {sortmerna_log}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    if non_rrna_r1:
        _tgt = os.path.dirname(non_rrna_r1) if os.path.splitext(non_rrna_r1)[1] else non_rrna_r1
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if qcResDir:
        _tgt = os.path.dirname(qcResDir) if os.path.splitext(qcResDir)[1] else qcResDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if FastqDir:
        _tgt = os.path.dirname(FastqDir) if os.path.splitext(FastqDir)[1] else FastqDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if non_rrna_r2:
        _tgt = os.path.dirname(non_rrna_r2) if os.path.splitext(non_rrna_r2)[1] else non_rrna_r2
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if sortmerna_log:
        _tgt = os.path.dirname(sortmerna_log) if os.path.splitext(sortmerna_log)[1] else sortmerna_log
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if IndexDir:
        _tgt = os.path.dirname(IndexDir) if os.path.splitext(IndexDir)[1] else IndexDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()