#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = SRAToolkit
# VERSION = 3.0.10
# THREADS = 8
# PROFILE = SRA Download, FASTQ Conversion, Rename, and Gzip Compression

"""
Tool: SRAToolkit (3.0.10)
Profile: SRA Download, FASTQ Conversion, Rename, and Gzip Compression
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="SRAToolkit Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SraAccession', required=True, default='', help='다운로드할 SRA Run accession (Default: )')
    parser.add_argument('--OutDir', required=True, default='', help='최종 FASTQ.GZ 저장 디렉토리 (Default: )')
    parser.add_argument('--TmpDir', required=False, default='[OutDir]/tmp', help='임시 작업 디렉토리 (Default: [OutDir]/tmp)')
    parser.add_argument('--sra_dir', required=False, default='[TmpDir]/sra/[SraAccession]', help='prefetch 다운로드 디렉토리 (Default: [TmpDir]/sra/[SraAccession])')
    parser.add_argument('--r1_fastq', required=False, default='[OutDir]/[SraAccession]_R1.fastq', help='Read 1 FASTQ (Default: [OutDir]/[SraAccession]_R1.fastq)')
    parser.add_argument('--r2_fastq', required=False, default='[OutDir]/[SraAccession]_R2.fastq', help='Read 2 FASTQ (Default: [OutDir]/[SraAccession]_R2.fastq)')
    parser.add_argument('--r1_fastq_gz', required=False, default='[OutDir]/[SraAccession]_R1.fastq.gz', help='압축된 Read 1 FASTQ (Default: [OutDir]/[SraAccession]_R1.fastq.gz)')
    parser.add_argument('--r2_fastq_gz', required=False, default='[OutDir]/[SraAccession]_R2.fastq.gz', help='압축된 Read 2 FASTQ (Default: [OutDir]/[SraAccession]_R2.fastq.gz)')
    parser.add_argument('--prefetch_bin', required=False, default='/storage/apps/sratoolkit.3.0.10/bin/prefetch', help='prefetch 실행 파일 (Default: /storage/apps/sratoolkit.3.0.10/bin/prefetch)')
    parser.add_argument('--fasterq_dump_bin', required=False, default='/storage/apps/sratoolkit.3.0.10/bin/fasterq-dump', help='fasterq-dump 실행 파일 (Default: /storage/apps/sratoolkit.3.0.10/bin/fasterq-dump)')
    parser.add_argument('--Threads', required=False, default='8', help='사용할 스레드 수 (Default: 8)')
    parser.add_argument('--MaxSize', required=False, default='100G', help='prefetch 최대 다운로드 크기 (Default: 100G)')
    parser.add_argument('--SplitMode', required=False, default='--split-files', help='paired-end FASTQ 분리 옵션 (Default: --split-files)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SraAccession = args.SraAccession
    OutDir = args.OutDir
    TmpDir = args.TmpDir
    sra_dir = args.sra_dir
    r1_fastq = args.r1_fastq
    r2_fastq = args.r2_fastq
    r1_fastq_gz = args.r1_fastq_gz
    r2_fastq_gz = args.r2_fastq_gz
    prefetch_bin = args.prefetch_bin
    fasterq_dump_bin = args.fasterq_dump_bin
    Threads = args.Threads
    MaxSize = args.MaxSize
    SplitMode = args.SplitMode

    # --- [Output Paths] ---
    if not sra_dir:
    sra_dir = f"{TmpDir}/sra/{SraAccession}"
    if not r1_fastq:
    r1_fastq = f"{OutDir}/{SraAccession}_R1.fastq"
    if not r2_fastq:
    r2_fastq = f"{OutDir}/{SraAccession}_R2.fastq"
    if not r1_fastq_gz:
    r1_fastq_gz = f"{OutDir}/{SraAccession}_R1.fastq.gz"
    if not r2_fastq_gz:
    r2_fastq_gz = f"{OutDir}/{SraAccession}_R2.fastq.gz"

    # --- [Command Execution] ---
    cmd = f"mkdir -p {OutDir} {TmpDir}/sra {TmpDir}/fasterq && {prefetch_bin} --max-size {MaxSize} --output-directory {TmpDir}/sra {SraAccession} && {fasterq_dump_bin} {SplitMode} --threads {Threads} --temp {TmpDir}/fasterq --outdir {OutDir} {sra_dir} && mv {OutDir}/{SraAccession}_1.fastq {r1_fastq} && mv {OutDir}/{SraAccession}_2.fastq {r2_fastq} && gzip -f {r1_fastq} && gzip -f {r2_fastq}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    if r1_fastq:
        _tgt = os.path.dirname(r1_fastq) if os.path.splitext(r1_fastq)[1] else r1_fastq
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if OutDir:
        _tgt = os.path.dirname(OutDir) if os.path.splitext(OutDir)[1] else OutDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if sra_dir:
        _tgt = os.path.dirname(sra_dir) if os.path.splitext(sra_dir)[1] else sra_dir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if TmpDir:
        _tgt = os.path.dirname(TmpDir) if os.path.splitext(TmpDir)[1] else TmpDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if r2_fastq_gz:
        _tgt = os.path.dirname(r2_fastq_gz) if os.path.splitext(r2_fastq_gz)[1] else r2_fastq_gz
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if r1_fastq_gz:
        _tgt = os.path.dirname(r1_fastq_gz) if os.path.splitext(r1_fastq_gz)[1] else r1_fastq_gz
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if r2_fastq:
        _tgt = os.path.dirname(r2_fastq) if os.path.splitext(r2_fastq)[1] else r2_fastq
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()