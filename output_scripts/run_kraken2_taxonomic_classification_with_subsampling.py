#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = Kraken2
# VERSION = 2.17.1
# THREADS = 8
# PROFILE = Taxonomic Classification with Subsampling

"""
Tool: Kraken2 (2.17.1)
Profile: Taxonomic Classification with Subsampling
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="Kraken2 Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='분석 고유 ID (Sample ID) (Default: )')
    parser.add_argument('--TrimFastqDir', required=True, default='', help='Trimmed FASTQ 파일이 위치한 디렉토리 경로 (Default: )')
    parser.add_argument('--KrakenDB', required=False, default='/storage/home/jhkim/Apps/kraken2/bin/kraken2_db', help='Kraken2-db 명령어 (Default: /storage/home/jhkim/Apps/kraken2/bin/kraken2_db)')
    parser.add_argument('--OutDir', required=True, default='', help='결과물이 저장될 디렉토리 경로 (Default: )')
    parser.add_argument('--sub_r1', required=False, default='[OutDir]/[SeqID].trimmed_sub_R1.fastq.gz', help='Subsampled Read 1 (Default: [OutDir]/[SeqID].trimmed_sub_R1.fastq.gz)')
    parser.add_argument('--sub_r2', required=False, default='[OutDir]/[SeqID].trimmed_sub_R2.fastq.gz', help='Subsampled Read 2 (Default: [OutDir]/[SeqID].trimmed_sub_R2.fastq.gz)')
    parser.add_argument('--kraken_out', required=False, default='[OutDir]/[SeqID].kraken', help='Kraken2 분류 결과 파일 (Default: [OutDir]/[SeqID].kraken)')
    parser.add_argument('--report_out', required=False, default='[OutDir]/[SeqID].kreport', help='Kraken2 계층적 요약 리포트 (Default: [OutDir]/[SeqID].kreport)')
    parser.add_argument('--plot_out', required=False, default='[OutDir]/[SeqID]_krona.html', help='생성된 Krona 대화형 HTML 플롯 (Default: [OutDir]/[SeqID]_krona.html)')
    parser.add_argument('--seqkit_bin', required=False, default='/storage/apps/seqkit-2.8.2/seqkit', help='Seqkit 실행 절대 경로 (Default: /storage/apps/seqkit-2.8.2/seqkit)')
    parser.add_argument('--SampleSeed', required=False, default='100', help='난수 시드 (PE 동기화를 위해 R1/R2 동일 적용) (Default: 100)')
    parser.add_argument('--SampleNum', required=False, default='100000', help='추출할 Read 개수 (-n) (Default: 100000)')
    parser.add_argument('--R1_suffix', required=False, default='.trimmed_R1.fastq.gz', help='Read 1 입력 파일의 Suffix (Default: .trimmed_R1.fastq.gz)')
    parser.add_argument('--R2_suffix', required=False, default='.trimmed_R2.fastq.gz', help='Read 2 입력 파일의 Suffix (Default: .trimmed_R2.fastq.gz)')
    parser.add_argument('--kraken2_bin', required=False, default='/storage/home/jhkim/Apps/kraken2/bin/kraken2', help='Kraken2 명령어 (Default: /storage/home/jhkim/Apps/kraken2/bin/kraken2)')
    parser.add_argument('--krona_bin', required=False, default='/storage/home/jhkim/Apps/Krona/KronaTools/bin/ktImportTaxonomy', help='Krona 명령어 (Default: /storage/home/jhkim/Apps/Krona/KronaTools/bin/ktImportTaxonomy)')
    parser.add_argument('--Threads', required=False, default='8', help='병렬 처리에 사용할 스레드 수 (Default: 8)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    TrimFastqDir = args.TrimFastqDir
    KrakenDB = args.KrakenDB
    OutDir = args.OutDir
    sub_r1 = args.sub_r1
    sub_r2 = args.sub_r2
    kraken_out = args.kraken_out
    report_out = args.report_out
    plot_out = args.plot_out
    seqkit_bin = args.seqkit_bin
    SampleSeed = args.SampleSeed
    SampleNum = args.SampleNum
    R1_suffix = args.R1_suffix
    R2_suffix = args.R2_suffix
    kraken2_bin = args.kraken2_bin
    krona_bin = args.krona_bin
    Threads = args.Threads

    # --- [Output Paths] ---
    if not sub_r1:
    sub_r1 = f"{OutDir}/{SeqID}.trimmed_sub_R1.fastq.gz"
    if not sub_r2:
    sub_r2 = f"{OutDir}/{SeqID}.trimmed_sub_R2.fastq.gz"
    if not kraken_out:
    kraken_out = f"{OutDir}/{SeqID}.kraken"
    if not report_out:
    report_out = f"{OutDir}/{SeqID}.kreport"
    if not plot_out:
    plot_out = f"{OutDir}/{SeqID}_krona.html"

    # --- [Command Execution] ---
    cmd = f"{seqkit_bin} sample -s {SampleSeed} -n {SampleNum} {TrimFastqDir}/{SeqID}{R1_suffix} | gzip > {sub_r1} && {seqkit_bin} sample -s {SampleSeed} -n {SampleNum} {TrimFastqDir}/{SeqID}{R2_suffix} | gzip > {sub_r2} && {kraken2_bin} --db {KrakenDB} --threads {Threads} --report {report_out} --paired {sub_r1} {sub_r2} > {kraken_out} && {krona_bin} -m 4 -t 3 {kraken_out} -o {plot_out}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    if sub_r1:
        _tgt = os.path.dirname(sub_r1) if os.path.splitext(sub_r1)[1] else sub_r1
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if kraken_out:
        _tgt = os.path.dirname(kraken_out) if os.path.splitext(kraken_out)[1] else kraken_out
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if report_out:
        _tgt = os.path.dirname(report_out) if os.path.splitext(report_out)[1] else report_out
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if plot_out:
        _tgt = os.path.dirname(plot_out) if os.path.splitext(plot_out)[1] else plot_out
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if OutDir:
        _tgt = os.path.dirname(OutDir) if os.path.splitext(OutDir)[1] else OutDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if TrimFastqDir:
        _tgt = os.path.dirname(TrimFastqDir) if os.path.splitext(TrimFastqDir)[1] else TrimFastqDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if sub_r2:
        _tgt = os.path.dirname(sub_r2) if os.path.splitext(sub_r2)[1] else sub_r2
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()