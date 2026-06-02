#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = fastq_screen
# VERSION = 0.15.3
# THREADS = 15
# PROFILE = singularity_fastq_screen

"""
Tool: fastq_screen (0.15.3)
Profile: singularity_fastq_screen
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="fastq_screen Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='샘플 식별자 (Sample ID) (Default: )')
    parser.add_argument('--RawFastqDir', required=True, default='', help='원본 FASTQ 파일이 위치한 경로 (Default: )')
    parser.add_argument('--qcResDir', required=True, default='', help='결과 리포트가 저장될 디렉토리 (QC_DIR) (Default: )')
    parser.add_argument('--screen_txt', required=False, default='[qcResDir]/[SeqID]_R1_screen.txt', help='텍스트 형태의 매핑 요약 결과 (Default: [qcResDir]/[SeqID]_R1_screen.txt)')
    parser.add_argument('--screen_png', required=False, default='[qcResDir]/[SeqID]_R1_screen.png', help='시각화된 매핑 비율 차트 (Default: [qcResDir]/[SeqID]_R1_screen.png)')
    parser.add_argument('--screen_html', required=False, default='[qcResDir]/[SeqID]_R1_screen.html', help='인터랙티브 HTML 리포트 (Default: [qcResDir]/[SeqID]_R1_screen.html)')
    parser.add_argument('--singularity_bin', required=False, default='singularity', help='Singularity 실행 경로 (Default: singularity)')
    parser.add_argument('--fastq_screen_bin', required=False, default='fastq_screen', help='컨테이너 내부 실행 바이너리 명 (Default: fastq_screen)')
    parser.add_argument('--sif', required=False, default='/storage/images/fastqScreen-0.15.3.sif', help='SIF 이미지 경로 (Default: /storage/images/fastqScreen-0.15.3.sif)')
    parser.add_argument('--bind', required=False, default='/storage,/data', help='Singularity 바운드 경로 (Default: /storage,/data)')
    parser.add_argument('--aligner', required=False, default='bowtie2', help='사용할 Aligner (기본: bowtie2) (Default: bowtie2)')
    parser.add_argument('--config_path', required=False, default='/storage/home/kangsm/runScripts/NGS_config.FastqScreen.conf', help='FastQ Screen 설정 파일 경로 (Default: /storage/home/kangsm/runScripts/NGS_config.FastqScreen.conf)')
    parser.add_argument('--Threads', required=False, default='15', help='분석에 사용할 스레드 수 (Default: 15)')
    parser.add_argument('--extra_args', required=False, default='--force', help='추가 플래그 (예: --force, --subset 100000) (Default: --force)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    RawFastqDir = args.RawFastqDir
    qcResDir = args.qcResDir
    screen_txt = args.screen_txt
    screen_png = args.screen_png
    screen_html = args.screen_html
    singularity_bin = args.singularity_bin
    fastq_screen_bin = args.fastq_screen_bin
    sif = args.sif
    bind = args.bind
    aligner = args.aligner
    config_path = args.config_path
    Threads = args.Threads
    extra_args = args.extra_args

    # --- [Output Paths] ---
    if not screen_txt:
    screen_txt = f"{qcResDir}/{SeqID}_R1_screen.txt"
    if not screen_png:
    screen_png = f"{qcResDir}/{SeqID}_R1_screen.png"
    if not screen_html:
    screen_html = f"{qcResDir}/{SeqID}_R1_screen.html"

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind} {sif} {fastq_screen_bin} --aligner {aligner} --conf {config_path} --outdir {qcResDir} --threads {Threads} {extra_args} {RawFastqDir}/{SeqID}_R1.fastq.gz"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    if screen_html:
        _tgt = os.path.dirname(screen_html) if os.path.splitext(screen_html)[1] else screen_html
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if screen_txt:
        _tgt = os.path.dirname(screen_txt) if os.path.splitext(screen_txt)[1] else screen_txt
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if screen_png:
        _tgt = os.path.dirname(screen_png) if os.path.splitext(screen_png)[1] else screen_png
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if qcResDir:
        _tgt = os.path.dirname(qcResDir) if os.path.splitext(qcResDir)[1] else qcResDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if RawFastqDir:
        _tgt = os.path.dirname(RawFastqDir) if os.path.splitext(RawFastqDir)[1] else RawFastqDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()