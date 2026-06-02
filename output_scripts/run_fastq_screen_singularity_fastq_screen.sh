#!/bin/bash
# [METADATA]
# TOOL_NAME = fastq_screen
# VERSION = 0.15.3
# THREADS = 15

# Tool Info: fastq_screen (0.15.3)
# Profile: singularity_fastq_screen

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           샘플 식별자 (Sample ID)"
    echo "  --RawFastqDir     원본 FASTQ 파일이 위치한 경로"
    echo "  --qcResDir        결과 리포트가 저장될 디렉토리 (QC_DIR)"
    echo ""
    echo "Optional Parameters:"
    echo "  --screen_txt      텍스트 형태의 매핑 요약 결과 (Default: [qcResDir]/[SeqID]_R1_screen.txt)"
    echo "  --screen_png      시각화된 매핑 비율 차트 (Default: [qcResDir]/[SeqID]_R1_screen.png)"
    echo "  --screen_html     인터랙티브 HTML 리포트 (Default: [qcResDir]/[SeqID]_R1_screen.html)"
    echo "  --singularity_bin Singularity 실행 경로 (Default: singularity)"
    echo "  --fastq_screen_bin 컨테이너 내부 실행 바이너리 명 (Default: fastq_screen)"
    echo "  --sif             SIF 이미지 경로 (Default: /storage/images/fastqScreen-0.15.3.sif)"
    echo "  --bind            Singularity 바운드 경로 (Default: /storage,/data)"
    echo "  --aligner         사용할 Aligner (기본: bowtie2) (Default: bowtie2)"
    echo "  --config_path     FastQ Screen 설정 파일 경로 (Default: /storage/home/kangsm/runScripts/NGS_config.FastqScreen.conf)"
    echo "  --Threads         분석에 사용할 스레드 수 (Default: 15)"
    echo "  --extra_args      추가 플래그 (예: --force, --subset 100000) (Default: --force)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
RawFastqDir=""
qcResDir=""
screen_txt="[qcResDir]/[SeqID]_R1_screen.txt"
screen_png="[qcResDir]/[SeqID]_R1_screen.png"
screen_html="[qcResDir]/[SeqID]_R1_screen.html"
singularity_bin="singularity"
fastq_screen_bin="fastq_screen"
sif="/storage/images/fastqScreen-0.15.3.sif"
bind="/storage,/data"
aligner="bowtie2"
config_path="/storage/home/kangsm/runScripts/NGS_config.FastqScreen.conf"
Threads="15"
extra_args="--force"

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --RawFastqDir) RawFastqDir="$2"; shift 2 ;;
        --qcResDir) qcResDir="$2"; shift 2 ;;
        --screen_txt) screen_txt="$2"; shift 2 ;;
        --screen_png) screen_png="$2"; shift 2 ;;
        --screen_html) screen_html="$2"; shift 2 ;;
        --singularity_bin) singularity_bin="$2"; shift 2 ;;
        --fastq_screen_bin) fastq_screen_bin="$2"; shift 2 ;;
        --sif) sif="$2"; shift 2 ;;
        --bind) bind="$2"; shift 2 ;;
        --aligner) aligner="$2"; shift 2 ;;
        --config_path) config_path="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        --extra_args) extra_args="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${RawFastqDir:-}" ]]; then echo "Error: --RawFastqDir is required"; usage; fi
if [[ -z "${qcResDir:-}" ]]; then echo "Error: --qcResDir is required"; usage; fi

# --- [Output Paths] ---
screen_txt="${qcResDir}/${SeqID}_R1_screen.txt"
screen_png="${qcResDir}/${SeqID}_R1_screen.png"
screen_html="${qcResDir}/${SeqID}_R1_screen.html"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${singularity_bin} exec -B ${bind} ${sif} ${fastq_screen_bin} --aligner ${aligner} --conf ${config_path} --outdir ${qcResDir} --threads ${Threads} ${extra_args} ${RawFastqDir}/${SeqID}_R1.fastq.gz"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
if [[ -n "${screen_html:-}" ]]; then
  if [[ "${screen_html}" == *.* ]]; then mkdir -p "$(dirname "${screen_html}")"; else mkdir -p "${screen_html}"; fi
fi
if [[ -n "${screen_txt:-}" ]]; then
  if [[ "${screen_txt}" == *.* ]]; then mkdir -p "$(dirname "${screen_txt}")"; else mkdir -p "${screen_txt}"; fi
fi
if [[ -n "${screen_png:-}" ]]; then
  if [[ "${screen_png}" == *.* ]]; then mkdir -p "$(dirname "${screen_png}")"; else mkdir -p "${screen_png}"; fi
fi
if [[ -n "${qcResDir:-}" ]]; then
  if [[ "${qcResDir}" == *.* ]]; then mkdir -p "$(dirname "${qcResDir}")"; else mkdir -p "${qcResDir}"; fi
fi
if [[ -n "${RawFastqDir:-}" ]]; then
  if [[ "${RawFastqDir}" == *.* ]]; then mkdir -p "$(dirname "${RawFastqDir}")"; else mkdir -p "${RawFastqDir}"; fi
fi

eval "$cmd"