#!/bin/bash
# [METADATA]
# TOOL_NAME = Cutadapt
# VERSION = 4.4
# THREADS = 4

# Tool Info: Cutadapt (4.4)
# Profile: Nextera Adapter Trimming (PE)

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           분석 고유 ID (Sample ID)"
    echo "  --RawFastqDir     원본 FASTQ 파일이 위치한 디렉토리 경로"
    echo "  --OutDir          Trimmed FASTQ 파일이 저장될 디렉토리 경로"
    echo ""
    echo "Optional Parameters:"
    echo "  --r1_clean        Trimmed Read 1 FASTQ (Default: [OutDir]/[SeqID]_R1.trimmed.fastq.gz)"
    echo "  --r2_clean        Trimmed Read 2 FASTQ (Default: [OutDir]/[SeqID]_R2.trimmed.fastq.gz)"
    echo "  --trim_log        Cutadapt 처리 통계 로그 파일 (Default: [OutDir]/[SeqID]_cutadapt.log)"
    echo "  --InputSuffix     입력 FASTQ 파일 확장자 (예: .fq.gz, .fastq.gz) (Default: .fastq.gz)"
    echo "  --AdapterR1       Read 1 3' 어댑터 서열 (-a) (Default: CTGTCTCTTATACACATCT)"
    echo "  --AdapterR2       Read 2 3' 어댑터 서열 (-A) (Default: CTGTCTCTTATACACATCT)"
    echo "  --QualityCutoff   3' 말단 Quality Trimming 임계값 (-q) (Default: 20)"
    echo "  --MinLength       Trimming 후 남길 최소 Read 길이 (-m) (Default: 30)"
    echo "  --Threads         병렬 처리에 사용할 스레드 수 (-j) (Default: 4)"
    echo "  --cutadapt_bin    Cutadapt 실행 명령어 또는 절대 경로 (Default: /storage/home/jhkim/Apps/Python-3.11.13/python -m cutadapt)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
RawFastqDir=""
OutDir=""
r1_clean="[OutDir]/[SeqID]_R1.trimmed.fastq.gz"
r2_clean="[OutDir]/[SeqID]_R2.trimmed.fastq.gz"
trim_log="[OutDir]/[SeqID]_cutadapt.log"
InputSuffix=".fastq.gz"
AdapterR1="CTGTCTCTTATACACATCT"
AdapterR2="CTGTCTCTTATACACATCT"
QualityCutoff="20"
MinLength="30"
Threads="4"
cutadapt_bin="/storage/home/jhkim/Apps/Python-3.11.13/python -m cutadapt"

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --RawFastqDir) RawFastqDir="$2"; shift 2 ;;
        --OutDir) OutDir="$2"; shift 2 ;;
        --r1_clean) r1_clean="$2"; shift 2 ;;
        --r2_clean) r2_clean="$2"; shift 2 ;;
        --trim_log) trim_log="$2"; shift 2 ;;
        --InputSuffix) InputSuffix="$2"; shift 2 ;;
        --AdapterR1) AdapterR1="$2"; shift 2 ;;
        --AdapterR2) AdapterR2="$2"; shift 2 ;;
        --QualityCutoff) QualityCutoff="$2"; shift 2 ;;
        --MinLength) MinLength="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        --cutadapt_bin) cutadapt_bin="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${RawFastqDir:-}" ]]; then echo "Error: --RawFastqDir is required"; usage; fi
if [[ -z "${OutDir:-}" ]]; then echo "Error: --OutDir is required"; usage; fi

# --- [Output Paths] ---
r1_clean="${OutDir}/${SeqID}_R1.trimmed.fastq.gz"
r2_clean="${OutDir}/${SeqID}_R2.trimmed.fastq.gz"
trim_log="${OutDir}/${SeqID}_cutadapt.log"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${cutadapt_bin} -j ${Threads} -a ${AdapterR1} -A ${AdapterR2} -q ${QualityCutoff} -m ${MinLength} -o ${r1_clean} -p ${r2_clean} ${RawFastqDir}/${SeqID}_R1${InputSuffix} ${RawFastqDir}/${SeqID}_R2${InputSuffix} > ${trim_log}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
if [[ -n "${trim_log:-}" ]]; then
  if [[ "${trim_log}" == *.* ]]; then mkdir -p "$(dirname "${trim_log}")"; else mkdir -p "${trim_log}"; fi
fi
if [[ -n "${RawFastqDir:-}" ]]; then
  if [[ "${RawFastqDir}" == *.* ]]; then mkdir -p "$(dirname "${RawFastqDir}")"; else mkdir -p "${RawFastqDir}"; fi
fi
if [[ -n "${OutDir:-}" ]]; then
  if [[ "${OutDir}" == *.* ]]; then mkdir -p "$(dirname "${OutDir}")"; else mkdir -p "${OutDir}"; fi
fi
if [[ -n "${r2_clean:-}" ]]; then
  if [[ "${r2_clean}" == *.* ]]; then mkdir -p "$(dirname "${r2_clean}")"; else mkdir -p "${r2_clean}"; fi
fi
if [[ -n "${r1_clean:-}" ]]; then
  if [[ "${r1_clean}" == *.* ]]; then mkdir -p "$(dirname "${r1_clean}")"; else mkdir -p "${r1_clean}"; fi
fi

eval "$cmd"