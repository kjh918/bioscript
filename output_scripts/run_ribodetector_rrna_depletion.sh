#!/bin/bash
# [METADATA]
# TOOL_NAME = ribodetector
# VERSION = 0.3.3
# THREADS = 10

# Tool Info: ribodetector (0.3.3)
# Profile: rrna_depletion

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           No description"
    echo "  --FastqDir        No description"
    echo "  --OutputDir       결과 FASTQ가 저장될 경로"
    echo ""
    echo "Optional Parameters:"
    echo "  --non_rrna_r1     No description (Default: [OutputDir]/[SeqID].nonrrna.1.fq)"
    echo "  --non_rrna_r2     No description (Default: [OutputDir]/[SeqID].nonrrna.2.fq)"
    echo "  --Threads         No description (Default: 10)"
    echo "  --ReadLen         시퀀싱 리드 길이 (Default: 100)"
    echo "  --ChunkSize       모델 로딩 및 처리 청크 크기 (Default: 256)"
    echo "  --ExcludeMode     제거할 타입 (rrna) (Default: rrna)"
    echo "  --InputSuffix     No description (Default: fastq.gz)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
FastqDir=""
OutputDir=""
non_rrna_r1="[OutputDir]/[SeqID].nonrrna.1.fq"
non_rrna_r2="[OutputDir]/[SeqID].nonrrna.2.fq"
Threads="10"
ReadLen="100"
ChunkSize="256"
ExcludeMode="rrna"
InputSuffix="fastq.gz"

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --FastqDir) FastqDir="$2"; shift 2 ;;
        --OutputDir) OutputDir="$2"; shift 2 ;;
        --non_rrna_r1) non_rrna_r1="$2"; shift 2 ;;
        --non_rrna_r2) non_rrna_r2="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        --ReadLen) ReadLen="$2"; shift 2 ;;
        --ChunkSize) ChunkSize="$2"; shift 2 ;;
        --ExcludeMode) ExcludeMode="$2"; shift 2 ;;
        --InputSuffix) InputSuffix="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${FastqDir:-}" ]]; then echo "Error: --FastqDir is required"; usage; fi
if [[ -z "${OutputDir:-}" ]]; then echo "Error: --OutputDir is required"; usage; fi

# --- [Output Paths] ---
non_rrna_r1="${OutputDir}/${SeqID}.nonrrna.1.fq"
non_rrna_r2="${OutputDir}/${SeqID}.nonrrna.2.fq"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="ribodetector_cpu -t ${Threads} -l ${ReadLen} -i ${FastqDir}/${SeqID}_R1.${InputSuffix} ${FastqDir}/${SeqID}_R2.${InputSuffix} -e ${ExcludeMode} --chunk_size ${ChunkSize} -o ${non_rrna_r1} ${non_rrna_r2}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
if [[ -n "${non_rrna_r1:-}" ]]; then
  if [[ "${non_rrna_r1}" == *.* ]]; then mkdir -p "$(dirname "${non_rrna_r1}")"; else mkdir -p "${non_rrna_r1}"; fi
fi
if [[ -n "${non_rrna_r2:-}" ]]; then
  if [[ "${non_rrna_r2}" == *.* ]]; then mkdir -p "$(dirname "${non_rrna_r2}")"; else mkdir -p "${non_rrna_r2}"; fi
fi
if [[ -n "${OutputDir:-}" ]]; then
  if [[ "${OutputDir}" == *.* ]]; then mkdir -p "$(dirname "${OutputDir}")"; else mkdir -p "${OutputDir}"; fi
fi
if [[ -n "${FastqDir:-}" ]]; then
  if [[ "${FastqDir}" == *.* ]]; then mkdir -p "$(dirname "${FastqDir}")"; else mkdir -p "${FastqDir}"; fi
fi

eval "$cmd"