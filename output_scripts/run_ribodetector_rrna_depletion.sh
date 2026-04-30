#!/bin/bash
# [METADATA]
# TOOL_NAME = ribodetector
# VERSION = 0.3.3
# THREADS = 1

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
mkdir -p "$(dirname "${FastqDir}")" 2>/dev/null || mkdir -p "${FastqDir}"
mkdir -p "$(dirname "${OutputDir}")" 2>/dev/null || mkdir -p "${OutputDir}"

eval "$cmd"