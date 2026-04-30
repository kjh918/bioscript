#!/bin/bash
# [METADATA]
# TOOL_NAME = salmon_quant
# VERSION = 1.10.1
# THREADS = 1

# Tool Info: salmon_quant (1.10.1)
# Profile: transcript_quantification

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           No description"
    echo "  --FastqDir        No description"
    echo "  --OutputDir       정량화 결과(quant 폴더)가 저장될 경로"
    echo "  --SalmonIndex     Path to Salmon index directory"
    echo ""
    echo "Optional Parameters:"
    echo "  --Threads         No description (Default: 12)"
    echo "  --libType         Library type (A: Auto-detect) (Default: A)"
    echo "  --InputSuffix     No description (Default: fastq.gz)"
    echo "  --validateMappings No description (Default: --validateMappings)"
    echo "  --extra_args      No description (Default: )"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
FastqDir=""
OutputDir=""
SalmonIndex=""
Threads="12"
libType="A"
InputSuffix="fastq.gz"
validateMappings="--validateMappings"
extra_args=""

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --FastqDir) FastqDir="$2"; shift 2 ;;
        --OutputDir) OutputDir="$2"; shift 2 ;;
        --SalmonIndex) SalmonIndex="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        --libType) libType="$2"; shift 2 ;;
        --InputSuffix) InputSuffix="$2"; shift 2 ;;
        --validateMappings) validateMappings="$2"; shift 2 ;;
        --extra_args) extra_args="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${FastqDir:-}" ]]; then echo "Error: --FastqDir is required"; usage; fi
if [[ -z "${OutputDir:-}" ]]; then echo "Error: --OutputDir is required"; usage; fi
if [[ -z "${SalmonIndex:-}" ]]; then echo "Error: --SalmonIndex is required"; usage; fi

# --- [Output Paths] ---
quant_sf="${OutputDir}/${SeqID}_quant/quant.sf"
quant_log="${OutputDir}/${SeqID}_quant/logs/salmon_quant.log"
lib_format="${OutputDir}/${SeqID}_quant/lib_format_counts.json"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="salmon quant --threads ${Threads} --index ${SalmonIndex} --libType ${libType} -1 ${FastqDir}/${SeqID}_R1.${InputSuffix} -2 ${FastqDir}/${SeqID}_R2.${InputSuffix} ${validateMappings} -o ${OutputDir}/${SeqID}_quant ${extra_args}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
mkdir -p "$(dirname "${FastqDir}")" 2>/dev/null || mkdir -p "${FastqDir}"
mkdir -p "$(dirname "${OutputDir}")" 2>/dev/null || mkdir -p "${OutputDir}"

eval "$cmd"