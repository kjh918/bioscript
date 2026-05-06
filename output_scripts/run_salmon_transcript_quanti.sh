#!/bin/bash
# [METADATA]
# TOOL_NAME = salmon
# VERSION = 1.10.0
# THREADS = 1

# Tool Info: salmon (1.10.0)
# Profile: transcript_quanti

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           No description"
    echo "  --BamDir          Directory containing input BAM files"
    echo "  --OutputDir       Path to output directory"
    echo "  --SalmonIndex     Path to Salmon index directory"
    echo ""
    echo "Optional Parameters:"
    echo "  --Threads         No description (Default: 12)"
    echo "  --libType         Library type (A: Auto-detect) (Default: A)"
    echo "  --InputSuffix     No description (Default: .Aligned.toTranscriptome.out.bam)"
    echo "  --extra_args      No description (Default: )"
    echo "  --samlon_bin      No description (Default: /storage/apps/salmon-1.10.0/bin/salmon)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
BamDir=""
OutputDir=""
SalmonIndex=""
Threads="12"
libType="A"
InputSuffix=".Aligned.toTranscriptome.out.bam"
extra_args=""
samlon_bin="/storage/apps/salmon-1.10.0/bin/salmon"

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --BamDir) BamDir="$2"; shift 2 ;;
        --OutputDir) OutputDir="$2"; shift 2 ;;
        --SalmonIndex) SalmonIndex="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        --libType) libType="$2"; shift 2 ;;
        --InputSuffix) InputSuffix="$2"; shift 2 ;;
        --extra_args) extra_args="$2"; shift 2 ;;
        --samlon_bin) samlon_bin="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${BamDir:-}" ]]; then echo "Error: --BamDir is required"; usage; fi
if [[ -z "${OutputDir:-}" ]]; then echo "Error: --OutputDir is required"; usage; fi
if [[ -z "${SalmonIndex:-}" ]]; then echo "Error: --SalmonIndex is required"; usage; fi

# --- [Output Paths] ---
quant_sf="${OutputDir}/${SeqID}_quant/quant.sf"
quant_log="${OutputDir}/${SeqID}_quant/logs/salmon_quant.log"
lib_format="${OutputDir}/${SeqID}_quant/lib_format_counts.json"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${samlon_bin} quant -t ${SalmonIndex} --threads ${Threads} --libType ${libType} -a ${BamDir}/${SeqID}${InputSuffix} -o ${OutputDir} ${extra_args}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
mkdir -p "$(dirname "${BamDir}")" 2>/dev/null || mkdir -p "${BamDir}"
mkdir -p "$(dirname "${OutputDir}")" 2>/dev/null || mkdir -p "${OutputDir}"

eval "$cmd"