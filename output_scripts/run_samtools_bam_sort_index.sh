#!/bin/bash
# [METADATA]
# TOOL_NAME = samtools
# VERSION = 1.10
# THREADS = 1

# Tool Info: samtools (1.10)
# Profile: bam_sort_index

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           Sequence identifier used for file naming"
    echo "  --BamDir          Directory containing the primary (unsorted) BAM file"
    echo ""
    echo "Optional Parameters:"
    echo "  --InputSuffix     Suffix of input BAM (e.g., primary, initial) (Default: primary)"
    echo "  --OutputSuffix    Suffix for output BAM (e.g., sorted, coord_sorted) (Default: sorted)"
    echo "  --samtools_bin    Path to samtools executable or binary name (Default: samtools)"
    echo "  --Threads         Number of threads for sorting and indexing (Default: 8)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
BamDir=""
InputSuffix="primary"
OutputSuffix="sorted"
samtools_bin="samtools"
Threads="8"

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --BamDir) BamDir="$2"; shift 2 ;;
        --InputSuffix) InputSuffix="$2"; shift 2 ;;
        --OutputSuffix) OutputSuffix="$2"; shift 2 ;;
        --samtools_bin) samtools_bin="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${BamDir:-}" ]]; then echo "Error: --BamDir is required"; usage; fi

# --- [Output Paths] ---
sorted_bam="${BamDir}/${SeqID}.${InputSuffix}.bam"
sorted_bai="${BamDir}/${SeqID}.${OutputSuffix}.bam.bai"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${samtools_bin} sort -@ ${Threads} -o ${sorted_bam} ${BamDir}/${SeqID}.bam && ${samtools_bin} index -b -@ ${Threads} ${sorted_bam}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
mkdir -p "$(dirname "${BamDir}")" 2>/dev/null || mkdir -p "${BamDir}"

eval "$cmd"