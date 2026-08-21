#!/bin/bash
# [METADATA]
# TOOL_NAME = samtools_flagstats
# VERSION = 1.10
# THREADS = 1

# Tool Info: samtools_flagstats (1.10)
# Profile: summary_flagstats

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           Sequence identifier"
    echo "  --BamDir          Directory containing BAM files"
    echo ""
    echo "Optional Parameters:"
    echo "  --flag_stats_txt  Flag Stats from BAM file (Default: [QcDir]/[SeqID][InputSuffix].[OutputSuffix].txt)"
    echo "  --InputSuffix     Suffix of the input BAM (e.g., primary, sorted, recal) (Default: .sorted)"
    echo "  --OutputSuffix    Suffix for the filtered output (Default: filtered)"
    echo "  --samtools_bin    No description (Default: samtools)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
BamDir=""
flag_stats_txt="[QcDir]/[SeqID][InputSuffix].[OutputSuffix].txt"
InputSuffix=".sorted"
OutputSuffix="filtered"
samtools_bin="samtools"

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --BamDir) BamDir="$2"; shift 2 ;;
        --flag_stats_txt) flag_stats_txt="$2"; shift 2 ;;
        --InputSuffix) InputSuffix="$2"; shift 2 ;;
        --OutputSuffix) OutputSuffix="$2"; shift 2 ;;
        --samtools_bin) samtools_bin="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${BamDir:-}" ]]; then echo "Error: --BamDir is required"; usage; fi

# --- [Output Paths] ---
flag_stats_txt="[QcDir]/${SeqID}${InputSuffix}.${OutputSuffix}.txt"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${samtools_bin} flagstats ${BamDir}/[SeqId]${InputSuffix}.bam > ${flag_stats_txt}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
if [[ -n "${flag_stats_txt:-}" ]]; then
  if [[ "${flag_stats_txt}" == *.* ]]; then mkdir -p "$(dirname "${flag_stats_txt}")"; else mkdir -p "${flag_stats_txt}"; fi
fi
if [[ -n "${BamDir:-}" ]]; then
  if [[ "${BamDir}" == *.* ]]; then mkdir -p "$(dirname "${BamDir}")"; else mkdir -p "${BamDir}"; fi
fi

eval "$cmd"