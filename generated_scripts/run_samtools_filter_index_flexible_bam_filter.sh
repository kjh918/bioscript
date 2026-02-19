#!/bin/bash
# [METADATA]
# TOOL_NAME = samtools_filter_index
# VERSION = 1.10
# THREADS = 1

# Tool Info: samtools_filter_index (1.10)
# Profile: flexible_bam_filter

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           Sequence identifier"
    echo "  --BamDir          Directory containing BAM files"
    echo ""
    echo "Optional Parameters:"
    echo "  --InputSuffix     Suffix of the input BAM (e.g., primary, sorted, recal) (Default: recal)"
    echo "  --OutputSuffix    Suffix for the filtered output (Default: filtered)"
    echo "  --samtools_bin    No description (Default: samtools)"
    echo "  --ln_bin          No description (Default: ln)"
    echo "  --Threads         No description (Default: 8)"
    echo "  --min_mapq        No description (Default: 20)"
    echo "  --include_flag    No description (Default: 0x2)"
    echo "  --exclude_flag    No description (Default: 0x104)"
    echo "  --expr_nm         No description (Default: [NM] < 12)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
BamDir=""
InputSuffix="recal"
OutputSuffix="filtered"
samtools_bin="samtools"
ln_bin="ln"
Threads="8"
min_mapq="20"
include_flag="0x2"
exclude_flag="0x104"
expr_nm="[NM] < 12"

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --BamDir) BamDir="$2"; shift 2 ;;
        --InputSuffix) InputSuffix="$2"; shift 2 ;;
        --OutputSuffix) OutputSuffix="$2"; shift 2 ;;
        --samtools_bin) samtools_bin="$2"; shift 2 ;;
        --ln_bin) ln_bin="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        --min_mapq) min_mapq="$2"; shift 2 ;;
        --include_flag) include_flag="$2"; shift 2 ;;
        --exclude_flag) exclude_flag="$2"; shift 2 ;;
        --expr_nm) expr_nm="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${BamDir:-}" ]]; then echo "Error: --BamDir is required"; usage; fi

# --- [Output Paths] ---
filtered_bam="${BamDir}/${SeqID}.${InputSuffix}.${OutputSuffix}.bam"
filtered_bai="${BamDir}/${SeqID}.${InputSuffix}.${OutputSuffix}.bam.bai"
analysis_ready_bam="${BamDir}/${SeqID}.analysisReady.bam"
analysis_ready_bai="${BamDir}/${SeqID}.analysisReady.bam.bai"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${samtools_bin} view -b -h -q ${min_mapq} -f ${include_flag} -F ${exclude_flag} -e '${expr_nm}' --threads ${Threads} ${BamDir}/${SeqID}.${InputSuffix}.bam > ${filtered_bam} && ${samtools_bin} index --threads ${Threads} ${filtered_bam} && ${ln_bin} -Tsf ${filtered_bam} ${analysis_ready_bam} && ${ln_bin} -Tsf ${filtered_bai} ${analysis_ready_bai}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
mkdir -p "$(dirname "${BamDir}")" 2>/dev/null || mkdir -p "${BamDir}"

eval "$cmd"