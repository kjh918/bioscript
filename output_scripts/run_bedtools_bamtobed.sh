#!/bin/bash
# [METADATA]
# TOOL_NAME = bedtools
# VERSION = 2.27.1
# THREADS = 8

# Tool Info: bedtools (2.27.1)
# Profile: bamtobed

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           Sequence identifier used for file naming"
    echo "  --BamDir          Directory containing the input BAM file"
    echo "  --BamToBedDir     Target directory for final delivery of BED files"
    echo ""
    echo "Optional Parameters:"
    echo "  --bed_gz          Compressed BED file in the working directory (Default: [BamDir]/[SeqID].[InputSuffix].bed.gz)"
    echo "  --final_bed       Final BED file delivered to the target directory (Default: [BamToBedDir]/[SeqID].bed.gz)"
    echo "  --InputSuffix     Suffix of the input BAM (e.g., analysisReady, recal, sorted) (Default: analysisReady)"
    echo "  --bedtools_bin    Path to bedtools binary (Default: bedtools)"
    echo "  --bgzip_bin       Path to bgzip binary (Default: bgzip)"
    echo "  --ln_bin          Path to ln binary (Default: ln)"
    echo "  --cp_bin          Path to cp binary (Default: cp)"
    echo "  --Threads         Number of threads for bgzip compression (Default: 8)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
BamDir=""
BamToBedDir=""
bed_gz="[BamDir]/[SeqID].[InputSuffix].bed.gz"
final_bed="[BamToBedDir]/[SeqID].bed.gz"
InputSuffix="analysisReady"
bedtools_bin="bedtools"
bgzip_bin="bgzip"
ln_bin="ln"
cp_bin="cp"
Threads="8"

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --BamDir) BamDir="$2"; shift 2 ;;
        --BamToBedDir) BamToBedDir="$2"; shift 2 ;;
        --bed_gz) bed_gz="$2"; shift 2 ;;
        --final_bed) final_bed="$2"; shift 2 ;;
        --InputSuffix) InputSuffix="$2"; shift 2 ;;
        --bedtools_bin) bedtools_bin="$2"; shift 2 ;;
        --bgzip_bin) bgzip_bin="$2"; shift 2 ;;
        --ln_bin) ln_bin="$2"; shift 2 ;;
        --cp_bin) cp_bin="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${BamDir:-}" ]]; then echo "Error: --BamDir is required"; usage; fi
if [[ -z "${BamToBedDir:-}" ]]; then echo "Error: --BamToBedDir is required"; usage; fi

# --- [Output Paths] ---
bed_gz="${BamDir}/${SeqID}.${InputSuffix}.bed.gz"
final_bed="${BamToBedDir}/${SeqID}.bed.gz"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${bedtools_bin} bamtobed -i ${BamDir}/${SeqID}.${InputSuffix}.bam > ${BamDir}/${SeqID}.${InputSuffix}.bed && ${bgzip_bin} -f -@ ${Threads} ${BamDir}/${SeqID}.${InputSuffix}.bed && ${ln_bin} -Tsf ${bed_gz} ${BamDir}/${SeqID}.bed.gz && ${cp_bin} ${bed_gz} ${final_bed}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
if [[ -n "${bed_gz:-}" ]]; then
  if [[ "${bed_gz}" == *.* ]]; then mkdir -p "$(dirname "${bed_gz}")"; else mkdir -p "${bed_gz}"; fi
fi
if [[ -n "${final_bed:-}" ]]; then
  if [[ "${final_bed}" == *.* ]]; then mkdir -p "$(dirname "${final_bed}")"; else mkdir -p "${final_bed}"; fi
fi
if [[ -n "${BamToBedDir:-}" ]]; then
  if [[ "${BamToBedDir}" == *.* ]]; then mkdir -p "$(dirname "${BamToBedDir}")"; else mkdir -p "${BamToBedDir}"; fi
fi
if [[ -n "${BamDir:-}" ]]; then
  if [[ "${BamDir}" == *.* ]]; then mkdir -p "$(dirname "${BamDir}")"; else mkdir -p "${BamDir}"; fi
fi

eval "$cmd"