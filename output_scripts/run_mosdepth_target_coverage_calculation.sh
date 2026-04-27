#!/bin/bash
# [METADATA]
# TOOL_NAME = mosdepth
# VERSION = 0.3.6
# THREADS = 1

# Tool Info: mosdepth (0.3.6)
# Profile: target_coverage_calculation

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           Sequence identifier used for file tracking"
    echo "  --BamDir          Directory containing the BAM file to be analyzed"
    echo "  --qcResDir        Directory for quality control results"
    echo "  --TargetBed       Target intervals (BED) for coverage calculation"
    echo ""
    echo "Optional Parameters:"
    echo "  --InputSuffix     No description (Default: analysisReady)"
    echo "  --singularity_bin No description (Default: singularity)"
    echo "  --mosdepth_sif    No description (Default: /storage/images/mosdepth-0.3.6.sif)"
    echo "  --mosdepth_bin    No description (Default: /opt/mosdepth)"
    echo "  --bind            No description (Default: /storage,/data)"
    echo "  --Threads         Number of threads for decompression (Default: 4)"
    echo "  --MapQ            Mapping quality threshold (default: 20) (Default: 20)"
    echo "  --extra_args      Additional mosdepth flags (e.g., --no-per-base) (Default: )"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
BamDir=""
qcResDir=""
TargetBed=""
InputSuffix="analysisReady"
singularity_bin="singularity"
mosdepth_sif="/storage/images/mosdepth-0.3.6.sif"
mosdepth_bin="/opt/mosdepth"
bind="/storage,/data"
Threads="4"
MapQ="20"
extra_args=""

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --BamDir) BamDir="$2"; shift 2 ;;
        --qcResDir) qcResDir="$2"; shift 2 ;;
        --TargetBed) TargetBed="$2"; shift 2 ;;
        --InputSuffix) InputSuffix="$2"; shift 2 ;;
        --singularity_bin) singularity_bin="$2"; shift 2 ;;
        --mosdepth_sif) mosdepth_sif="$2"; shift 2 ;;
        --mosdepth_bin) mosdepth_bin="$2"; shift 2 ;;
        --bind) bind="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        --MapQ) MapQ="$2"; shift 2 ;;
        --extra_args) extra_args="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${BamDir:-}" ]]; then echo "Error: --BamDir is required"; usage; fi
if [[ -z "${qcResDir:-}" ]]; then echo "Error: --qcResDir is required"; usage; fi
if [[ -z "${TargetBed:-}" ]]; then echo "Error: --TargetBed is required"; usage; fi

# --- [Output Paths] ---
global_dist="${qcResDir}/${SeqID}.mosdepth.global.dist.txt"
region_dist="${qcResDir}/${SeqID}.mosdepth.region.dist.txt"
summary_txt="${qcResDir}/${SeqID}.mosdepth.summary.txt"
regions_bed="${qcResDir}/${SeqID}.regions.bed.gz"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${singularity_bin} exec -B ${bind} ${mosdepth_sif} ${mosdepth_bin} --threads ${Threads} --by ${TargetBed} --mapq ${MapQ} ${extra_args} ${qcResDir}/${SeqID} ${BamDir}/${SeqID}.${InputSuffix}.bam"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
mkdir -p "$(dirname "${BamDir}")" 2>/dev/null || mkdir -p "${BamDir}"
mkdir -p "$(dirname "${qcResDir}")" 2>/dev/null || mkdir -p "${qcResDir}"

eval "$cmd"