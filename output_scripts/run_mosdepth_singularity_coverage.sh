#!/bin/bash
# [METADATA]
# TOOL_NAME = mosdepth
# VERSION = 0.3.6
# THREADS = 1

# Tool Info: mosdepth (0.3.6)
# Profile: singularity_coverage

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           Sequence identifier used for file tracking"
    echo "  --BamDir          Directory containing the BAM file to be analyzed"
    echo "  --qcResDir        Output directory for coverage results"
    echo ""
    echo "Optional Parameters:"
    echo "  --InputSuffix     Suffix of the input BAM (e.g., analysisReady, recal, sorted) (Default: analysisReady)"
    echo "  --singularity_bin Path to singularity executable (Default: singularity)"
    echo "  --mosdepth_bin    Path to mosdepth binary inside container (Default: /opt/mosdepth)"
    echo "  --sif             mosdepth singularity image file (Default: /storage/images/mosdepth-0.3.6.sif)"
    echo "  --bind            Mount paths for singularity (Default: /storage,/data)"
    echo "  --Threads         Number of threads for decompression and processing (Default: 8)"
    echo "  --bin             Window size for quantized output (used with --by) (Default: 100000)"
    echo "  --mapq            Minimum mapping quality for a read to be counted (Default: 20)"
    echo "  --mosdepth_args   Additional mosdepth arguments (Default: --no-per-base --fast-mode)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
BamDir=""
qcResDir=""
InputSuffix="analysisReady"
singularity_bin="singularity"
mosdepth_bin="/opt/mosdepth"
sif="/storage/images/mosdepth-0.3.6.sif"
bind="/storage,/data"
Threads="8"
bin="100000"
mapq="20"
mosdepth_args="--no-per-base --fast-mode"

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --BamDir) BamDir="$2"; shift 2 ;;
        --qcResDir) qcResDir="$2"; shift 2 ;;
        --InputSuffix) InputSuffix="$2"; shift 2 ;;
        --singularity_bin) singularity_bin="$2"; shift 2 ;;
        --mosdepth_bin) mosdepth_bin="$2"; shift 2 ;;
        --sif) sif="$2"; shift 2 ;;
        --bind) bind="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        --bin) bin="$2"; shift 2 ;;
        --mapq) mapq="$2"; shift 2 ;;
        --mosdepth_args) mosdepth_args="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${BamDir:-}" ]]; then echo "Error: --BamDir is required"; usage; fi
if [[ -z "${qcResDir:-}" ]]; then echo "Error: --qcResDir is required"; usage; fi

# --- [Output Paths] ---
prefix="${qcResDir}/${SeqID}.${InputSuffix}"
regions_bed_gz="${prefix}.regions.bed.gz"
regions_bed_gz_csi="${prefix}.regions.bed.gz.csi"
summary_txt="${prefix}.mosdepth.summary.txt"
global_dist_txt="${prefix}.global.dist.txt"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${singularity_bin} exec -B ${bind} ${sif} ${mosdepth_bin} --threads ${Threads} ${mosdepth_args} --by ${bin} --mapq ${mapq} ${prefix} ${BamDir}/${SeqID}.${InputSuffix}.bam"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
mkdir -p "$(dirname "${BamDir}")" 2>/dev/null || mkdir -p "${BamDir}"
mkdir -p "$(dirname "${qcResDir}")" 2>/dev/null || mkdir -p "${qcResDir}"

eval "$cmd"