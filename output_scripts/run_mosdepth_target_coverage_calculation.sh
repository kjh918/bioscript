#!/bin/bash
# [METADATA]
# TOOL_NAME = mosdepth
# VERSION = 0.3.6
# THREADS = 4

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
    echo "  --global_dist     No description (Default: [qcResDir]/[SeqID].mosdepth.global.dist.txt)"
    echo "  --region_dist     No description (Default: [qcResDir]/[SeqID].mosdepth.region.dist.txt)"
    echo "  --summary_txt     No description (Default: [qcResDir]/[SeqID].mosdepth.summary.txt)"
    echo "  --regions_bed     No description (Default: [qcResDir]/[SeqID].regions.bed.gz)"
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
global_dist="[qcResDir]/[SeqID].mosdepth.global.dist.txt"
region_dist="[qcResDir]/[SeqID].mosdepth.region.dist.txt"
summary_txt="[qcResDir]/[SeqID].mosdepth.summary.txt"
regions_bed="[qcResDir]/[SeqID].regions.bed.gz"
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
        --global_dist) global_dist="$2"; shift 2 ;;
        --region_dist) region_dist="$2"; shift 2 ;;
        --summary_txt) summary_txt="$2"; shift 2 ;;
        --regions_bed) regions_bed="$2"; shift 2 ;;
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
if [[ -n "${qcResDir:-}" ]]; then
  if [[ "${qcResDir}" == *.* ]]; then mkdir -p "$(dirname "${qcResDir}")"; else mkdir -p "${qcResDir}"; fi
fi
if [[ -n "${region_dist:-}" ]]; then
  if [[ "${region_dist}" == *.* ]]; then mkdir -p "$(dirname "${region_dist}")"; else mkdir -p "${region_dist}"; fi
fi
if [[ -n "${summary_txt:-}" ]]; then
  if [[ "${summary_txt}" == *.* ]]; then mkdir -p "$(dirname "${summary_txt}")"; else mkdir -p "${summary_txt}"; fi
fi
if [[ -n "${BamDir:-}" ]]; then
  if [[ "${BamDir}" == *.* ]]; then mkdir -p "$(dirname "${BamDir}")"; else mkdir -p "${BamDir}"; fi
fi
if [[ -n "${regions_bed:-}" ]]; then
  if [[ "${regions_bed}" == *.* ]]; then mkdir -p "$(dirname "${regions_bed}")"; else mkdir -p "${regions_bed}"; fi
fi
if [[ -n "${global_dist:-}" ]]; then
  if [[ "${global_dist}" == *.* ]]; then mkdir -p "$(dirname "${global_dist}")"; else mkdir -p "${global_dist}"; fi
fi

eval "$cmd"