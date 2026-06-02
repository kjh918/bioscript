#!/bin/bash
# [METADATA]
# TOOL_NAME = gatk4
# VERSION = 4.4.0.0
# THREADS = 14

# Tool Info: gatk4 (4.4.0.0)
# Profile: dedup_using_singularity

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           Sequence identifier used for file naming"
    echo "  --BamDir          Directory containing the input BAM file"
    echo "  --qcResDir        Directory for duplicate metrics output"
    echo ""
    echo "Optional Parameters:"
    echo "  --dedup_bam       BAM file with marked/removed duplicates (Default: [BamDir]/[SeqID].[InputSuffix].[OutputSuffix].bam)"
    echo "  --dedup_bai       Index file for the deduplicated BAM (Default: [BamDir]/[SeqID].[InputSuffix].[OutputSuffix].bam.bai)"
    echo "  --metrics_txt     Text file containing duplication metrics (Default: [qcResDir]/[SeqID].[InputSuffix].[OutputSuffix].metrics.txt)"
    echo "  --InputSuffix     Suffix of the input BAM (e.g., sorted, primary) (Default: sorted)"
    echo "  --OutputSuffix    Suffix for the output file (e.g., dedup, md) (Default: dedup)"
    echo "  --singularity_bin Path to singularity executable (Default: singularity)"
    echo "  --gatk_bin        GATK binary path inside container (Default: gatk)"
    echo "  --sif             GATK singularity image (.sif) (Default: /storage/images/gatk-4.4.0.0.sif)"
    echo "  --bind            Mount paths for singularity (Default: /storage,/data)"
    echo "  --Threads         Number of threads for ParallelGC (Default: 14)"
    echo "  --xmx_mb          Maximum Java heap memory (MB) (Default: 16384)"
    echo "  --remove_seq_dups Whether to remove sequencing duplicates (Default: true)"
    echo "  --create_index    Whether to create a BAM index for the output (Default: true)"
    echo "  --optical_duplicate_pixel_distance Pixel distance for optical duplicate detection (Default: 2500)"
    echo "  --other_md_args   Other additional GATK MarkDuplicates arguments (Default: )"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
BamDir=""
qcResDir=""
dedup_bam="[BamDir]/[SeqID].[InputSuffix].[OutputSuffix].bam"
dedup_bai="[BamDir]/[SeqID].[InputSuffix].[OutputSuffix].bam.bai"
metrics_txt="[qcResDir]/[SeqID].[InputSuffix].[OutputSuffix].metrics.txt"
InputSuffix="sorted"
OutputSuffix="dedup"
singularity_bin="singularity"
gatk_bin="gatk"
sif="/storage/images/gatk-4.4.0.0.sif"
bind="/storage,/data"
Threads="14"
xmx_mb="16384"
remove_seq_dups="true"
create_index="true"
optical_duplicate_pixel_distance="2500"
other_md_args=""

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --BamDir) BamDir="$2"; shift 2 ;;
        --qcResDir) qcResDir="$2"; shift 2 ;;
        --dedup_bam) dedup_bam="$2"; shift 2 ;;
        --dedup_bai) dedup_bai="$2"; shift 2 ;;
        --metrics_txt) metrics_txt="$2"; shift 2 ;;
        --InputSuffix) InputSuffix="$2"; shift 2 ;;
        --OutputSuffix) OutputSuffix="$2"; shift 2 ;;
        --singularity_bin) singularity_bin="$2"; shift 2 ;;
        --gatk_bin) gatk_bin="$2"; shift 2 ;;
        --sif) sif="$2"; shift 2 ;;
        --bind) bind="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        --xmx_mb) xmx_mb="$2"; shift 2 ;;
        --remove_seq_dups) remove_seq_dups="$2"; shift 2 ;;
        --create_index) create_index="$2"; shift 2 ;;
        --optical_duplicate_pixel_distance) optical_duplicate_pixel_distance="$2"; shift 2 ;;
        --other_md_args) other_md_args="$2"; shift 2 ;;
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
dedup_bam="${BamDir}/${SeqID}.${InputSuffix}.${OutputSuffix}.bam"
dedup_bai="${BamDir}/${SeqID}.${InputSuffix}.${OutputSuffix}.bam.bai"
metrics_txt="${qcResDir}/${SeqID}.${InputSuffix}.${OutputSuffix}.metrics.txt"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${singularity_bin} exec -B ${bind} ${sif} ${gatk_bin} MarkDuplicates --java-options '-XX:ParallelGCThreads=${Threads} -Xmx${xmx_mb}m' --INPUT ${BamDir}/${SeqID}.${InputSuffix}.bam --OUTPUT ${dedup_bam} --METRICS_FILE ${metrics_txt} --CREATE_INDEX ${create_index} --REMOVE_SEQUENCING_DUPLICATES ${remove_seq_dups} --OPTICAL_DUPLICATE_PIXEL_DISTANCE ${optical_duplicate_pixel_distance} ${other_md_args}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
if [[ -n "${dedup_bai:-}" ]]; then
  if [[ "${dedup_bai}" == *.* ]]; then mkdir -p "$(dirname "${dedup_bai}")"; else mkdir -p "${dedup_bai}"; fi
fi
if [[ -n "${dedup_bam:-}" ]]; then
  if [[ "${dedup_bam}" == *.* ]]; then mkdir -p "$(dirname "${dedup_bam}")"; else mkdir -p "${dedup_bam}"; fi
fi
if [[ -n "${BamDir:-}" ]]; then
  if [[ "${BamDir}" == *.* ]]; then mkdir -p "$(dirname "${BamDir}")"; else mkdir -p "${BamDir}"; fi
fi
if [[ -n "${metrics_txt:-}" ]]; then
  if [[ "${metrics_txt}" == *.* ]]; then mkdir -p "$(dirname "${metrics_txt}")"; else mkdir -p "${metrics_txt}"; fi
fi
if [[ -n "${qcResDir:-}" ]]; then
  if [[ "${qcResDir}" == *.* ]]; then mkdir -p "$(dirname "${qcResDir}")"; else mkdir -p "${qcResDir}"; fi
fi

eval "$cmd"