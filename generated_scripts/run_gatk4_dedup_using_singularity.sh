#!/bin/bash
# [METADATA]
# TOOL_NAME = gatk4
# VERSION = 4.4.0.0
# THREADS = 1

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
    echo "  --other_md_args   Other additional GATK MarkDuplicates arguments (Default: )"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
BamDir=""
qcResDir=""
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
other_md_args=""

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --BamDir) BamDir="$2"; shift 2 ;;
        --qcResDir) qcResDir="$2"; shift 2 ;;
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
cmd="${singularity_bin} exec -B ${bind} ${sif} ${gatk_bin} MarkDuplicates --java-options '-XX:ParallelGCThreads=${Threads} -Xmx${xmx_mb}m' --INPUT ${BamDir}/${SeqID}.${InputSuffix}.bam --OUTPUT ${dedup_bam} --METRICS_FILE ${metrics_txt} --CREATE_INDEX ${create_index} --REMOVE_SEQUENCING_DUPLICATES ${remove_seq_dups} ${other_md_args}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
mkdir -p "$(dirname "${BamDir}")" 2>/dev/null || mkdir -p "${BamDir}"
mkdir -p "$(dirname "${qcResDir}")" 2>/dev/null || mkdir -p "${qcResDir}"

eval "$cmd"