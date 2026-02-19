#!/bin/bash
# [METADATA]
# TOOL_NAME = gatk4
# VERSION = 4.4.0.0
# THREADS = 1

# Tool Info: gatk4 (4.4.0.0)
# Profile: qc_artifacts_using_singularity

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           Sequence identifier used for file naming"
    echo "  --BamDir          Directory containing the BAM file to be analyzed"
    echo "  --qcResDir        Output directory for artifact metrics"
    echo "  --ReferenceFasta  Path to the reference genome FASTA file"
    echo ""
    echo "Optional Parameters:"
    echo "  --InputSuffix     Suffix of the input BAM file (default : analysisReady / e.g., recal, sorted) (Default: analysisReady)"
    echo "  --singularity_bin Path to singularity executable (Default: singularity)"
    echo "  --gatk_bin        GATK wrapper script or binary path inside container (Default: gatk)"
    echo "  --sif             Path to GATK singularity image (.sif) (Default: /storage/images/gatk-4.4.0.0.sif)"
    echo "  --bind            Mount paths for singularity (Default: /storage,/data)"
    echo "  --Threads         Number of threads for ParallelGC (Default: 14)"
    echo "  --xmx_mb          Maximum Java heap memory (MB) (Default: 16384)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
BamDir=""
qcResDir=""
ReferenceFasta=""
InputSuffix="analysisReady"
singularity_bin="singularity"
gatk_bin="gatk"
sif="/storage/images/gatk-4.4.0.0.sif"
bind="/storage,/data"
Threads="14"
xmx_mb="16384"

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --BamDir) BamDir="$2"; shift 2 ;;
        --qcResDir) qcResDir="$2"; shift 2 ;;
        --ReferenceFasta) ReferenceFasta="$2"; shift 2 ;;
        --InputSuffix) InputSuffix="$2"; shift 2 ;;
        --singularity_bin) singularity_bin="$2"; shift 2 ;;
        --gatk_bin) gatk_bin="$2"; shift 2 ;;
        --sif) sif="$2"; shift 2 ;;
        --bind) bind="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        --xmx_mb) xmx_mb="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${BamDir:-}" ]]; then echo "Error: --BamDir is required"; usage; fi
if [[ -z "${qcResDir:-}" ]]; then echo "Error: --qcResDir is required"; usage; fi
if [[ -z "${ReferenceFasta:-}" ]]; then echo "Error: --ReferenceFasta is required"; usage; fi

# --- [Output Paths] ---
artifacts_txt="${qcResDir}/${SeqID}.artifacts.txt"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${singularity_bin} exec -B ${bind} ${sif} ${gatk_bin} CollectSequencingArtifactMetrics --java-options '-XX:ParallelGCThreads=${Threads} -Xmx${xmx_mb}m' --INPUT ${BamDir}/${SeqID}.${InputSuffix}.bam --OUTPUT ${artifacts_txt} --FILE_EXTENSION .txt --REFERENCE_SEQUENCE ${ReferenceFasta}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
mkdir -p "$(dirname "${BamDir}")" 2>/dev/null || mkdir -p "${BamDir}"
mkdir -p "$(dirname "${qcResDir}")" 2>/dev/null || mkdir -p "${qcResDir}"

eval "$cmd"