#!/bin/bash
# [METADATA]
# TOOL_NAME = gatk4_hs_metrics
# VERSION = 4.4.0.0
# THREADS = 1

# Tool Info: gatk4_hs_metrics (4.4.0.0)
# Profile: wes_target_efficiency

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           No description"
    echo "  --BamDir          No description"
    echo "  --qcResDir        No description"
    echo "  --ReferenceFasta  No description"
    echo "  --BaitIntervals   Target capture kit bait intervals (Picard-style)"
    echo "  --TargetIntervals Target capture kit target intervals (Picard-style)"
    echo "  --InputSuffix     Target BAM suffix"
    echo ""
    echo "Optional Parameters:"
    echo "  --singularity_bin No description (Default: singularity)"
    echo "  --gatk4_sif       No description (Default: /storage/images/gatk-4.4.0.0.sif)"
    echo "  --bind            No description (Default: /storage,/data)"
    echo "  --xmx_mb          No description (Default: 16384)"
    echo "  --Threads         No description (Default: 14)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
BamDir=""
qcResDir=""
ReferenceFasta=""
BaitIntervals=""
TargetIntervals=""
InputSuffix="analysisReady"
singularity_bin="singularity"
gatk4_sif="/storage/images/gatk-4.4.0.0.sif"
bind="/storage,/data"
xmx_mb="16384"
Threads="14"

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --BamDir) BamDir="$2"; shift 2 ;;
        --qcResDir) qcResDir="$2"; shift 2 ;;
        --ReferenceFasta) ReferenceFasta="$2"; shift 2 ;;
        --BaitIntervals) BaitIntervals="$2"; shift 2 ;;
        --TargetIntervals) TargetIntervals="$2"; shift 2 ;;
        --InputSuffix) InputSuffix="$2"; shift 2 ;;
        --singularity_bin) singularity_bin="$2"; shift 2 ;;
        --gatk4_sif) gatk4_sif="$2"; shift 2 ;;
        --bind) bind="$2"; shift 2 ;;
        --xmx_mb) xmx_mb="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
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
if [[ -z "${BaitIntervals:-}" ]]; then echo "Error: --BaitIntervals is required"; usage; fi
if [[ -z "${TargetIntervals:-}" ]]; then echo "Error: --TargetIntervals is required"; usage; fi
if [[ -z "${InputSuffix:-}" ]]; then echo "Error: --InputSuffix is required"; usage; fi

# --- [Output Paths] ---
hs_metrics="${qcResDir}/${SeqID}.${InputSuffix}.HS.metrics.txt"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${singularity_bin} exec -B ${bind} ${gatk4_sif} gatk CollectHsMetrics --java-options '-XX:ParallelGCThreads=${Threads} -Xmx${xmx_mb}m' --INPUT ${BamDir}/${SeqID}.${InputSuffix}.bam --OUTPUT ${hs_metrics} --REFERENCE_SEQUENCE ${ReferenceFasta} --BAIT_INTERVALS ${BaitIntervals} --TARGET_INTERVALS ${TargetIntervals}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
mkdir -p "$(dirname "${BamDir}")" 2>/dev/null || mkdir -p "${BamDir}"
mkdir -p "$(dirname "${qcResDir}")" 2>/dev/null || mkdir -p "${qcResDir}"

eval "$cmd"