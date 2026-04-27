#!/bin/bash
# [METADATA]
# TOOL_NAME = gatk4_learn_read_orientation_model
# VERSION = 4.4.0.0
# THREADS = 1

# Tool Info: gatk4_learn_read_orientation_model (4.4.0.0)
# Profile: orientation_bias_learning

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           No description"
    echo "  --qcResDir        f1r2 파일이 있고 모델을 저장할 디렉토리"
    echo ""
    echo "Optional Parameters:"
    echo "  --InputSuffix     Input f1r2 suffix (with dot if exists) (Default: )"
    echo "  --OutputSuffix    Output model suffix (Default: )"
    echo "  --singularity_bin No description (Default: singularity)"
    echo "  --gatk4_sif       No description (Default: /storage/images/gatk-4.4.0.0.sif)"
    echo "  --bind            No description (Default: /storage,/data)"
    echo "  --xmx_mb          8GB memory from original script (Default: 8192)"
    echo "  --Threads         ParallelGCThreads (Default: 14)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
qcResDir=""
InputSuffix=""
OutputSuffix=""
singularity_bin="singularity"
gatk4_sif="/storage/images/gatk-4.4.0.0.sif"
bind="/storage,/data"
xmx_mb="8192"
Threads="14"

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --qcResDir) qcResDir="$2"; shift 2 ;;
        --InputSuffix) InputSuffix="$2"; shift 2 ;;
        --OutputSuffix) OutputSuffix="$2"; shift 2 ;;
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
if [[ -z "${qcResDir:-}" ]]; then echo "Error: --qcResDir is required"; usage; fi

# --- [Output Paths] ---
read_orientation_model="${qcResDir}/${SeqID}.${OutputSuffix}.read-orientation-model.tar.gz"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${singularity_bin} exec -B ${bind} ${gatk4_sif} gatk LearnReadOrientationModel --java-options '-XX:ParallelGCThreads=${Threads} -Xmx${xmx_mb}m' -I ${qcResDir}/${SeqID}.${InputSuffix}f1r2.tar.gz -O ${read_orientation_model}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
mkdir -p "$(dirname "${qcResDir}")" 2>/dev/null || mkdir -p "${qcResDir}"

eval "$cmd"