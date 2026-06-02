#!/bin/bash
# [METADATA]
# TOOL_NAME = TOOL_NAME
# VERSION = TASK_VERSION
# THREADS = 8

# Tool Info: TOOL_NAME (TASK_VERSION)
# Profile: TOOL_FUNCTION

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           Unique sample identifier used to generate the output file names"
    echo "  --InputFileDir    Path of Input sample data."
    echo "  --qcResDir        Path of QC Results."
    echo "  --OutputDir       Path of Output files."
    echo "  --OutputSummary   Path of Output files."
    echo ""
    echo "Optional Parameters:"
    echo "  --InputSuffix     File name input-suffix (e.g., '.bam') (Default: )"
    echo "  --OutputSuffix    File name output-suffix (e.g., '.summary.txt') (Default: )"
    echo "  --singularity_bin No description (Default: singularity)"
    echo "  --bind            No description (Default: /storage,/data)"
    echo "  --xmx_mb          16GB memory from original script (Default: 16384)"
    echo "  --Threads         Threads (default : 8) (Default: 8)"
    echo "  --extra_args      Extra arguments (Default: )"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
InputFileDir=""
qcResDir=""
OutputDir=""
OutputSummary=""
InputSuffix=""
OutputSuffix=""
singularity_bin="singularity"
bind="/storage,/data"
xmx_mb="16384"
Threads="8"
extra_args=""

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --InputFileDir) InputFileDir="$2"; shift 2 ;;
        --qcResDir) qcResDir="$2"; shift 2 ;;
        --OutputDir) OutputDir="$2"; shift 2 ;;
        --OutputSummary) OutputSummary="$2"; shift 2 ;;
        --InputSuffix) InputSuffix="$2"; shift 2 ;;
        --OutputSuffix) OutputSuffix="$2"; shift 2 ;;
        --singularity_bin) singularity_bin="$2"; shift 2 ;;
        --bind) bind="$2"; shift 2 ;;
        --xmx_mb) xmx_mb="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        --extra_args) extra_args="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${InputFileDir:-}" ]]; then echo "Error: --InputFileDir is required"; usage; fi
if [[ -z "${qcResDir:-}" ]]; then echo "Error: --qcResDir is required"; usage; fi
if [[ -z "${OutputDir:-}" ]]; then echo "Error: --OutputDir is required"; usage; fi
if [[ -z "${OutputSummary:-}" ]]; then echo "Error: --OutputSummary is required"; usage; fi

# --- [Output Paths] ---
OutputSummary=""

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${singularity_bin} exec -B ${bind}  ${extra_args}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
if [[ -n "${InputFileDir:-}" ]]; then
  if [[ "${InputFileDir}" == *.* ]]; then mkdir -p "$(dirname "${InputFileDir}")"; else mkdir -p "${InputFileDir}"; fi
fi
if [[ -n "${OutputDir:-}" ]]; then
  if [[ "${OutputDir}" == *.* ]]; then mkdir -p "$(dirname "${OutputDir}")"; else mkdir -p "${OutputDir}"; fi
fi
if [[ -n "${OutputSummary:-}" ]]; then
  if [[ "${OutputSummary}" == *.* ]]; then mkdir -p "$(dirname "${OutputSummary}")"; else mkdir -p "${OutputSummary}"; fi
fi
if [[ -n "${qcResDir:-}" ]]; then
  if [[ "${qcResDir}" == *.* ]]; then mkdir -p "$(dirname "${qcResDir}")"; else mkdir -p "${qcResDir}"; fi
fi

eval "$cmd"