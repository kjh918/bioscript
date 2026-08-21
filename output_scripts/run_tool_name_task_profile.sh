#!/bin/bash
# [METADATA]
# TOOL_NAME = TOOL_NAME
# VERSION = VERSION
# THREADS = 8

# Tool Info: TOOL_NAME (VERSION)
# Profile: TASK_PROFILE

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           분석 고유 ID"
    echo "  --InputDir        입력 파일 디렉토리"
    echo "  --OutDir          결과 저장 디렉토리"
    echo ""
    echo "Optional Parameters:"
    echo "  --output_file     주요 결과 파일 (Default: [OutDir]/[SeqID].out)"
    echo "  --tool_bin        Tool path (Default: TOOL_COMMAND)"
    echo "  --InputSuffix     Input bam suffix (Default: .UseAnalysis.bam)"
    echo "  --OutputSuffix    Output txt suffix (Default: .txt)"
    echo "  --Threads         사용할 스레드 수 (Default: 8)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
InputDir=""
OutDir=""
output_file="[OutDir]/[SeqID].out"
tool_bin="TOOL_COMMAND"
InputSuffix=".UseAnalysis.bam"
OutputSuffix=".txt"
Threads="8"

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --InputDir) InputDir="$2"; shift 2 ;;
        --OutDir) OutDir="$2"; shift 2 ;;
        --output_file) output_file="$2"; shift 2 ;;
        --tool_bin) tool_bin="$2"; shift 2 ;;
        --InputSuffix) InputSuffix="$2"; shift 2 ;;
        --OutputSuffix) OutputSuffix="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${InputDir:-}" ]]; then echo "Error: --InputDir is required"; usage; fi
if [[ -z "${OutDir:-}" ]]; then echo "Error: --OutDir is required"; usage; fi

# --- [Output Paths] ---
output_file="${OutDir}/${SeqID}.out"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${tool_bin} INPUT_OPTIONS OUTPUT_OPTIONS"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
if [[ -n "${output_file:-}" ]]; then
  if [[ "${output_file}" == *.* ]]; then mkdir -p "$(dirname "${output_file}")"; else mkdir -p "${output_file}"; fi
fi
if [[ -n "${OutDir:-}" ]]; then
  if [[ "${OutDir}" == *.* ]]; then mkdir -p "$(dirname "${OutDir}")"; else mkdir -p "${OutDir}"; fi
fi
if [[ -n "${InputDir:-}" ]]; then
  if [[ "${InputDir}" == *.* ]]; then mkdir -p "$(dirname "${InputDir}")"; else mkdir -p "${InputDir}"; fi
fi

eval "$cmd"