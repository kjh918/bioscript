#!/bin/bash
# [METADATA]
# TOOL_NAME = cbNIPT_CallCNV
# VERSION = 1.0.0
# THREADS = 4

# Tool Info: cbNIPT_CallCNV (1.0.0)
# Profile: Calculate CNV and variant profile

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           분석 고유 ID"
    echo "  --BamPath         분석할 단일 BAM 파일 경로"
    echo "  --AnnotatedBins   make-bins로 생성한 사전 계산 빈 파일(.bed.gz)"
    echo "  --OutDir          결과 저장 디렉토리"
    echo "  --script_path     Path to the CNV analyzer script"
    echo ""
    echo "Optional Parameters:"
    echo "  --MinMapQ         최소 Mapping Quality (Default: 20)"
    echo "  --LowessFrac      GC 보정 LOWESS 비율 (Default: 0.2)"
    echo "  --BaselinePloidy  기준 Ploidy (Default: 2)"
    echo "  --SmoothWindow    Rolling median window size (Default: 5)"
    echo "  --MinDepth        최소 Depth (Default: 2.0)"
    echo "  --MinCoverage     최소 Coverage 비율 (Default: 0.5)"
    echo "  --MaskLowerQ      하위 Coverage 필터 (Default: 0.01)"
    echo "  --MaskUpperQ      상위 Coverage 필터 (Default: 0.99)"
    echo "  --SegPenalty      Segmentation 민감도 (Default: 10.0)"
    echo "  --Threads         사용할 스레드 수 (Default: 4)"
    echo "  --python_bin      Python executable (Default: python)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
BamPath=""
AnnotatedBins=""
OutDir=""
MinMapQ="20"
LowessFrac="0.2"
BaselinePloidy="2"
SmoothWindow="5"
MinDepth="2.0"
MinCoverage="0.5"
MaskLowerQ="0.01"
MaskUpperQ="0.99"
SegPenalty="10.0"
Threads="4"
python_bin="python"
script_path=""

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --BamPath) BamPath="$2"; shift 2 ;;
        --AnnotatedBins) AnnotatedBins="$2"; shift 2 ;;
        --OutDir) OutDir="$2"; shift 2 ;;
        --MinMapQ) MinMapQ="$2"; shift 2 ;;
        --LowessFrac) LowessFrac="$2"; shift 2 ;;
        --BaselinePloidy) BaselinePloidy="$2"; shift 2 ;;
        --SmoothWindow) SmoothWindow="$2"; shift 2 ;;
        --MinDepth) MinDepth="$2"; shift 2 ;;
        --MinCoverage) MinCoverage="$2"; shift 2 ;;
        --MaskLowerQ) MaskLowerQ="$2"; shift 2 ;;
        --MaskUpperQ) MaskUpperQ="$2"; shift 2 ;;
        --SegPenalty) SegPenalty="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        --python_bin) python_bin="$2"; shift 2 ;;
        --script_path) script_path="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${BamPath:-}" ]]; then echo "Error: --BamPath is required"; usage; fi
if [[ -z "${AnnotatedBins:-}" ]]; then echo "Error: --AnnotatedBins is required"; usage; fi
if [[ -z "${OutDir:-}" ]]; then echo "Error: --OutDir is required"; usage; fi
if [[ -z "${script_path:-}" ]]; then echo "Error: --script_path is required"; usage; fi

# --- [Output Paths] ---
OutDir=""

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${python_bin} ${script_path} call-cnv --SeqID ${SeqID} --BamPath ${BamPath} --AnnotatedBins ${AnnotatedBins} --OutDir ${OutDir} --MinMapQ ${MinMapQ} --LowessFrac ${LowessFrac} --BaselinePloidy ${BaselinePloidy} --SmoothWindow ${SmoothWindow} --MinDepth ${MinDepth} --MinCoverage ${MinCoverage} --MaskLowerQ ${MaskLowerQ} --MaskUpperQ ${MaskUpperQ} --SegPenalty ${SegPenalty} --Threads ${Threads}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
if [[ -n "${OutDir:-}" ]]; then
  if [[ "${OutDir}" == *.* ]]; then mkdir -p "$(dirname "${OutDir}")"; else mkdir -p "${OutDir}"; fi
fi

eval "$cmd"