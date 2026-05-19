#!/bin/bash
# [METADATA]
# TOOL_NAME = cbnipt_cnv_call_cli
# VERSION = 1.0.0
# THREADS = 1

# Tool Info: cbnipt_cnv_call_cli (1.0.0)
# Profile: cnv_calling

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           The parsed identifier extracted from the BAM filename."
    echo "  --BamPath         Absolute path to the specific analysis-ready BAM file."
    echo "  --ReferenceFasta  Absolute path to the hg38 reference genome FASTA file."
    echo "  --OutDir          Absolute path to the specific run output directory."
    echo ""
    echo "Optional Parameters:"
    echo "  --PythonBin       Absolute path to the target Python interpreter environment. (Default: /storage/home/jhkim/Apps/Python-3.11.13/python)"
    echo "  --CliScript       Absolute path to the target cli.py script tool. (Default: /storage/home/jhkim/scripts/bioscript/manual/cbnipt_manual_cnv_call/scripts/cli.py)"
    echo "  --Threads         The parsed identifier extracted from the BAM filename. (Default: 4)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
BamPath=""
ReferenceFasta=""
OutDir=""
PythonBin="/storage/home/jhkim/Apps/Python-3.11.13/python"
CliScript="/storage/home/jhkim/scripts/bioscript/manual/cbnipt_manual_cnv_call/scripts/cli.py"
Threads="4"

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --BamPath) BamPath="$2"; shift 2 ;;
        --ReferenceFasta) ReferenceFasta="$2"; shift 2 ;;
        --OutDir) OutDir="$2"; shift 2 ;;
        --PythonBin) PythonBin="$2"; shift 2 ;;
        --CliScript) CliScript="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${BamPath:-}" ]]; then echo "Error: --BamPath is required"; usage; fi
if [[ -z "${ReferenceFasta:-}" ]]; then echo "Error: --ReferenceFasta is required"; usage; fi
if [[ -z "${OutDir:-}" ]]; then echo "Error: --OutDir is required"; usage; fi

# --- [Output Paths] ---
OutDir=""

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${PythonBin} ${CliScript} --SeqID ${SeqID} --BamPath ${BamPath} --ReferenceFasta ${ReferenceFasta} --OutDir ${OutDir} --Threads ${Threads}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
mkdir -p "$(dirname "${OutDir}")" 2>/dev/null || mkdir -p "${OutDir}"

eval "$cmd"