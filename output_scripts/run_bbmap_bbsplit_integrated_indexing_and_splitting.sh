#!/bin/bash
# [METADATA]
# TOOL_NAME = bbmap_bbsplit
# VERSION = 39.81
# THREADS = 1

# Tool Info: bbmap_bbsplit (39.81)
# Profile: integrated_indexing_and_splitting

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           Sample identifier for naming outputs"
    echo "  --Read1           No description"
    echo "  --Read2           No description"
    echo "  --PrimaryRef      Path to primary reference fasta"
    echo "  --OtherRefsArgs   Formatted string of additional references"
    echo "  --OutputDir       Directory for results and statistics"
    echo "  --MemoryMB        Java heap memory (e.g., 32768)"
    echo "  --Threads         Number of CPU threads"
    echo "  --BBMapPath       Directory path containing bbsplit.sh"
    echo "  --BaseNamePattern No description"
    echo ""
    echo "Optional Parameters:"
    echo "  --ExtraArgs       No description (Default: ambiguous2=all)"
    echo "  --singularity_bin No description (Default: singularity)"
    echo "  --bbmap_sif       No description (Default: /storage/home/jhkim/Apps/bbmap_39.81.sif)"
    echo "  --bind            No description (Default: /storage,/data)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
Read1=""
Read2=""
PrimaryRef=""
OtherRefsArgs=""
OutputDir=""
MemoryMB=""
Threads=""
BBMapPath=""
BaseNamePattern=""
ExtraArgs="ambiguous2=all"
singularity_bin="singularity"
bbmap_sif="/storage/home/jhkim/Apps/bbmap_39.81.sif"
bind="/storage,/data"

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --Read1) Read1="$2"; shift 2 ;;
        --Read2) Read2="$2"; shift 2 ;;
        --PrimaryRef) PrimaryRef="$2"; shift 2 ;;
        --OtherRefsArgs) OtherRefsArgs="$2"; shift 2 ;;
        --OutputDir) OutputDir="$2"; shift 2 ;;
        --MemoryMB) MemoryMB="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        --BBMapPath) BBMapPath="$2"; shift 2 ;;
        --BaseNamePattern) BaseNamePattern="$2"; shift 2 ;;
        --ExtraArgs) ExtraArgs="$2"; shift 2 ;;
        --singularity_bin) singularity_bin="$2"; shift 2 ;;
        --bbmap_sif) bbmap_sif="$2"; shift 2 ;;
        --bind) bind="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${Read1:-}" ]]; then echo "Error: --Read1 is required"; usage; fi
if [[ -z "${Read2:-}" ]]; then echo "Error: --Read2 is required"; usage; fi
if [[ -z "${PrimaryRef:-}" ]]; then echo "Error: --PrimaryRef is required"; usage; fi
if [[ -z "${OtherRefsArgs:-}" ]]; then echo "Error: --OtherRefsArgs is required"; usage; fi
if [[ -z "${OutputDir:-}" ]]; then echo "Error: --OutputDir is required"; usage; fi
if [[ -z "${MemoryMB:-}" ]]; then echo "Error: --MemoryMB is required"; usage; fi
if [[ -z "${Threads:-}" ]]; then echo "Error: --Threads is required"; usage; fi
if [[ -z "${BBMapPath:-}" ]]; then echo "Error: --BBMapPath is required"; usage; fi
if [[ -z "${BaseNamePattern:-}" ]]; then echo "Error: --BaseNamePattern is required"; usage; fi

# --- [Output Paths] ---

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${singularity_bin} exec -B ${bind} ${bbmap_sif} bbsplit.sh -Xmx${MemoryMB}M ref_primary=${PrimaryRef} ${OtherRefsArgs} in1=${Read1} in2=${Read2} basename=${BaseNamePattern} refstats=${OutputDir}/${SeqID}.refstats.txt threads=${Threads} ${ExtraArgs}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
mkdir -p "$(dirname "${OutputDir}")" 2>/dev/null || mkdir -p "${OutputDir}"

eval "$cmd"