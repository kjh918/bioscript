#!/bin/bash
# [METADATA]
# TOOL_NAME = sortmerna
# VERSION = 4.3.7
# THREADS = 1

# Tool Info: sortmerna (4.3.7)
# Profile: rrna_filtering

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           No description"
    echo "  --FastqDir        No description"
    echo "  --qcResDir        Log 및 결과 파일 저장 경로"
    echo "  --RefArgs         List of --ref references"
    echo ""
    echo "Optional Parameters:"
    echo "  --IndexDir        Pre-built index directory (--idx-dir) (Default: )"
    echo "  --Threads         No description (Default: 8)"
    echo "  --InputSuffix     No description (Default: fastq.gz)"
    echo "  --paired_cmd      No description (Default: --paired_in --out2)"
    echo "  --extra_args      No description (Default: )"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
FastqDir=""
qcResDir=""
RefArgs=""
IndexDir=""
Threads="8"
InputSuffix="fastq.gz"
paired_cmd="--paired_in --out2"
extra_args=""

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --FastqDir) FastqDir="$2"; shift 2 ;;
        --qcResDir) qcResDir="$2"; shift 2 ;;
        --RefArgs) RefArgs="$2"; shift 2 ;;
        --IndexDir) IndexDir="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        --InputSuffix) InputSuffix="$2"; shift 2 ;;
        --paired_cmd) paired_cmd="$2"; shift 2 ;;
        --extra_args) extra_args="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${FastqDir:-}" ]]; then echo "Error: --FastqDir is required"; usage; fi
if [[ -z "${qcResDir:-}" ]]; then echo "Error: --qcResDir is required"; usage; fi
if [[ -z "${RefArgs:-}" ]]; then echo "Error: --RefArgs is required"; usage; fi

# --- [Output Paths] ---
non_rrna_r1="${qcResDir}/${SeqID}_1.non_rRNA.fastq.gz"
non_rrna_r2="${qcResDir}/${SeqID}_2.non_rRNA.fastq.gz"
sortmerna_log="${qcResDir}/${SeqID}.sortmerna.log"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="sortmerna ${RefArgs} --reads ${FastqDir}/${SeqID}_1.${InputSuffix} --reads ${FastqDir}/${SeqID}_2.${InputSuffix} --threads ${Threads} --workdir . --aligned ${qcResDir}/rRNA_reads --fastx --other ${qcResDir}/non_rRNA_reads ${paired_cmd} ${extra_args} && mv ${qcResDir}/non_rRNA_reads_fwd.f*q.gz ${non_rrna_r1} && mv ${qcResDir}/non_rRNA_reads_rev.f*q.gz ${non_rrna_r2} && mv ${qcResDir}/rRNA_reads.log ${sortmerna_log}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
mkdir -p "$(dirname "${FastqDir}")" 2>/dev/null || mkdir -p "${FastqDir}"
mkdir -p "$(dirname "${qcResDir}")" 2>/dev/null || mkdir -p "${qcResDir}"
mkdir -p "$(dirname "${IndexDir}")" 2>/dev/null || mkdir -p "${IndexDir}"

eval "$cmd"