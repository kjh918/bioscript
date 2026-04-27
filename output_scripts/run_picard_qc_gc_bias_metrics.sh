#!/bin/bash
# [METADATA]
# TOOL_NAME = picard
# VERSION = 3.1.0
# THREADS = 1

# Tool Info: picard (3.1.0)
# Profile: qc_gc_bias_metrics

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           Sequence identifier used for file naming"
    echo "  --BamDir          Directory containing the sorted input BAM files"
    echo "  --qcResDir        Directory where duplication metrics will be saved"
    echo ""
    echo "Optional Parameters:"
    echo "  --java_bin        Path to java executable (Default: java)"
    echo "  --picard_jar      Path to picard.jar file (Default: /storage/apps/bin/picard.jar)"
    echo "  --Threads         Number of parallel GC threads (-XX:ParallelGCThreads) (Default: 14)"
    echo "  --Memory          Maximum Java heap size (-Xmx) (Default: 16384m)"
    echo "  --TmpDir          Directory for Picard temporary files (Default: /tmp)"
    echo "  --create_index    Automatically create BAM index (.bai) (Default: true)"
    echo "  --remove_duplicates If true, removes duplicates completely rather than just marking them (Default: true)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
BamDir=""
qcResDir=""
java_bin="java"
picard_jar="/storage/apps/bin/picard.jar"
Threads="14"
Memory="16384m"
TmpDir="/tmp"
create_index="true"
remove_duplicates="true"

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --BamDir) BamDir="$2"; shift 2 ;;
        --qcResDir) qcResDir="$2"; shift 2 ;;
        --java_bin) java_bin="$2"; shift 2 ;;
        --picard_jar) picard_jar="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        --Memory) Memory="$2"; shift 2 ;;
        --TmpDir) TmpDir="$2"; shift 2 ;;
        --create_index) create_index="$2"; shift 2 ;;
        --remove_duplicates) remove_duplicates="$2"; shift 2 ;;
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
out_bam="${BamDir}/${SeqID}.sorted.dedup.bam"
metrics="${qcResDir}/${SeqID}.mark.duplicates.metrics.txt"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${java_bin} -XX:ParallelGCThreads=${Threads} -Xmx${Memory} -jar ${picard_jar} MarkDuplicates  INPUT=${BamDir}/${SeqID}.sorted.bam  OUTPUT=${out_bam}  METRICS_FILE=${metrics}  CREATE_INDEX=${create_index}  REMOVE_DUPLICATES=${remove_duplicates}  TMP_DIR=${TmpDir}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
mkdir -p "$(dirname "${BamDir}")" 2>/dev/null || mkdir -p "${BamDir}"
mkdir -p "$(dirname "${qcResDir}")" 2>/dev/null || mkdir -p "${qcResDir}"
mkdir -p "$(dirname "${TmpDir}")" 2>/dev/null || mkdir -p "${TmpDir}"

eval "$cmd"