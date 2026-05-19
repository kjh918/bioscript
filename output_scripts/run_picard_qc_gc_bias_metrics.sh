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
    echo "  --BamDir          Directory containing the input BAM and results"
    echo "  --ReferenceFasta  Path to the reference genome FASTA file"
    echo ""
    echo "Optional Parameters:"
    echo "  --InputSuffix     Suffix of the input BAM (e.g., analysisReady, recal, dedup) (Default: analysisReady)"
    echo "  --java_bin        Path to java executable (Default: java)"
    echo "  --picard_jar      Path to Picard.jar file (Default: /storage/apps/bin/picard-3.1.0.jar)"
    echo "  --gc_threads      Number of threads for Java ParallelGC (Default: 14)"
    echo "  --xmx_mb          Maximum Java heap memory (MB) (Default: 16384)"
    echo "  --min_gc          Minimum GC content to include in metrics (Default: 0)"
    echo "  --max_gc          Maximum GC content to include in metrics (Default: 100)"
    echo "  --window_size     Window size for GC content calculation (Default: 100)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
BamDir=""
ReferenceFasta=""
InputSuffix="analysisReady"
java_bin="java"
picard_jar="/storage/apps/bin/picard-3.1.0.jar"
gc_threads="14"
xmx_mb="16384"
min_gc="0"
max_gc="100"
window_size="100"

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --BamDir) BamDir="$2"; shift 2 ;;
        --ReferenceFasta) ReferenceFasta="$2"; shift 2 ;;
        --InputSuffix) InputSuffix="$2"; shift 2 ;;
        --java_bin) java_bin="$2"; shift 2 ;;
        --picard_jar) picard_jar="$2"; shift 2 ;;
        --gc_threads) gc_threads="$2"; shift 2 ;;
        --xmx_mb) xmx_mb="$2"; shift 2 ;;
        --min_gc) min_gc="$2"; shift 2 ;;
        --max_gc) max_gc="$2"; shift 2 ;;
        --window_size) window_size="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${BamDir:-}" ]]; then echo "Error: --BamDir is required"; usage; fi
if [[ -z "${ReferenceFasta:-}" ]]; then echo "Error: --ReferenceFasta is required"; usage; fi

# --- [Output Paths] ---
gc_bias_metrics_txt="${BamDir}/${SeqID}.${InputSuffix}.gc_bias_metrics.txt"
gc_bias_summary_txt="${BamDir}/${SeqID}.${InputSuffix}.gc_bias_summary.txt"
gc_bias_chart_pdf="${BamDir}/${SeqID}.${InputSuffix}.gc_bias_metrics.pdf"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${java_bin} -XX:ParallelGCThreads=${gc_threads} -Xmx${xmx_mb}m -jar ${picard_jar} CollectGcBiasMetrics I=${BamDir}/${SeqID}.${InputSuffix}.bam O=${gc_bias_metrics_txt} S=${gc_bias_summary_txt} CHART=${gc_bias_chart_pdf} R=${ReferenceFasta} MINIMUM_GC=${min_gc} MAXIMUM_GC=${max_gc} WINDOW_SIZE=${window_size}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
mkdir -p "$(dirname "${BamDir}")" 2>/dev/null || mkdir -p "${BamDir}"

eval "$cmd"