#!/usr/bin/env bash
set -euo pipefail

# --- Defaults ---
BamDir=''
ReferenceFasta=''
SeqID=''
gc_bias_chart_pdf='[BamDir]/[SeqID].gc.bias.metrics.pdf'
gc_bias_metrics_txt='[BamDir]/[SeqID].gc.bias.metrics.txt'
gc_bias_summary_txt='[BamDir]/[SeqID].gc.bias.summary.txt'
gc_threads=14
java_bin=java
max_gc=100
min_gc=0
picard_jar=/storage/apps/bin/picard-3.1.0.jar
window_size=100
xmx_mb=16384
CWD="."

# --- CLI Argument Parsing ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    --BamDir) BamDir="$2"; shift 2;;
    --ReferenceFasta) ReferenceFasta="$2"; shift 2;;
    --SeqID) SeqID="$2"; shift 2;;
    --gc_bias_chart_pdf) gc_bias_chart_pdf="$2"; shift 2;;
    --gc_bias_metrics_txt) gc_bias_metrics_txt="$2"; shift 2;;
    --gc_bias_summary_txt) gc_bias_summary_txt="$2"; shift 2;;
    --gc_threads) gc_threads="$2"; shift 2;;
    --java_bin) java_bin="$2"; shift 2;;
    --max_gc) max_gc="$2"; shift 2;;
    --min_gc) min_gc="$2"; shift 2;;
    --picard_jar) picard_jar="$2"; shift 2;;
    --window_size) window_size="$2"; shift 2;;
    --xmx_mb) xmx_mb="$2"; shift 2;;
    --cwd) CWD="$2"; shift 2;;
    -h|--help) echo "Usage: $0 [options]"; exit 0;;
    *) echo "Unknown argument: $1" >&2; exit 1;;
  esac
done

# --- Template Engine ---
render() {
  local s="$1"
  # 중첩된 변수 치환을 위해 3회 반복 (예: [A] -> [B] -> value)
  for i in {1..3}; do
    s="${s//\\[BamDir\\]/${BamDir}}"
    s="${s//\\[ReferenceFasta\\]/${ReferenceFasta}}"
    s="${s//\\[SeqID\\]/${SeqID}}"
    s="${s//\\[gc_bias_chart_pdf\\]/${gc_bias_chart_pdf}}"
    s="${s//\\[gc_bias_metrics_txt\\]/${gc_bias_metrics_txt}}"
    s="${s//\\[gc_bias_summary_txt\\]/${gc_bias_summary_txt}}"
    s="${s//\\[gc_threads\\]/${gc_threads}}"
    s="${s//\\[java_bin\\]/${java_bin}}"
    s="${s//\\[max_gc\\]/${max_gc}}"
    s="${s//\\[min_gc\\]/${min_gc}}"
    s="${s//\\[picard_jar\\]/${picard_jar}}"
    s="${s//\\[window_size\\]/${window_size}}"
    s="${s//\\[xmx_mb\\]/${xmx_mb}}"
  done
  # 앞뒤 공백 제거 및 연속 공백 축소
  echo "$s" | xargs
}

# --- Finalize Outputs & Command ---
gc_bias_metrics_txt=$(render "${gc_bias_metrics_txt}")
gc_bias_summary_txt=$(render "${gc_bias_summary_txt}")
gc_bias_chart_pdf=$(render "${gc_bias_chart_pdf}")
CMD_LINE='[java_bin] -XX:ParallelGCThreads=[gc_threads] -Xmx[xmx_mb]m -jar [picard_jar] CollectGcBiasMetrics I=[BamDir]/[SeqID].analysisReady.bam O=[BamDir]/[SeqID].gc.bias.metrics.txt S=[BamDir]/[SeqID].gc.bias.summary.txt CHART=[BamDir]/[SeqID].gc.bias.metrics.pdf R=[ReferenceFasta] MINIMUM_GC=[min_gc] MAXIMUM_GC=[max_gc] WINDOW_SIZE=[window_size]'
CMD=$(render "$CMD_LINE")

echo -e "\n[RUNNING CMD]\n$CMD\n"

# --- Execution ---
mkdir -p "$CWD"
cd "$CWD" && eval "$CMD"
