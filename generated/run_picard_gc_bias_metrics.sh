#!/usr/bin/env bash
set -euo pipefail

# Runner: picard 3.1.0 (gc_bias_metrics)

# defaults (override by CLI)
SeqID=''
BamDir=''
ReferenceFasta=''
java_bin=java
picard_jar=/storage/apps/bin/picard-3.1.0.jar
gc_threads=14
xmx_mb=16384
min_gc=0
max_gc=100
window_size=100
gc_bias_metrics_txt='[BamDir]/[SeqID].gc.bias.metrics.txt'
gc_bias_summary_txt='[BamDir]/[SeqID].gc.bias.summary.txt'
gc_bias_chart_pdf='[BamDir]/[SeqID].gc.bias.metrics.pdf'

CWD="."
PRINT_CMD=0
PRINT_OUTPUTS=0
EMIT_OUTPUTS=""

OUTPUT_KEYS=(gc_bias_metrics_txt gc_bias_summary_txt gc_bias_chart_pdf)

usage() {
  echo "Usage: $0 --<Key> <val> ... [--cwd <dir>] [--print-cmd] [--print-outputs] [--emit-outputs <file>]" >&2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --SeqID) SeqID="$2"; shift 2;;
    --BamDir) BamDir="$2"; shift 2;;
    --ReferenceFasta) ReferenceFasta="$2"; shift 2;;
    --java_bin) java_bin="$2"; shift 2;;
    --picard_jar) picard_jar="$2"; shift 2;;
    --gc_threads) gc_threads="$2"; shift 2;;
    --xmx_mb) xmx_mb="$2"; shift 2;;
    --min_gc) min_gc="$2"; shift 2;;
    --max_gc) max_gc="$2"; shift 2;;
    --window_size) window_size="$2"; shift 2;;
    --gc_bias_metrics_txt) gc_bias_metrics_txt="$2"; shift 2;;
    --gc_bias_summary_txt) gc_bias_summary_txt="$2"; shift 2;;
    --gc_bias_chart_pdf) gc_bias_chart_pdf="$2"; shift 2;;
    --cwd) CWD="$2"; shift 2;;
    --print-cmd) PRINT_CMD=1; shift;;
    --print-outputs) PRINT_OUTPUTS=1; shift;;
    --emit-outputs) EMIT_OUTPUTS="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1" >&2; usage; exit 2;;
  esac
done

# required inputs
[[ -n "${SeqID}" ]] || { echo "Missing required --SeqID" >&2; exit 2; }
[[ -n "${BamDir}" ]] || { echo "Missing required --BamDir" >&2; exit 2; }
[[ -n "${ReferenceFasta}" ]] || { echo "Missing required --ReferenceFasta" >&2; exit 2; }

render() {
  local s="$1"
  s="${s//[SeqID]/${SeqID}}"
  s="${s//[BamDir]/${BamDir}}"
  s="${s//[ReferenceFasta]/${ReferenceFasta}}"
  s="${s//[java_bin]/${java_bin}}"
  s="${s//[picard_jar]/${picard_jar}}"
  s="${s//[gc_threads]/${gc_threads}}"
  s="${s//[xmx_mb]/${xmx_mb}}"
  s="${s//[min_gc]/${min_gc}}"
  s="${s//[max_gc]/${max_gc}}"
  s="${s//[window_size]/${window_size}}"
  s="${s//[gc_bias_metrics_txt]/${gc_bias_metrics_txt}}"
  s="${s//[gc_bias_summary_txt]/${gc_bias_summary_txt}}"
  s="${s//[gc_bias_chart_pdf]/${gc_bias_chart_pdf}}"
  s="$(echo "$s" | tr -s ' ' | sed 's/^ *//; s/ *$//')"
  echo "$s"
}

emit_outputs() {
  local out=""
  for k in "${OUTPUT_KEYS[@]}"; do
    out+="${k}\t${!k}\n"
  done
  if [[ "$PRINT_OUTPUTS" -eq 1 ]]; then
    printf "%b" "$out"
  fi
  if [[ -n "$EMIT_OUTPUTS" ]]; then
    printf "%b" "$out" > "$EMIT_OUTPUTS"
  fi
}

# finalize outputs
if [[ -z "${gc_bias_metrics_txt}" ]]; then
  __tmp='[BamDir]/[SeqID].gc.bias.metrics.txt'
  gc_bias_metrics_txt="$(render "$__tmp")"
fi
if [[ -z "${gc_bias_summary_txt}" ]]; then
  __tmp='[BamDir]/[SeqID].gc.bias.summary.txt'
  gc_bias_summary_txt="$(render "$__tmp")"
fi
if [[ -z "${gc_bias_chart_pdf}" ]]; then
  __tmp='[BamDir]/[SeqID].gc.bias.metrics.pdf'
  gc_bias_chart_pdf="$(render "$__tmp")"
fi

CMD_LINE='[java_bin] -XX:ParallelGCThreads=[gc_threads] -Xmx[xmx_mb]m -jar [picard_jar] CollectGcBiasMetrics I=[BamDir]/[SeqID].analysisReady.bam O=[BamDir]/[SeqID].gc.bias.metrics.txt S=[BamDir]/[SeqID].gc.bias.summary.txt CHART=[BamDir]/[SeqID].gc.bias.metrics.pdf R=[ReferenceFasta] MINIMUM_GC=[min_gc] MAXIMUM_GC=[max_gc] WINDOW_SIZE=[window_size]'
CMD="$(render "$CMD_LINE")"

if [[ "$PRINT_CMD" -eq 1 ]]; then
  echo "$CMD"
  exit 0
fi

mkdir -p "$CWD"
bash -lc "cd \"$CWD\" && $CMD"

emit_outputs
