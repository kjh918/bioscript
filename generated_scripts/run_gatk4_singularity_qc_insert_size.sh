#!/usr/bin/env bash
set -euo pipefail

# --- Defaults ---
BamDir=''
ReferenceFasta=''
SeqID=''
bind=/storage,/data
gatk_jar=/gatk/gatk-package-4.4.0.0-local.jar
gc_threads=14
insert_size_hist_pdf='[qcResDir]/[SeqID].insert.size.histogram.pdf'
insert_size_metrics_txt='[qcResDir]/[SeqID].insert.size.metrics.txt'
java_bin=java
lc_all=en_US.UTF-8
qcResDir=''
sif=/storage/images/gatk-4.4.0.0.sif
xmx_mb=16384
CWD="."

# --- CLI Argument Parsing ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    --BamDir) BamDir="$2"; shift 2;;
    --ReferenceFasta) ReferenceFasta="$2"; shift 2;;
    --SeqID) SeqID="$2"; shift 2;;
    --bind) bind="$2"; shift 2;;
    --gatk_jar) gatk_jar="$2"; shift 2;;
    --gc_threads) gc_threads="$2"; shift 2;;
    --insert_size_hist_pdf) insert_size_hist_pdf="$2"; shift 2;;
    --insert_size_metrics_txt) insert_size_metrics_txt="$2"; shift 2;;
    --java_bin) java_bin="$2"; shift 2;;
    --lc_all) lc_all="$2"; shift 2;;
    --qcResDir) qcResDir="$2"; shift 2;;
    --sif) sif="$2"; shift 2;;
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
    s="${s//\\[bind\\]/${bind}}"
    s="${s//\\[gatk_jar\\]/${gatk_jar}}"
    s="${s//\\[gc_threads\\]/${gc_threads}}"
    s="${s//\\[insert_size_hist_pdf\\]/${insert_size_hist_pdf}}"
    s="${s//\\[insert_size_metrics_txt\\]/${insert_size_metrics_txt}}"
    s="${s//\\[java_bin\\]/${java_bin}}"
    s="${s//\\[lc_all\\]/${lc_all}}"
    s="${s//\\[qcResDir\\]/${qcResDir}}"
    s="${s//\\[sif\\]/${sif}}"
    s="${s//\\[xmx_mb\\]/${xmx_mb}}"
  done
  # 앞뒤 공백 제거 및 연속 공백 축소
  echo "$s" | xargs
}

# --- Finalize Outputs & Command ---
insert_size_metrics_txt=$(render "${insert_size_metrics_txt}")
insert_size_hist_pdf=$(render "${insert_size_hist_pdf}")
CMD_LINE='export LC_ALL=[lc_all] && singularity exec -B [bind] [sif] [java_bin] -XX:ParallelGCThreads=[gc_threads] -Xmx[xmx_mb]m -jar [gatk_jar] CollectInsertSizeMetrics --INPUT [BamDir]/[SeqID].analysisReady.bam --OUTPUT [qcResDir]/[SeqID].insert.size.metrics.txt --Histogram_FILE [qcResDir]/[SeqID].insert.size.histogram.pdf --REFERENCE_SEQUENCE [ReferenceFasta]'
CMD=$(render "$CMD_LINE")

echo -e "\n[RUNNING CMD]\n$CMD\n"

# --- Execution ---
mkdir -p "$CWD"
cd "$CWD" && eval "$CMD"
