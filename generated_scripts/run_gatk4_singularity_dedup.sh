#!/usr/bin/env bash
set -euo pipefail

# --- Defaults ---
BamDir=''
SeqID=''
bind=/storage,/data
dedup_bai='[BamDir]/[SeqID].sorted.dedup.bam.bai'
dedup_bam='[BamDir]/[SeqID].sorted.dedup.bam'
gc_threads=14
md_args='--CREATE_INDEX true --REMOVE_SEQUENCING_DUPLICATES true'
metrics_txt='[qcResDir]/[SeqID].mark.duplicates.metrics.txt'
qcResDir=''
sif=/storage/images/gatk-4.4.0.0.sif
xmx_mb=16384
CWD="."

# --- CLI Argument Parsing ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    --BamDir) BamDir="$2"; shift 2;;
    --SeqID) SeqID="$2"; shift 2;;
    --bind) bind="$2"; shift 2;;
    --dedup_bai) dedup_bai="$2"; shift 2;;
    --dedup_bam) dedup_bam="$2"; shift 2;;
    --gc_threads) gc_threads="$2"; shift 2;;
    --md_args) md_args="$2"; shift 2;;
    --metrics_txt) metrics_txt="$2"; shift 2;;
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
    s="${s//\\[SeqID\\]/${SeqID}}"
    s="${s//\\[bind\\]/${bind}}"
    s="${s//\\[dedup_bai\\]/${dedup_bai}}"
    s="${s//\\[dedup_bam\\]/${dedup_bam}}"
    s="${s//\\[gc_threads\\]/${gc_threads}}"
    s="${s//\\[md_args\\]/${md_args}}"
    s="${s//\\[metrics_txt\\]/${metrics_txt}}"
    s="${s//\\[qcResDir\\]/${qcResDir}}"
    s="${s//\\[sif\\]/${sif}}"
    s="${s//\\[xmx_mb\\]/${xmx_mb}}"
  done
  # 앞뒤 공백 제거 및 연속 공백 축소
  echo "$s" | xargs
}

# --- Finalize Outputs & Command ---
dedup_bam=$(render "${dedup_bam}")
dedup_bai=$(render "${dedup_bai}")
metrics_txt=$(render "${metrics_txt}")
CMD_LINE='singularity exec -B [bind] [sif] gatk MarkDuplicates --java-options "-XX:ParallelGCThreads=[gc_threads] -Xmx[xmx_mb]m" --INPUT [BamDir]/[SeqID].sorted.bam --OUTPUT [BamDir]/[SeqID].sorted.dedup.bam --METRICS_FILE [qcResDir]/[SeqID].mark.duplicates.metrics.txt [md_args]'
CMD=$(render "$CMD_LINE")

echo -e "\n[RUNNING CMD]\n$CMD\n"

# --- Execution ---
mkdir -p "$CWD"
cd "$CWD" && eval "$CMD"
