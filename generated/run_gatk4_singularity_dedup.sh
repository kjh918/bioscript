#!/usr/bin/env bash
set -euo pipefail

# Runner: gatk4 4.4.0.0 (singularity_dedup)

# defaults (override by CLI)
SeqID=''
BamDir=''
qcResDir=''
bind=/storage,/data
sif=/storage/images/gatk-4.4.0.0.sif
gc_threads=14
xmx_mb=16384
md_args='--CREATE_INDEX true --REMOVE_SEQUENCING_DUPLICATES true'
dedup_bam='[BamDir]/[SeqID].sorted.dedup.bam'
dedup_bai='[BamDir]/[SeqID].sorted.dedup.bam.bai'
metrics_txt='[qcResDir]/[SeqID].mark.duplicates.metrics.txt'

CWD="."
PRINT_CMD=0
PRINT_OUTPUTS=0
EMIT_OUTPUTS=""

OUTPUT_KEYS=(dedup_bam dedup_bai metrics_txt)

usage() {
  echo "Usage: $0 --<Key> <val> ... [--cwd <dir>] [--print-cmd] [--print-outputs] [--emit-outputs <file>]" >&2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --SeqID) SeqID="$2"; shift 2;;
    --BamDir) BamDir="$2"; shift 2;;
    --qcResDir) qcResDir="$2"; shift 2;;
    --bind) bind="$2"; shift 2;;
    --sif) sif="$2"; shift 2;;
    --gc_threads) gc_threads="$2"; shift 2;;
    --xmx_mb) xmx_mb="$2"; shift 2;;
    --md_args) md_args="$2"; shift 2;;
    --dedup_bam) dedup_bam="$2"; shift 2;;
    --dedup_bai) dedup_bai="$2"; shift 2;;
    --metrics_txt) metrics_txt="$2"; shift 2;;
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
[[ -n "${qcResDir}" ]] || { echo "Missing required --qcResDir" >&2; exit 2; }

render() {
  local s="$1"
  s="${s//[SeqID]/${SeqID}}"
  s="${s//[BamDir]/${BamDir}}"
  s="${s//[qcResDir]/${qcResDir}}"
  s="${s//[bind]/${bind}}"
  s="${s//[sif]/${sif}}"
  s="${s//[gc_threads]/${gc_threads}}"
  s="${s//[xmx_mb]/${xmx_mb}}"
  s="${s//[md_args]/${md_args}}"
  s="${s//[dedup_bam]/${dedup_bam}}"
  s="${s//[dedup_bai]/${dedup_bai}}"
  s="${s//[metrics_txt]/${metrics_txt}}"
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
if [[ -z "${dedup_bam}" ]]; then
  __tmp='[BamDir]/[SeqID].sorted.dedup.bam'
  dedup_bam="$(render "$__tmp")"
fi
if [[ -z "${dedup_bai}" ]]; then
  __tmp='[BamDir]/[SeqID].sorted.dedup.bam.bai'
  dedup_bai="$(render "$__tmp")"
fi
if [[ -z "${metrics_txt}" ]]; then
  __tmp='[qcResDir]/[SeqID].mark.duplicates.metrics.txt'
  metrics_txt="$(render "$__tmp")"
fi

CMD_LINE='singularity exec -B [bind] [sif] gatk MarkDuplicates --java-options "-XX:ParallelGCThreads=[gc_threads] -Xmx[xmx_mb]m" --INPUT [BamDir]/[SeqID].sorted.bam --OUTPUT [BamDir]/[SeqID].sorted.dedup.bam --METRICS_FILE [qcResDir]/[SeqID].mark.duplicates.metrics.txt [md_args]'
CMD="$(render "$CMD_LINE")"

if [[ "$PRINT_CMD" -eq 1 ]]; then
  echo "$CMD"
  exit 0
fi

mkdir -p "$CWD"
bash -lc "cd \"$CWD\" && $CMD"

emit_outputs
