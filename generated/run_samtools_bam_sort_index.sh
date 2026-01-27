#!/usr/bin/env bash
set -euo pipefail

# Runner: samtools 1.10 (bam_sort_index)

# defaults (override by CLI)
SeqID=''
BamDir=''
Threads=8
sorted_bam='[BamDir]/[SeqID].sorted.bam'
sorted_bai='[BamDir]/[SeqID].sorted.bam.bai'

CWD="."
PRINT_CMD=0
PRINT_OUTPUTS=0
EMIT_OUTPUTS=""

OUTPUT_KEYS=(sorted_bam sorted_bai)

usage() {
  echo "Usage: $0 --<Key> <val> ... [--cwd <dir>] [--print-cmd] [--print-outputs] [--emit-outputs <file>]" >&2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --SeqID) SeqID="$2"; shift 2;;
    --BamDir) BamDir="$2"; shift 2;;
    --Threads) Threads="$2"; shift 2;;
    --sorted_bam) sorted_bam="$2"; shift 2;;
    --sorted_bai) sorted_bai="$2"; shift 2;;
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

render() {
  local s="$1"
  s="${s//[SeqID]/${SeqID}}"
  s="${s//[BamDir]/${BamDir}}"
  s="${s//[Threads]/${Threads}}"
  s="${s//[sorted_bam]/${sorted_bam}}"
  s="${s//[sorted_bai]/${sorted_bai}}"
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
if [[ -z "${sorted_bam}" ]]; then
  __tmp='[BamDir]/[SeqID].sorted.bam'
  sorted_bam="$(render "$__tmp")"
fi
if [[ -z "${sorted_bai}" ]]; then
  __tmp='[BamDir]/[SeqID].sorted.bam.bai'
  sorted_bai="$(render "$__tmp")"
fi

CMD_LINE='samtools sort -@ [Threads] -o [BamDir]/[SeqID].sorted.bam [BamDir]/[SeqID].primary.bam && samtools index -b -@ [Threads] [BamDir]/[SeqID].sorted.bam'
CMD="$(render "$CMD_LINE")"

if [[ "$PRINT_CMD" -eq 1 ]]; then
  echo "$CMD"
  exit 0
fi

mkdir -p "$CWD"
bash -lc "cd \"$CWD\" && $CMD"

emit_outputs
