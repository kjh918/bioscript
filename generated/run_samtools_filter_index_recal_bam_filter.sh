#!/usr/bin/env bash
set -euo pipefail

# Runner: samtools_filter_index 1.10 (recal_bam_filter)

# defaults (override by CLI)
SeqID=''
BamDir=''
Threads=8
min_mapq=20
include_flag=0x2
exclude_flag1=0x100
exclude_flag2=0x4
expr_sclen='sclen < 20'
expr_nm='[NM] < 12'
filtered_bam='[BamDir]/[SeqID].recal.filtered.bam'
filtered_bai='[BamDir]/[SeqID].recal.filtered.bam.bai'

CWD="."
PRINT_CMD=0
PRINT_OUTPUTS=0
EMIT_OUTPUTS=""

OUTPUT_KEYS=(filtered_bam filtered_bai)

usage() {
  echo "Usage: $0 --<Key> <val> ... [--cwd <dir>] [--print-cmd] [--print-outputs] [--emit-outputs <file>]" >&2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --SeqID) SeqID="$2"; shift 2;;
    --BamDir) BamDir="$2"; shift 2;;
    --Threads) Threads="$2"; shift 2;;
    --min_mapq) min_mapq="$2"; shift 2;;
    --include_flag) include_flag="$2"; shift 2;;
    --exclude_flag1) exclude_flag1="$2"; shift 2;;
    --exclude_flag2) exclude_flag2="$2"; shift 2;;
    --expr_sclen) expr_sclen="$2"; shift 2;;
    --expr_nm) expr_nm="$2"; shift 2;;
    --filtered_bam) filtered_bam="$2"; shift 2;;
    --filtered_bai) filtered_bai="$2"; shift 2;;
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
  s="${s//[min_mapq]/${min_mapq}}"
  s="${s//[include_flag]/${include_flag}}"
  s="${s//[exclude_flag1]/${exclude_flag1}}"
  s="${s//[exclude_flag2]/${exclude_flag2}}"
  s="${s//[expr_sclen]/${expr_sclen}}"
  s="${s//[expr_nm]/${expr_nm}}"
  s="${s//[filtered_bam]/${filtered_bam}}"
  s="${s//[filtered_bai]/${filtered_bai}}"
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
if [[ -z "${filtered_bam}" ]]; then
  __tmp='[BamDir]/[SeqID].recal.filtered.bam'
  filtered_bam="$(render "$__tmp")"
fi
if [[ -z "${filtered_bai}" ]]; then
  __tmp='[BamDir]/[SeqID].recal.filtered.bam.bai'
  filtered_bai="$(render "$__tmp")"
fi

CMD_LINE='samtools view -b -h -q [min_mapq] -f [include_flag] -F [exclude_flag1] -F [exclude_flag2] -e "[expr_sclen]" -e "[expr_nm]" --threads [Threads] [BamDir]/[SeqID].recal.bam > [BamDir]/[SeqID].recal.filtered.bam && samtools index --threads [Threads] [BamDir]/[SeqID].recal.filtered.bam'
CMD="$(render "$CMD_LINE")"

if [[ "$PRINT_CMD" -eq 1 ]]; then
  echo "$CMD"
  exit 0
fi

mkdir -p "$CWD"
bash -lc "cd \"$CWD\" && $CMD"

emit_outputs
