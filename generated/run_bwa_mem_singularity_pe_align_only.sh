#!/usr/bin/env bash
set -euo pipefail

# Runner: bwa_mem 0.7.17 (singularity_pe_align_only)

# defaults (override by CLI)
SeqID=''
TrimFastqDir=''
BamDir=''
BwaIndex=''
ReadGroupID=''
ReadGroupPlatform=''
ReadGroupLibrary=''
ReadGroupCenter=''
bind=/storage,/data
sif=/storage/images/bwa-0.7.17.sif
Threads=8
bwa_args='-M -Y -L 50,50'
primary_bam='[BamDir]/[SeqID].primary.bam'

CWD="."
PRINT_CMD=0
PRINT_OUTPUTS=0
EMIT_OUTPUTS=""

OUTPUT_KEYS=(primary_bam)

usage() {
  echo "Usage: $0 --<Key> <val> ... [--cwd <dir>] [--print-cmd] [--print-outputs] [--emit-outputs <file>]" >&2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --SeqID) SeqID="$2"; shift 2;;
    --TrimFastqDir) TrimFastqDir="$2"; shift 2;;
    --BamDir) BamDir="$2"; shift 2;;
    --BwaIndex) BwaIndex="$2"; shift 2;;
    --ReadGroupID) ReadGroupID="$2"; shift 2;;
    --ReadGroupPlatform) ReadGroupPlatform="$2"; shift 2;;
    --ReadGroupLibrary) ReadGroupLibrary="$2"; shift 2;;
    --ReadGroupCenter) ReadGroupCenter="$2"; shift 2;;
    --bind) bind="$2"; shift 2;;
    --sif) sif="$2"; shift 2;;
    --Threads) Threads="$2"; shift 2;;
    --bwa_args) bwa_args="$2"; shift 2;;
    --primary_bam) primary_bam="$2"; shift 2;;
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
[[ -n "${TrimFastqDir}" ]] || { echo "Missing required --TrimFastqDir" >&2; exit 2; }
[[ -n "${BamDir}" ]] || { echo "Missing required --BamDir" >&2; exit 2; }
[[ -n "${BwaIndex}" ]] || { echo "Missing required --BwaIndex" >&2; exit 2; }
[[ -n "${ReadGroupID}" ]] || { echo "Missing required --ReadGroupID" >&2; exit 2; }
[[ -n "${ReadGroupPlatform}" ]] || { echo "Missing required --ReadGroupPlatform" >&2; exit 2; }
[[ -n "${ReadGroupLibrary}" ]] || { echo "Missing required --ReadGroupLibrary" >&2; exit 2; }
[[ -n "${ReadGroupCenter}" ]] || { echo "Missing required --ReadGroupCenter" >&2; exit 2; }

render() {
  local s="$1"
  s="${s//[SeqID]/${SeqID}}"
  s="${s//[TrimFastqDir]/${TrimFastqDir}}"
  s="${s//[BamDir]/${BamDir}}"
  s="${s//[BwaIndex]/${BwaIndex}}"
  s="${s//[ReadGroupID]/${ReadGroupID}}"
  s="${s//[ReadGroupPlatform]/${ReadGroupPlatform}}"
  s="${s//[ReadGroupLibrary]/${ReadGroupLibrary}}"
  s="${s//[ReadGroupCenter]/${ReadGroupCenter}}"
  s="${s//[bind]/${bind}}"
  s="${s//[sif]/${sif}}"
  s="${s//[Threads]/${Threads}}"
  s="${s//[bwa_args]/${bwa_args}}"
  s="${s//[primary_bam]/${primary_bam}}"
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
if [[ -z "${primary_bam}" ]]; then
  __tmp='[BamDir]/[SeqID].primary.bam'
  primary_bam="$(render "$__tmp")"
fi

CMD_LINE='singularity exec -B [bind] [sif] bwa mem [bwa_args] -t [Threads] -R "@RG\tID:[ReadGroupID]\tPL:[ReadGroupPlatform]\tLB:[ReadGroupLibrary]\tSM:[SeqID]\tCN:[ReadGroupCenter]" [BwaIndex] [TrimFastqDir]/[SeqID].trimmed_R1.fastq.gz [TrimFastqDir]/[SeqID].trimmed_R2.fastq.gz | samtools view -bS -o [BamDir]/[SeqID].primary.bam -'
CMD="$(render "$CMD_LINE")"

if [[ "$PRINT_CMD" -eq 1 ]]; then
  echo "$CMD"
  exit 0
fi

mkdir -p "$CWD"
bash -lc "cd \"$CWD\" && $CMD"

emit_outputs
