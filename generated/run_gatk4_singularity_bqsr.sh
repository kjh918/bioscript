#!/usr/bin/env bash
set -euo pipefail

# Runner: gatk4 4.4.0.0 (singularity_bqsr)

# defaults (override by CLI)
SeqID=''
RecalInputBAM=''
ReferenceFasta=''
BamDir=''
qcResDir=''
KnownSnp=''
KnownIndel1=''
KnownIndel2=''
bind=/storage,/data
sif=/storage/images/gatk-4.4.0.0.sif
gc_threads=14
xmx_mb=16384
recal_table='[qcResDir]/[SeqID].recal.table.txt'
recal_bam='[BamDir]/[SeqID].recal.bam'
recal_bai='[BamDir]/[SeqID].recal.bam.bai'

CWD="."
PRINT_CMD=0
PRINT_OUTPUTS=0
EMIT_OUTPUTS=""

OUTPUT_KEYS=(recal_table recal_bam recal_bai)

usage() {
  echo "Usage: $0 --<Key> <val> ... [--cwd <dir>] [--print-cmd] [--print-outputs] [--emit-outputs <file>]" >&2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --SeqID) SeqID="$2"; shift 2;;
    --RecalInputBAM) RecalInputBAM="$2"; shift 2;;
    --ReferenceFasta) ReferenceFasta="$2"; shift 2;;
    --BamDir) BamDir="$2"; shift 2;;
    --qcResDir) qcResDir="$2"; shift 2;;
    --KnownSnp) KnownSnp="$2"; shift 2;;
    --KnownIndel1) KnownIndel1="$2"; shift 2;;
    --KnownIndel2) KnownIndel2="$2"; shift 2;;
    --bind) bind="$2"; shift 2;;
    --sif) sif="$2"; shift 2;;
    --gc_threads) gc_threads="$2"; shift 2;;
    --xmx_mb) xmx_mb="$2"; shift 2;;
    --recal_table) recal_table="$2"; shift 2;;
    --recal_bam) recal_bam="$2"; shift 2;;
    --recal_bai) recal_bai="$2"; shift 2;;
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
[[ -n "${RecalInputBAM}" ]] || { echo "Missing required --RecalInputBAM" >&2; exit 2; }
[[ -n "${ReferenceFasta}" ]] || { echo "Missing required --ReferenceFasta" >&2; exit 2; }
[[ -n "${BamDir}" ]] || { echo "Missing required --BamDir" >&2; exit 2; }
[[ -n "${qcResDir}" ]] || { echo "Missing required --qcResDir" >&2; exit 2; }
[[ -n "${KnownSnp}" ]] || { echo "Missing required --KnownSnp" >&2; exit 2; }
[[ -n "${KnownIndel1}" ]] || { echo "Missing required --KnownIndel1" >&2; exit 2; }
[[ -n "${KnownIndel2}" ]] || { echo "Missing required --KnownIndel2" >&2; exit 2; }

render() {
  local s="$1"
  s="${s//[SeqID]/${SeqID}}"
  s="${s//[RecalInputBAM]/${RecalInputBAM}}"
  s="${s//[ReferenceFasta]/${ReferenceFasta}}"
  s="${s//[BamDir]/${BamDir}}"
  s="${s//[qcResDir]/${qcResDir}}"
  s="${s//[KnownSnp]/${KnownSnp}}"
  s="${s//[KnownIndel1]/${KnownIndel1}}"
  s="${s//[KnownIndel2]/${KnownIndel2}}"
  s="${s//[bind]/${bind}}"
  s="${s//[sif]/${sif}}"
  s="${s//[gc_threads]/${gc_threads}}"
  s="${s//[xmx_mb]/${xmx_mb}}"
  s="${s//[recal_table]/${recal_table}}"
  s="${s//[recal_bam]/${recal_bam}}"
  s="${s//[recal_bai]/${recal_bai}}"
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
if [[ -z "${recal_table}" ]]; then
  __tmp='[qcResDir]/[SeqID].recal.table.txt'
  recal_table="$(render "$__tmp")"
fi
if [[ -z "${recal_bam}" ]]; then
  __tmp='[BamDir]/[SeqID].recal.bam'
  recal_bam="$(render "$__tmp")"
fi
if [[ -z "${recal_bai}" ]]; then
  __tmp='[BamDir]/[SeqID].recal.bam.bai'
  recal_bai="$(render "$__tmp")"
fi

CMD_LINE='singularity exec -B [bind] [sif] gatk BaseRecalibrator --java-options "-XX:ParallelGCThreads=[gc_threads] -Xmx[xmx_mb]m" --input [RecalInputBAM] --reference [ReferenceFasta] --output [qcResDir]/[SeqID].recal.table.txt --known-sites [KnownSnp] --known-sites [KnownIndel1] --known-sites [KnownIndel2] && singularity exec -B [bind] [sif] gatk ApplyBQSR --java-options "-XX:ParallelGCThreads=[gc_threads] -Xmx[xmx_mb]m" --input [RecalInputBAM] --bqsr-recal-file [qcResDir]/[SeqID].recal.table.txt --output [BamDir]/[SeqID].recal.bam'
CMD="$(render "$CMD_LINE")"

if [[ "$PRINT_CMD" -eq 1 ]]; then
  echo "$CMD"
  exit 0
fi

mkdir -p "$CWD"
bash -lc "cd \"$CWD\" && $CMD"

emit_outputs
