#!/usr/bin/env bash
set -euo pipefail

# Runner: fastp 0.23.4 (PicoplexGold_singularity_pe_trim)

# defaults (override by CLI)
SeqID=''
RawFastqDir=''
TrimFastqDir=''
qcResDir=''
Threads=8
sif=/storage/images/fastp-0.23.4.sif
bind=/storage,/data
length_required=100
average_qual=10
qualified_quality_phred=15
trim_front1=14
trim_front2=14
adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
out_read1='[TrimFastqDir]/[SeqID].trimmed_R1.fastq.gz'
out_read2='[TrimFastqDir]/[SeqID].trimmed_R2.fastq.gz'
json='[qcResDir]/[SeqID].fastp.json'
html='[qcResDir]/[SeqID].fastp.html'

CWD="."
PRINT_CMD=0
PRINT_OUTPUTS=0
EMIT_OUTPUTS=""

OUTPUT_KEYS=(out_read1 out_read2 json html)

usage() {
  echo "Usage: $0 --<Key> <val> ... [--cwd <dir>] [--print-cmd] [--print-outputs] [--emit-outputs <file>]" >&2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --SeqID) SeqID="$2"; shift 2;;
    --RawFastqDir) RawFastqDir="$2"; shift 2;;
    --TrimFastqDir) TrimFastqDir="$2"; shift 2;;
    --qcResDir) qcResDir="$2"; shift 2;;
    --Threads) Threads="$2"; shift 2;;
    --sif) sif="$2"; shift 2;;
    --bind) bind="$2"; shift 2;;
    --length_required) length_required="$2"; shift 2;;
    --average_qual) average_qual="$2"; shift 2;;
    --qualified_quality_phred) qualified_quality_phred="$2"; shift 2;;
    --trim_front1) trim_front1="$2"; shift 2;;
    --trim_front2) trim_front2="$2"; shift 2;;
    --adapter_sequence) adapter_sequence="$2"; shift 2;;
    --adapter_sequence_r2) adapter_sequence_r2="$2"; shift 2;;
    --out_read1) out_read1="$2"; shift 2;;
    --out_read2) out_read2="$2"; shift 2;;
    --json) json="$2"; shift 2;;
    --html) html="$2"; shift 2;;
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
[[ -n "${RawFastqDir}" ]] || { echo "Missing required --RawFastqDir" >&2; exit 2; }
[[ -n "${TrimFastqDir}" ]] || { echo "Missing required --TrimFastqDir" >&2; exit 2; }
[[ -n "${qcResDir}" ]] || { echo "Missing required --qcResDir" >&2; exit 2; }

render() {
  local s="$1"
  s="${s//[SeqID]/${SeqID}}"
  s="${s//[RawFastqDir]/${RawFastqDir}}"
  s="${s//[TrimFastqDir]/${TrimFastqDir}}"
  s="${s//[qcResDir]/${qcResDir}}"
  s="${s//[Threads]/${Threads}}"
  s="${s//[sif]/${sif}}"
  s="${s//[bind]/${bind}}"
  s="${s//[length_required]/${length_required}}"
  s="${s//[average_qual]/${average_qual}}"
  s="${s//[qualified_quality_phred]/${qualified_quality_phred}}"
  s="${s//[trim_front1]/${trim_front1}}"
  s="${s//[trim_front2]/${trim_front2}}"
  s="${s//[adapter_sequence]/${adapter_sequence}}"
  s="${s//[adapter_sequence_r2]/${adapter_sequence_r2}}"
  s="${s//[out_read1]/${out_read1}}"
  s="${s//[out_read2]/${out_read2}}"
  s="${s//[json]/${json}}"
  s="${s//[html]/${html}}"
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
if [[ -z "${out_read1}" ]]; then
  __tmp='[TrimFastqDir]/[SeqID].trimmed_R1.fastq.gz'
  out_read1="$(render "$__tmp")"
fi
if [[ -z "${out_read2}" ]]; then
  __tmp='[TrimFastqDir]/[SeqID].trimmed_R2.fastq.gz'
  out_read2="$(render "$__tmp")"
fi
if [[ -z "${json}" ]]; then
  __tmp='[qcResDir]/[SeqID].fastp.json'
  json="$(render "$__tmp")"
fi
if [[ -z "${html}" ]]; then
  __tmp='[qcResDir]/[SeqID].fastp.html'
  html="$(render "$__tmp")"
fi

CMD_LINE='singularity exec -B [bind] [sif] fastp --thread [Threads] --in1 [RawFastqDir]/[SeqID]_R1.fastq.gz --in2 [RawFastqDir]/[SeqID]_R2.fastq.gz --out1 [out_read1] --out2 [out_read2] --json [json] --html [html] --trim_poly_g --detect_adapter_for_pe --adapter_sequence [adapter_sequence] --adapter_sequence_r2 [adapter_sequence_r2] --length_required [length_required] --average_qual [average_qual] --qualified_quality_phred [qualified_quality_phred] --trim_front1 [trim_front1] --trim_front2 [trim_front2]'
CMD="$(render "$CMD_LINE")"

if [[ "$PRINT_CMD" -eq 1 ]]; then
  echo "$CMD"
  exit 0
fi

mkdir -p "$CWD"
bash -lc "cd \"$CWD\" && $CMD"

emit_outputs
