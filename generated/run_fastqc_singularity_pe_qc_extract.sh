#!/usr/bin/env bash
set -euo pipefail

# Runner: fastqc 0.12.1 (singularity_pe_qc_extract)

# defaults (override by CLI)
SeqID=''
RawFastqDir=''
qcResDir=''
bind=/storage,/data
sif=/storage/images/fastqc-0.12.1.sif
threads=8
fastqc_args=--extract
r1_html='[qcResDir]/[SeqID]_R1_fastqc.html'
r2_html='[qcResDir]/[SeqID]_R2_fastqc.html'
r1_zip='[qcResDir]/[SeqID]_R1_fastqc.zip'
r2_zip='[qcResDir]/[SeqID]_R2_fastqc.zip'
r1_dir='[qcResDir]/[SeqID]_R1_fastqc'
r2_dir='[qcResDir]/[SeqID]_R2_fastqc'

CWD="."
PRINT_CMD=0
PRINT_OUTPUTS=0
EMIT_OUTPUTS=""

OUTPUT_KEYS=(r1_html r2_html r1_zip r2_zip r1_dir r2_dir)

usage() {
  echo "Usage: $0 --<Key> <val> ... [--cwd <dir>] [--print-cmd] [--print-outputs] [--emit-outputs <file>]" >&2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --SeqID) SeqID="$2"; shift 2;;
    --RawFastqDir) RawFastqDir="$2"; shift 2;;
    --qcResDir) qcResDir="$2"; shift 2;;
    --bind) bind="$2"; shift 2;;
    --sif) sif="$2"; shift 2;;
    --threads) threads="$2"; shift 2;;
    --fastqc_args) fastqc_args="$2"; shift 2;;
    --r1_html) r1_html="$2"; shift 2;;
    --r2_html) r2_html="$2"; shift 2;;
    --r1_zip) r1_zip="$2"; shift 2;;
    --r2_zip) r2_zip="$2"; shift 2;;
    --r1_dir) r1_dir="$2"; shift 2;;
    --r2_dir) r2_dir="$2"; shift 2;;
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
[[ -n "${qcResDir}" ]] || { echo "Missing required --qcResDir" >&2; exit 2; }

render() {
  local s="$1"
  s="${s//[SeqID]/${SeqID}}"
  s="${s//[RawFastqDir]/${RawFastqDir}}"
  s="${s//[qcResDir]/${qcResDir}}"
  s="${s//[bind]/${bind}}"
  s="${s//[sif]/${sif}}"
  s="${s//[threads]/${threads}}"
  s="${s//[fastqc_args]/${fastqc_args}}"
  s="${s//[r1_html]/${r1_html}}"
  s="${s//[r2_html]/${r2_html}}"
  s="${s//[r1_zip]/${r1_zip}}"
  s="${s//[r2_zip]/${r2_zip}}"
  s="${s//[r1_dir]/${r1_dir}}"
  s="${s//[r2_dir]/${r2_dir}}"
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
if [[ -z "${r1_html}" ]]; then
  __tmp='[qcResDir]/[SeqID]_R1_fastqc.html'
  r1_html="$(render "$__tmp")"
fi
if [[ -z "${r2_html}" ]]; then
  __tmp='[qcResDir]/[SeqID]_R2_fastqc.html'
  r2_html="$(render "$__tmp")"
fi
if [[ -z "${r1_zip}" ]]; then
  __tmp='[qcResDir]/[SeqID]_R1_fastqc.zip'
  r1_zip="$(render "$__tmp")"
fi
if [[ -z "${r2_zip}" ]]; then
  __tmp='[qcResDir]/[SeqID]_R2_fastqc.zip'
  r2_zip="$(render "$__tmp")"
fi
if [[ -z "${r1_dir}" ]]; then
  __tmp='[qcResDir]/[SeqID]_R1_fastqc'
  r1_dir="$(render "$__tmp")"
fi
if [[ -z "${r2_dir}" ]]; then
  __tmp='[qcResDir]/[SeqID]_R2_fastqc'
  r2_dir="$(render "$__tmp")"
fi

CMD_LINE='singularity exec -B [bind] [sif] fastqc [fastqc_args] --threads [threads] --outdir [qcResDir] [RawFastqDir]/[SeqID]_R1.fastq.gz [RawFastqDir]/[SeqID]_R2.fastq.gz'
CMD="$(render "$CMD_LINE")"

if [[ "$PRINT_CMD" -eq 1 ]]; then
  echo "$CMD"
  exit 0
fi

mkdir -p "$CWD"
bash -lc "cd \"$CWD\" && $CMD"

emit_outputs
