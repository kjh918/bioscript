#!/usr/bin/env bash
set -euo pipefail

# --- Defaults ---
RawFastqDir=''
SeqID=''
bind=/storage,/data
fastqc_args=--extract
qcResDir=''
r1_dir='[qcResDir]/[SeqID]_R1_fastqc'
r1_html='[qcResDir]/[SeqID]_R1_fastqc.html'
r1_zip='[qcResDir]/[SeqID]_R1_fastqc.zip'
r2_dir='[qcResDir]/[SeqID]_R2_fastqc'
r2_html='[qcResDir]/[SeqID]_R2_fastqc.html'
r2_zip='[qcResDir]/[SeqID]_R2_fastqc.zip'
sif=/storage/images/fastqc-0.12.1.sif
threads=8
CWD="."

# --- CLI Argument Parsing ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    --RawFastqDir) RawFastqDir="$2"; shift 2;;
    --SeqID) SeqID="$2"; shift 2;;
    --bind) bind="$2"; shift 2;;
    --fastqc_args) fastqc_args="$2"; shift 2;;
    --qcResDir) qcResDir="$2"; shift 2;;
    --r1_dir) r1_dir="$2"; shift 2;;
    --r1_html) r1_html="$2"; shift 2;;
    --r1_zip) r1_zip="$2"; shift 2;;
    --r2_dir) r2_dir="$2"; shift 2;;
    --r2_html) r2_html="$2"; shift 2;;
    --r2_zip) r2_zip="$2"; shift 2;;
    --sif) sif="$2"; shift 2;;
    --threads) threads="$2"; shift 2;;
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
    s="${s//\\[RawFastqDir\\]/${RawFastqDir}}"
    s="${s//\\[SeqID\\]/${SeqID}}"
    s="${s//\\[bind\\]/${bind}}"
    s="${s//\\[fastqc_args\\]/${fastqc_args}}"
    s="${s//\\[qcResDir\\]/${qcResDir}}"
    s="${s//\\[r1_dir\\]/${r1_dir}}"
    s="${s//\\[r1_html\\]/${r1_html}}"
    s="${s//\\[r1_zip\\]/${r1_zip}}"
    s="${s//\\[r2_dir\\]/${r2_dir}}"
    s="${s//\\[r2_html\\]/${r2_html}}"
    s="${s//\\[r2_zip\\]/${r2_zip}}"
    s="${s//\\[sif\\]/${sif}}"
    s="${s//\\[threads\\]/${threads}}"
  done
  # 앞뒤 공백 제거 및 연속 공백 축소
  echo "$s" | xargs
}

# --- Finalize Outputs & Command ---
r1_html=$(render "${r1_html}")
r2_html=$(render "${r2_html}")
r1_zip=$(render "${r1_zip}")
r2_zip=$(render "${r2_zip}")
r1_dir=$(render "${r1_dir}")
r2_dir=$(render "${r2_dir}")
CMD_LINE='singularity exec -B [bind] [sif] fastqc [fastqc_args] --threads [threads] --outdir [qcResDir] [RawFastqDir]/[SeqID]_R1.fastq.gz [RawFastqDir]/[SeqID]_R2.fastq.gz'
CMD=$(render "$CMD_LINE")

echo -e "\n[RUNNING CMD]\n$CMD\n"

# --- Execution ---
mkdir -p "$CWD"
cd "$CWD" && eval "$CMD"
