#!/usr/bin/env bash
set -euo pipefail

# --- Defaults ---
BamDir=''
SeqID=''
Threads=8
sorted_bai='[BamDir]/[SeqID].sorted.bam.bai'
sorted_bam='[BamDir]/[SeqID].sorted.bam'
CWD="."

# --- CLI Argument Parsing ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    --BamDir) BamDir="$2"; shift 2;;
    --SeqID) SeqID="$2"; shift 2;;
    --Threads) Threads="$2"; shift 2;;
    --sorted_bai) sorted_bai="$2"; shift 2;;
    --sorted_bam) sorted_bam="$2"; shift 2;;
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
    s="${s//\\[Threads\\]/${Threads}}"
    s="${s//\\[sorted_bai\\]/${sorted_bai}}"
    s="${s//\\[sorted_bam\\]/${sorted_bam}}"
  done
  # 앞뒤 공백 제거 및 연속 공백 축소
  echo "$s" | xargs
}

# --- Finalize Outputs & Command ---
sorted_bam=$(render "${sorted_bam}")
sorted_bai=$(render "${sorted_bai}")
CMD_LINE='samtools sort -@ [Threads] -o [BamDir]/[SeqID].sorted.bam [BamDir]/[SeqID].primary.bam && samtools index -b -@ [Threads] [BamDir]/[SeqID].sorted.bam'
CMD=$(render "$CMD_LINE")

echo -e "\n[RUNNING CMD]\n$CMD\n"

# --- Execution ---
mkdir -p "$CWD"
cd "$CWD" && eval "$CMD"
