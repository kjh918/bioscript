#!/usr/bin/env bash
set -euo pipefail

# --- Defaults ---
BamDir=''
SeqID=''
Threads=8
analysis_ready_bai='[BamDir]/[SeqID].analysisReady.bam.bai'
analysis_ready_bam='[BamDir]/[SeqID].analysisReady.bam'
exclude_flag1=0x100
exclude_flag2=0x4
expr_nm='[NM] < 12'
expr_sclen='sclen < 20'
filtered_bai='[BamDir]/[SeqID].recal.filtered.bam.bai'
filtered_bam='[BamDir]/[SeqID].recal.filtered.bam'
include_flag=0x2
min_mapq=20
CWD="."

# --- CLI Argument Parsing ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    --BamDir) BamDir="$2"; shift 2;;
    --SeqID) SeqID="$2"; shift 2;;
    --Threads) Threads="$2"; shift 2;;
    --analysis_ready_bai) analysis_ready_bai="$2"; shift 2;;
    --analysis_ready_bam) analysis_ready_bam="$2"; shift 2;;
    --exclude_flag1) exclude_flag1="$2"; shift 2;;
    --exclude_flag2) exclude_flag2="$2"; shift 2;;
    --expr_nm) expr_nm="$2"; shift 2;;
    --expr_sclen) expr_sclen="$2"; shift 2;;
    --filtered_bai) filtered_bai="$2"; shift 2;;
    --filtered_bam) filtered_bam="$2"; shift 2;;
    --include_flag) include_flag="$2"; shift 2;;
    --min_mapq) min_mapq="$2"; shift 2;;
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
    s="${s//\\[analysis_ready_bai\\]/${analysis_ready_bai}}"
    s="${s//\\[analysis_ready_bam\\]/${analysis_ready_bam}}"
    s="${s//\\[exclude_flag1\\]/${exclude_flag1}}"
    s="${s//\\[exclude_flag2\\]/${exclude_flag2}}"
    s="${s//\\[expr_nm\\]/${expr_nm}}"
    s="${s//\\[expr_sclen\\]/${expr_sclen}}"
    s="${s//\\[filtered_bai\\]/${filtered_bai}}"
    s="${s//\\[filtered_bam\\]/${filtered_bam}}"
    s="${s//\\[include_flag\\]/${include_flag}}"
    s="${s//\\[min_mapq\\]/${min_mapq}}"
  done
  # 앞뒤 공백 제거 및 연속 공백 축소
  echo "$s" | xargs
}

# --- Finalize Outputs & Command ---
filtered_bam=$(render "${filtered_bam}")
filtered_bai=$(render "${filtered_bai}")
analysis_ready_bam=$(render "${analysis_ready_bam}")
analysis_ready_bai=$(render "${analysis_ready_bai}")
CMD_LINE='samtools view -b -h -q [min_mapq] -f [include_flag] -F [exclude_flag1] -F [exclude_flag2] -e "[expr_sclen]" -e "[expr_nm]" --threads [Threads] [BamDir]/[SeqID].recal.bam > [BamDir]/[SeqID].recal.filtered.bam && samtools index --threads [Threads] [BamDir]/[SeqID].recal.filtered.bam && ln -Tsf  [BamDir]/[SeqID].recal.filtered.bam [BamDir]/[SeqID].analysisReady.bam && ln -Tsf  [BamDir]/[SeqID].recal.filtered.bam.bai [BamDir]/[SeqID].analysisReady.bam.bai'
CMD=$(render "$CMD_LINE")

echo -e "\n[RUNNING CMD]\n$CMD\n"

# --- Execution ---
mkdir -p "$CWD"
cd "$CWD" && eval "$CMD"
