#!/usr/bin/env bash
set -euo pipefail

# --- Defaults ---
BamDir=''
BwaIndex=''
ReadGroupCenter=''
ReadGroupID=''
ReadGroupLibrary=''
ReadGroupPlatform=''
SeqID=''
Threads=8
TrimFastqDir=''
bind=/storage,/data
bwa_args='-M -Y -L 50,50'
primary_bam='[BamDir]/[SeqID].primary.bam'
sif=/storage/images/bwa-0.7.17.sif
CWD="."

# --- CLI Argument Parsing ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    --BamDir) BamDir="$2"; shift 2;;
    --BwaIndex) BwaIndex="$2"; shift 2;;
    --ReadGroupCenter) ReadGroupCenter="$2"; shift 2;;
    --ReadGroupID) ReadGroupID="$2"; shift 2;;
    --ReadGroupLibrary) ReadGroupLibrary="$2"; shift 2;;
    --ReadGroupPlatform) ReadGroupPlatform="$2"; shift 2;;
    --SeqID) SeqID="$2"; shift 2;;
    --Threads) Threads="$2"; shift 2;;
    --TrimFastqDir) TrimFastqDir="$2"; shift 2;;
    --bind) bind="$2"; shift 2;;
    --bwa_args) bwa_args="$2"; shift 2;;
    --primary_bam) primary_bam="$2"; shift 2;;
    --sif) sif="$2"; shift 2;;
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
    s="${s//\\[BwaIndex\\]/${BwaIndex}}"
    s="${s//\\[ReadGroupCenter\\]/${ReadGroupCenter}}"
    s="${s//\\[ReadGroupID\\]/${ReadGroupID}}"
    s="${s//\\[ReadGroupLibrary\\]/${ReadGroupLibrary}}"
    s="${s//\\[ReadGroupPlatform\\]/${ReadGroupPlatform}}"
    s="${s//\\[SeqID\\]/${SeqID}}"
    s="${s//\\[Threads\\]/${Threads}}"
    s="${s//\\[TrimFastqDir\\]/${TrimFastqDir}}"
    s="${s//\\[bind\\]/${bind}}"
    s="${s//\\[bwa_args\\]/${bwa_args}}"
    s="${s//\\[primary_bam\\]/${primary_bam}}"
    s="${s//\\[sif\\]/${sif}}"
  done
  # 앞뒤 공백 제거 및 연속 공백 축소
  echo "$s" | xargs
}

# --- Finalize Outputs & Command ---
primary_bam=$(render "${primary_bam}")
CMD_LINE='singularity exec -B [bind] [sif] bwa mem [bwa_args] -t [Threads] -R "@RG\tID:[ReadGroupID]\tPL:[ReadGroupPlatform]\tLB:[ReadGroupLibrary]\tSM:[SeqID]\tCN:[ReadGroupCenter]" [BwaIndex] [TrimFastqDir]/[SeqID].trimmed_R1.fastq.gz [TrimFastqDir]/[SeqID].trimmed_R2.fastq.gz | samtools view -bS -o [BamDir]/[SeqID].primary.bam -'
CMD=$(render "$CMD_LINE")

echo -e "\n[RUNNING CMD]\n$CMD\n"

# --- Execution ---
mkdir -p "$CWD"
cd "$CWD" && eval "$CMD"
