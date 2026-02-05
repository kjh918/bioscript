#!/usr/bin/env bash
set -euo pipefail

# --- Defaults ---
BamDir=''
KnownIndel1=''
KnownIndel2=''
KnownSnp=''
RecalInputBAM=''
ReferenceFasta=''
SeqID=''
bind=/storage,/data
gc_threads=14
qcResDir=''
recal_bai='[BamDir]/[SeqID].recal.bam.bai'
recal_bam='[BamDir]/[SeqID].recal.bam'
recal_table='[qcResDir]/[SeqID].recal.table.txt'
sif=/storage/images/gatk-4.4.0.0.sif
xmx_mb=16384
CWD="."

# --- CLI Argument Parsing ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    --BamDir) BamDir="$2"; shift 2;;
    --KnownIndel1) KnownIndel1="$2"; shift 2;;
    --KnownIndel2) KnownIndel2="$2"; shift 2;;
    --KnownSnp) KnownSnp="$2"; shift 2;;
    --RecalInputBAM) RecalInputBAM="$2"; shift 2;;
    --ReferenceFasta) ReferenceFasta="$2"; shift 2;;
    --SeqID) SeqID="$2"; shift 2;;
    --bind) bind="$2"; shift 2;;
    --gc_threads) gc_threads="$2"; shift 2;;
    --qcResDir) qcResDir="$2"; shift 2;;
    --recal_bai) recal_bai="$2"; shift 2;;
    --recal_bam) recal_bam="$2"; shift 2;;
    --recal_table) recal_table="$2"; shift 2;;
    --sif) sif="$2"; shift 2;;
    --xmx_mb) xmx_mb="$2"; shift 2;;
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
    s="${s//\\[KnownIndel1\\]/${KnownIndel1}}"
    s="${s//\\[KnownIndel2\\]/${KnownIndel2}}"
    s="${s//\\[KnownSnp\\]/${KnownSnp}}"
    s="${s//\\[RecalInputBAM\\]/${RecalInputBAM}}"
    s="${s//\\[ReferenceFasta\\]/${ReferenceFasta}}"
    s="${s//\\[SeqID\\]/${SeqID}}"
    s="${s//\\[bind\\]/${bind}}"
    s="${s//\\[gc_threads\\]/${gc_threads}}"
    s="${s//\\[qcResDir\\]/${qcResDir}}"
    s="${s//\\[recal_bai\\]/${recal_bai}}"
    s="${s//\\[recal_bam\\]/${recal_bam}}"
    s="${s//\\[recal_table\\]/${recal_table}}"
    s="${s//\\[sif\\]/${sif}}"
    s="${s//\\[xmx_mb\\]/${xmx_mb}}"
  done
  # 앞뒤 공백 제거 및 연속 공백 축소
  echo "$s" | xargs
}

# --- Finalize Outputs & Command ---
recal_table=$(render "${recal_table}")
recal_bam=$(render "${recal_bam}")
recal_bai=$(render "${recal_bai}")
CMD_LINE='singularity exec -B [bind] [sif] gatk BaseRecalibrator --java-options "-XX:ParallelGCThreads=[gc_threads] -Xmx[xmx_mb]m" --input [RecalInputBAM] --reference [ReferenceFasta] --output [qcResDir]/[SeqID].recal.table.txt --known-sites [KnownSnp] --known-sites [KnownIndel1] --known-sites [KnownIndel2] && singularity exec -B [bind] [sif] gatk ApplyBQSR --java-options "-XX:ParallelGCThreads=[gc_threads] -Xmx[xmx_mb]m" --input [RecalInputBAM] --bqsr-recal-file [qcResDir]/[SeqID].recal.table.txt --output [BamDir]/[SeqID].recal.bam'
CMD=$(render "$CMD_LINE")

echo -e "\n[RUNNING CMD]\n$CMD\n"

# --- Execution ---
mkdir -p "$CWD"
cd "$CWD" && eval "$CMD"
