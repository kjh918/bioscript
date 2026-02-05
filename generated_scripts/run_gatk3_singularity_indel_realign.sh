#!/usr/bin/env bash
set -euo pipefail

# --- Defaults ---
BamDir=''
KnownIndel1=''
KnownIndel2=''
ReferenceFasta=''
SeqID=''
Threads=8
bind=/storage,/data
qcResDir=''
realigned_bam='[BamDir]/[SeqID].sorted.dedup.realign.bam'
sif=/storage/images/gatk-3.8-1.sif
target_intervals='[qcResDir]/[SeqID].realign.target.intervals'
xmx_mb=16384
CWD="."

# --- CLI Argument Parsing ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    --BamDir) BamDir="$2"; shift 2;;
    --KnownIndel1) KnownIndel1="$2"; shift 2;;
    --KnownIndel2) KnownIndel2="$2"; shift 2;;
    --ReferenceFasta) ReferenceFasta="$2"; shift 2;;
    --SeqID) SeqID="$2"; shift 2;;
    --Threads) Threads="$2"; shift 2;;
    --bind) bind="$2"; shift 2;;
    --qcResDir) qcResDir="$2"; shift 2;;
    --realigned_bam) realigned_bam="$2"; shift 2;;
    --sif) sif="$2"; shift 2;;
    --target_intervals) target_intervals="$2"; shift 2;;
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
    s="${s//\\[ReferenceFasta\\]/${ReferenceFasta}}"
    s="${s//\\[SeqID\\]/${SeqID}}"
    s="${s//\\[Threads\\]/${Threads}}"
    s="${s//\\[bind\\]/${bind}}"
    s="${s//\\[qcResDir\\]/${qcResDir}}"
    s="${s//\\[realigned_bam\\]/${realigned_bam}}"
    s="${s//\\[sif\\]/${sif}}"
    s="${s//\\[target_intervals\\]/${target_intervals}}"
    s="${s//\\[xmx_mb\\]/${xmx_mb}}"
  done
  # 앞뒤 공백 제거 및 연속 공백 축소
  echo "$s" | xargs
}

# --- Finalize Outputs & Command ---
target_intervals=$(render "${target_intervals}")
realigned_bam=$(render "${realigned_bam}")
CMD_LINE='singularity exec -B [bind] [sif] java -Xmx[xmx_mb]m -jar /usr/GenomeAnalysisTK.jar -T RealignerTargetCreator -R [ReferenceFasta] -I [BamDir]/[SeqID].sorted.dedup.bam -o [qcResDir]/[SeqID].realign.target.intervals -known [KnownIndel1] -known [KnownIndel2] -nt [Threads] && singularity exec -B [bind] [sif] java -Xmx[xmx_mb]m -jar /usr/GenomeAnalysisTK.jar -T IndelRealigner -R [ReferenceFasta] -targetIntervals [qcResDir]/[SeqID].realign.target.intervals -known [KnownIndel1] -known [KnownIndel2] -I [BamDir]/[SeqID].sorted.dedup.bam -o [BamDir]/[SeqID].sorted.dedup.realign.bam'
CMD=$(render "$CMD_LINE")

echo -e "\n[RUNNING CMD]\n$CMD\n"

# --- Execution ---
mkdir -p "$CWD"
cd "$CWD" && eval "$CMD"
