#!/usr/bin/env bash
set -euo pipefail

# --- Defaults ---
BamDir=''
SeqID=''
Threads=8
bin=100000
bind=/storage,/data
global_dist_txt='[qcResDir]/[SeqID].global.dist.txt'
mapq=20
mosdepth_args='--no-per-base --fast-mode'
prefix='[qcResDir]/[SeqID]'
qcResDir=''
regions_bed_gz='[qcResDir]/[SeqID].regions.bed.gz'
regions_bed_gz_csi='[qcResDir]/[SeqID].regions.bed.gz.csi'
sif=/storage/images/mosdepth-0.3.6.sif
summary_txt='[qcResDir]/[SeqID].mosdepth.summary.txt'
CWD="."

# --- CLI Argument Parsing ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    --BamDir) BamDir="$2"; shift 2;;
    --SeqID) SeqID="$2"; shift 2;;
    --Threads) Threads="$2"; shift 2;;
    --bin) bin="$2"; shift 2;;
    --bind) bind="$2"; shift 2;;
    --global_dist_txt) global_dist_txt="$2"; shift 2;;
    --mapq) mapq="$2"; shift 2;;
    --mosdepth_args) mosdepth_args="$2"; shift 2;;
    --prefix) prefix="$2"; shift 2;;
    --qcResDir) qcResDir="$2"; shift 2;;
    --regions_bed_gz) regions_bed_gz="$2"; shift 2;;
    --regions_bed_gz_csi) regions_bed_gz_csi="$2"; shift 2;;
    --sif) sif="$2"; shift 2;;
    --summary_txt) summary_txt="$2"; shift 2;;
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
    s="${s//\\[bin\\]/${bin}}"
    s="${s//\\[bind\\]/${bind}}"
    s="${s//\\[global_dist_txt\\]/${global_dist_txt}}"
    s="${s//\\[mapq\\]/${mapq}}"
    s="${s//\\[mosdepth_args\\]/${mosdepth_args}}"
    s="${s//\\[prefix\\]/${prefix}}"
    s="${s//\\[qcResDir\\]/${qcResDir}}"
    s="${s//\\[regions_bed_gz\\]/${regions_bed_gz}}"
    s="${s//\\[regions_bed_gz_csi\\]/${regions_bed_gz_csi}}"
    s="${s//\\[sif\\]/${sif}}"
    s="${s//\\[summary_txt\\]/${summary_txt}}"
  done
  # 앞뒤 공백 제거 및 연속 공백 축소
  echo "$s" | xargs
}

# --- Finalize Outputs & Command ---
prefix=$(render "${prefix}")
regions_bed_gz=$(render "${regions_bed_gz}")
regions_bed_gz_csi=$(render "${regions_bed_gz_csi}")
summary_txt=$(render "${summary_txt}")
global_dist_txt=$(render "${global_dist_txt}")
CMD_LINE='singularity exec -B [bind] [sif] /opt/mosdepth --threads [Threads] [mosdepth_args] --by [bin] --mapq [mapq] [qcResDir]/[SeqID] [BamDir]/[SeqID].analysisReady.bam'
CMD=$(render "$CMD_LINE")

echo -e "\n[RUNNING CMD]\n$CMD\n"

# --- Execution ---
mkdir -p "$CWD"
cd "$CWD" && eval "$CMD"
