#!/usr/bin/env bash
set -euo pipefail

# Runner: mosdepth 0.3.6 (singularity_coverage_100kb)

# defaults (override by CLI)
SeqID=''
BamDir=''
qcResDir=''
bind=/storage,/data
sif=/storage/images/mosdepth-0.3.6.sif
Threads=8
by_bp=100000
mapq=20
mosdepth_args='--no-per-base --fast-mode'
prefix='[qcResDir]/[SeqID]'
regions_bed_gz='[qcResDir]/[SeqID].regions.bed.gz'
regions_bed_gz_csi='[qcResDir]/[SeqID].regions.bed.gz.csi'
summary_txt='[qcResDir]/[SeqID].mosdepth.summary.txt'
global_dist_txt='[qcResDir]/[SeqID].global.dist.txt'

CWD="."
PRINT_CMD=0
PRINT_OUTPUTS=0
EMIT_OUTPUTS=""

OUTPUT_KEYS=(prefix regions_bed_gz regions_bed_gz_csi summary_txt global_dist_txt)

usage() {
  echo "Usage: $0 --<Key> <val> ... [--cwd <dir>] [--print-cmd] [--print-outputs] [--emit-outputs <file>]" >&2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --SeqID) SeqID="$2"; shift 2;;
    --BamDir) BamDir="$2"; shift 2;;
    --qcResDir) qcResDir="$2"; shift 2;;
    --bind) bind="$2"; shift 2;;
    --sif) sif="$2"; shift 2;;
    --Threads) Threads="$2"; shift 2;;
    --by_bp) by_bp="$2"; shift 2;;
    --mapq) mapq="$2"; shift 2;;
    --mosdepth_args) mosdepth_args="$2"; shift 2;;
    --prefix) prefix="$2"; shift 2;;
    --regions_bed_gz) regions_bed_gz="$2"; shift 2;;
    --regions_bed_gz_csi) regions_bed_gz_csi="$2"; shift 2;;
    --summary_txt) summary_txt="$2"; shift 2;;
    --global_dist_txt) global_dist_txt="$2"; shift 2;;
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
[[ -n "${BamDir}" ]] || { echo "Missing required --BamDir" >&2; exit 2; }
[[ -n "${qcResDir}" ]] || { echo "Missing required --qcResDir" >&2; exit 2; }

render() {
  local s="$1"
  s="${s//[SeqID]/${SeqID}}"
  s="${s//[BamDir]/${BamDir}}"
  s="${s//[qcResDir]/${qcResDir}}"
  s="${s//[bind]/${bind}}"
  s="${s//[sif]/${sif}}"
  s="${s//[Threads]/${Threads}}"
  s="${s//[by_bp]/${by_bp}}"
  s="${s//[mapq]/${mapq}}"
  s="${s//[mosdepth_args]/${mosdepth_args}}"
  s="${s//[prefix]/${prefix}}"
  s="${s//[regions_bed_gz]/${regions_bed_gz}}"
  s="${s//[regions_bed_gz_csi]/${regions_bed_gz_csi}}"
  s="${s//[summary_txt]/${summary_txt}}"
  s="${s//[global_dist_txt]/${global_dist_txt}}"
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
if [[ -z "${prefix}" ]]; then
  __tmp='[qcResDir]/[SeqID]'
  prefix="$(render "$__tmp")"
fi
if [[ -z "${regions_bed_gz}" ]]; then
  __tmp='[qcResDir]/[SeqID].regions.bed.gz'
  regions_bed_gz="$(render "$__tmp")"
fi
if [[ -z "${regions_bed_gz_csi}" ]]; then
  __tmp='[qcResDir]/[SeqID].regions.bed.gz.csi'
  regions_bed_gz_csi="$(render "$__tmp")"
fi
if [[ -z "${summary_txt}" ]]; then
  __tmp='[qcResDir]/[SeqID].mosdepth.summary.txt'
  summary_txt="$(render "$__tmp")"
fi
if [[ -z "${global_dist_txt}" ]]; then
  __tmp='[qcResDir]/[SeqID].global.dist.txt'
  global_dist_txt="$(render "$__tmp")"
fi

CMD_LINE='singularity exec -B [bind] [sif] /opt/mosdepth --threads [Threads] [mosdepth_args] --by [by_bp] --mapq [mapq] [qcResDir]/[SeqID] [BamDir]/[SeqID].analysisReady.bam'
CMD="$(render "$CMD_LINE")"

if [[ "$PRINT_CMD" -eq 1 ]]; then
  echo "$CMD"
  exit 0
fi

mkdir -p "$CWD"
bash -lc "cd \"$CWD\" && $CMD"

emit_outputs
