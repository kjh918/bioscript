#!/usr/bin/env bash
set -euo pipefail

# Runner: CopyKitAnalysis 1.0.0 (cbNIPT_CopyNumber_Analysis)

# defaults (override by CLI)
SeqID=''
SampleID=''
SamplePloidy=''
NGS_DataBaseDir=''
ResultBaseDir=''
threads=2
BinSizeList='110kb 220kb'
GenomeVersion=hg38
Rscript_path=/storage/home/kangsm/runScripts/cbNIPT/GCX.cbNIPT_Run.CopyKit.Analysis.R
AnalysisRunDir='[ResultBaseDir]/[SampleID]'

CWD="."
PRINT_CMD=0
PRINT_OUTPUTS=0
EMIT_OUTPUTS=""

OUTPUT_KEYS=(AnalysisRunDir)

usage() {
  echo "Usage: $0 --<Key> <val> ... [--cwd <dir>] [--print-cmd] [--print-outputs] [--emit-outputs <file>]" >&2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --SeqID) SeqID="$2"; shift 2;;
    --SampleID) SampleID="$2"; shift 2;;
    --SamplePloidy) SamplePloidy="$2"; shift 2;;
    --NGS_DataBaseDir) NGS_DataBaseDir="$2"; shift 2;;
    --ResultBaseDir) ResultBaseDir="$2"; shift 2;;
    --threads) threads="$2"; shift 2;;
    --BinSizeList) BinSizeList="$2"; shift 2;;
    --GenomeVersion) GenomeVersion="$2"; shift 2;;
    --Rscript_path) Rscript_path="$2"; shift 2;;
    --AnalysisRunDir) AnalysisRunDir="$2"; shift 2;;
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
[[ -n "${SampleID}" ]] || { echo "Missing required --SampleID" >&2; exit 2; }
[[ -n "${SamplePloidy}" ]] || { echo "Missing required --SamplePloidy" >&2; exit 2; }
[[ -n "${ResultBaseDir}" ]] || { echo "Missing required --ResultBaseDir" >&2; exit 2; }

render() {
  local s="$1"
  s="${s//[SeqID]/${SeqID}}"
  s="${s//[SampleID]/${SampleID}}"
  s="${s//[SamplePloidy]/${SamplePloidy}}"
  s="${s//[NGS_DataBaseDir]/${NGS_DataBaseDir}}"
  s="${s//[ResultBaseDir]/${ResultBaseDir}}"
  s="${s//[threads]/${threads}}"
  s="${s//[BinSizeList]/${BinSizeList}}"
  s="${s//[GenomeVersion]/${GenomeVersion}}"
  s="${s//[Rscript_path]/${Rscript_path}}"
  s="${s//[AnalysisRunDir]/${AnalysisRunDir}}"
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
if [[ -z "${AnalysisRunDir}" ]]; then
  __tmp='[ResultBaseDir]/[SampleID]'
  AnalysisRunDir="$(render "$__tmp")"
fi

CMD_LINE='mkdir -p [AnalysisRunDir] && ln -Tsf [NGS_DataBaseDir]/[SeqID]/bam/[SeqID].analysisReady.bam [AnalysisRunDir]/[SampleID].bam && ln -Tsf [NGS_DataBaseDir]/[SeqID]/bam/[SeqID].analysisReady.bam.bai [AnalysisRunDir]/[SampleID].bam.bai && parallel -j [threads] -k "Rscript [Rscript_path]  --SeqID [SampleID]  --AnalysisRunDir [AnalysisRunDir]  --BinSize {}  --Ploidy [SamplePloidy]  --GenomeVersion [GenomeVersion]" ::: [BinSizeList]'
CMD="$(render "$CMD_LINE")"

if [[ "$PRINT_CMD" -eq 1 ]]; then
  echo "$CMD"
  exit 0
fi

mkdir -p "$CWD"
bash -lc "cd \"$CWD\" && $CMD"

emit_outputs
