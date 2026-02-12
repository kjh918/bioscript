#!/usr/bin/env bash
set -euo pipefail

# --- Defaults ---
AnalysisRunDir='[ResultBaseDir]/[SampleID]'
BinSize=220kb
GenomeVersion=hg38
NGS_DataBaseDir=''
Rscript_path=/storage/home/kangsm/runScripts/cbNIPT/GCX.cbNIPT_Run.CopyKit.Analysis.R
SampleID=''
SamplePloidy=2
SeqID=''
Threads=1
CWD="."

# --- CLI Argument Parsing ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    --AnalysisRunDir) AnalysisRunDir="$2"; shift 2;;
    --BinSize) BinSize="$2"; shift 2;;
    --GenomeVersion) GenomeVersion="$2"; shift 2;;
    --NGS_DataBaseDir) NGS_DataBaseDir="$2"; shift 2;;
    --Rscript_path) Rscript_path="$2"; shift 2;;
    --SampleID) SampleID="$2"; shift 2;;
    --SamplePloidy) SamplePloidy="$2"; shift 2;;
    --SeqID) SeqID="$2"; shift 2;;
    --Threads) Threads="$2"; shift 2;;
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
    s="${s//\\[AnalysisRunDir\\]/${AnalysisRunDir}}"
    s="${s//\\[BinSize\\]/${BinSize}}"
    s="${s//\\[GenomeVersion\\]/${GenomeVersion}}"
    s="${s//\\[NGS_DataBaseDir\\]/${NGS_DataBaseDir}}"
    s="${s//\\[Rscript_path\\]/${Rscript_path}}"
    s="${s//\\[SampleID\\]/${SampleID}}"
    s="${s//\\[SamplePloidy\\]/${SamplePloidy}}"
    s="${s//\\[SeqID\\]/${SeqID}}"
    s="${s//\\[Threads\\]/${Threads}}"
  done
  # 앞뒤 공백 제거 및 연속 공백 축소
  echo "$s" | xargs
}

# --- Finalize Outputs & Command ---
AnalysisRunDir=$(render "${AnalysisRunDir}")
CMD_LINE='mkdir -p [AnalysisRunDir] && ln -Tsf [NGS_DataBaseDir]/[SeqID].analysisReady.bam [AnalysisRunDir]/[SampleID].bam && ln -Tsf [NGS_DataBaseDir]/[SeqID].analysisReady.bam.bai [AnalysisRunDir]/[SampleID].bam.bai && Rscript [Rscript_path]  --SeqID [SampleID]  --AnalysisRunDir [AnalysisRunDir]  --BinSize [BinSize]  --Ploidy [SamplePloidy]  --GenomeVersion [GenomeVersion]'
CMD=$(render "$CMD_LINE")

echo -e "\n[RUNNING CMD]\n$CMD\n"

# --- Execution ---
mkdir -p "$CWD"
cd "$CWD" && eval "$CMD"
