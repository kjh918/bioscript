#!/usr/bin/env bash
set -euo pipefail

# --- Defaults ---
BamToBedDir=''
BinCount='$(wc -l < [BinFile])'
BinFile='[GinkgoHomeDir]/genomes/[Genome]/original/[BinMeth]'
BinMeth='variable_[BinSize]000_[ReadLength]_bwa'
BinSize=''
BinUnsorted='[GinkgoHomeDir]/scripts/binUnsorted'
CNVCaller='[GinkgoHomeDir]/scripts/CNVcaller'
Genome=''
GinkgoHomeDir=''
NormalDir=''
NormalSamples='Normal_01.bed.gz Normal_02.bed.gz Normal_03.bed.gz'
ProcessR='[GinkgoHomeDir]/scripts/process.R'
R_Args_Base='status.xml data 2 [BinMeth] ward euclidean 1 refDummy.bed_mapped 0 ploidyDummy.txt 0 1'
R_Args_Reclust='status.xml [BinMeth] ward euclidean 0 ploidyDummy.txt 0'
ReadLength=''
ReclustR='[GinkgoHomeDir]/scripts/reclust.R'
ResultBaseDir=''
SeqID=''
Threads=''
WorkDir='[ResultBaseDir]/[SeqID]/[SeqID]_[BinSize]kb'
CWD="."

# --- CLI Argument Parsing ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    --BamToBedDir) BamToBedDir="$2"; shift 2;;
    --BinCount) BinCount="$2"; shift 2;;
    --BinFile) BinFile="$2"; shift 2;;
    --BinMeth) BinMeth="$2"; shift 2;;
    --BinSize) BinSize="$2"; shift 2;;
    --BinUnsorted) BinUnsorted="$2"; shift 2;;
    --CNVCaller) CNVCaller="$2"; shift 2;;
    --Genome) Genome="$2"; shift 2;;
    --GinkgoHomeDir) GinkgoHomeDir="$2"; shift 2;;
    --NormalDir) NormalDir="$2"; shift 2;;
    --NormalSamples) NormalSamples="$2"; shift 2;;
    --ProcessR) ProcessR="$2"; shift 2;;
    --R_Args_Base) R_Args_Base="$2"; shift 2;;
    --R_Args_Reclust) R_Args_Reclust="$2"; shift 2;;
    --ReadLength) ReadLength="$2"; shift 2;;
    --ReclustR) ReclustR="$2"; shift 2;;
    --ResultBaseDir) ResultBaseDir="$2"; shift 2;;
    --SeqID) SeqID="$2"; shift 2;;
    --Threads) Threads="$2"; shift 2;;
    --WorkDir) WorkDir="$2"; shift 2;;
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
    s="${s//\\[BamToBedDir\\]/${BamToBedDir}}"
    s="${s//\\[BinCount\\]/${BinCount}}"
    s="${s//\\[BinFile\\]/${BinFile}}"
    s="${s//\\[BinMeth\\]/${BinMeth}}"
    s="${s//\\[BinSize\\]/${BinSize}}"
    s="${s//\\[BinUnsorted\\]/${BinUnsorted}}"
    s="${s//\\[CNVCaller\\]/${CNVCaller}}"
    s="${s//\\[Genome\\]/${Genome}}"
    s="${s//\\[GinkgoHomeDir\\]/${GinkgoHomeDir}}"
    s="${s//\\[NormalDir\\]/${NormalDir}}"
    s="${s//\\[NormalSamples\\]/${NormalSamples}}"
    s="${s//\\[ProcessR\\]/${ProcessR}}"
    s="${s//\\[R_Args_Base\\]/${R_Args_Base}}"
    s="${s//\\[R_Args_Reclust\\]/${R_Args_Reclust}}"
    s="${s//\\[ReadLength\\]/${ReadLength}}"
    s="${s//\\[ReclustR\\]/${ReclustR}}"
    s="${s//\\[ResultBaseDir\\]/${ResultBaseDir}}"
    s="${s//\\[SeqID\\]/${SeqID}}"
    s="${s//\\[Threads\\]/${Threads}}"
    s="${s//\\[WorkDir\\]/${WorkDir}}"
  done
  # 앞뒤 공백 제거 및 연속 공백 축소
  echo "$s" | xargs
}

# --- Finalize Outputs & Command ---
BinMeth=$(render "${BinMeth}")
WorkDir=$(render "${WorkDir}")
BinUnsorted=$(render "${BinUnsorted}")
BinFile=$(render "${BinFile}")
BinCount=$(render "${BinCount}")
ProcessR=$(render "${ProcessR}")
ReclustR=$(render "${ReclustR}")
CNVCaller=$(render "${CNVCaller}")
CMD_LINE='mkdir -p [WorkDir] && ln -sf [NormalDir]/*.bed.gz [WorkDir]/ && ln -sf [BamToBedDir]/[SeqID].bed.gz [WorkDir]/ &&
find [WorkDir] -maxdepth 1 -name "*.gz" | xargs -I {] bash -c '"'"'N=$(basename "{]" | sed "s/\.bed.*//"); [BinUnsorted] [BinFile] [BinCount] <(zcat "{]" | awk "{if(\$1!~/^chr/) print \"chr\"\$0; else print \$0]") $N "{]"_mapped'"'"' &&
find [WorkDir] -maxdepth 1 -name "*.bed" ! -name "*.gz" | xargs -I {] bash -c '"'"'N=$(basename "{]" | sed "s/\.bed.*//"); [BinUnsorted] [BinFile] [BinCount] <(awk "{if(\$1!~/^chr/) print \"chr\"\$0; else print \$0]" "{]") $N "{]"_mapped'"'"' &&
paste [WorkDir]/*_mapped > [WorkDir]/data &&
[ProcessR] [GinkgoHomeDir]/genomes/[Genome]/original [WorkDir] [R_Args_Base] && [ReclustR] [GinkgoHomeDir]/genomes/[Genome]/original [WorkDir] [R_Args_Reclust] && [CNVCaller] [WorkDir]/SegCopy [WorkDir]/CNV1 [WorkDir]/CNV2'
CMD=$(render "$CMD_LINE")

echo -e "\n[RUNNING CMD]\n$CMD\n"

# --- Execution ---
mkdir -p "$CWD"
cd "$CWD" && eval "$CMD"
