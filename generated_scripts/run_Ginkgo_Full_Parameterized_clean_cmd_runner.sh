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
FinalDir='[ResultBaseDir]/[SeqID]/Binsize_[BinSize]kb'
Genome=''
GinkgoHomeDir=''
NormalDir=''
NormalSamples='Normal_01.bed.gz Normal_02.bed.gz Normal_03.bed.gz'
Normalization=''
ProcessR='[GinkgoHomeDir]/scripts/process.R'
R_Args_Base='status.xml data 2 [BinMeth] ward euclidean 1 refDummy.bed_mapped 0 ploidyDummy.txt 0 1'
R_Args_Reclust='status.xml [BinMeth] ward euclidean 0 ploidyDummy.txt 0'
ReadLength=''
ReclustR='[GinkgoHomeDir]/scripts/reclust.R'
ResultBaseDir=''
SampleName=''
SeqID=''
WorkDir='[GinkgoHomeDir]/uploads/[SeqID]_[BinSize]kb'
f=''
n=''
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
    --FinalDir) FinalDir="$2"; shift 2;;
    --Genome) Genome="$2"; shift 2;;
    --GinkgoHomeDir) GinkgoHomeDir="$2"; shift 2;;
    --NormalDir) NormalDir="$2"; shift 2;;
    --NormalSamples) NormalSamples="$2"; shift 2;;
    --Normalization) Normalization="$2"; shift 2;;
    --ProcessR) ProcessR="$2"; shift 2;;
    --R_Args_Base) R_Args_Base="$2"; shift 2;;
    --R_Args_Reclust) R_Args_Reclust="$2"; shift 2;;
    --ReadLength) ReadLength="$2"; shift 2;;
    --ReclustR) ReclustR="$2"; shift 2;;
    --ResultBaseDir) ResultBaseDir="$2"; shift 2;;
    --SampleName) SampleName="$2"; shift 2;;
    --SeqID) SeqID="$2"; shift 2;;
    --WorkDir) WorkDir="$2"; shift 2;;
    --f) f="$2"; shift 2;;
    --n) n="$2"; shift 2;;
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
    s="${s//\\[FinalDir\\]/${FinalDir}}"
    s="${s//\\[Genome\\]/${Genome}}"
    s="${s//\\[GinkgoHomeDir\\]/${GinkgoHomeDir}}"
    s="${s//\\[NormalDir\\]/${NormalDir}}"
    s="${s//\\[NormalSamples\\]/${NormalSamples}}"
    s="${s//\\[Normalization\\]/${Normalization}}"
    s="${s//\\[ProcessR\\]/${ProcessR}}"
    s="${s//\\[R_Args_Base\\]/${R_Args_Base}}"
    s="${s//\\[R_Args_Reclust\\]/${R_Args_Reclust}}"
    s="${s//\\[ReadLength\\]/${ReadLength}}"
    s="${s//\\[ReclustR\\]/${ReclustR}}"
    s="${s//\\[ResultBaseDir\\]/${ResultBaseDir}}"
    s="${s//\\[SampleName\\]/${SampleName}}"
    s="${s//\\[SeqID\\]/${SeqID}}"
    s="${s//\\[WorkDir\\]/${WorkDir}}"
    s="${s//\\[f\\]/${f}}"
    s="${s//\\[n\\]/${n}}"
  done
  # 앞뒤 공백 제거 및 연속 공백 축소
  echo "$s" | xargs
}

# --- Finalize Outputs & Command ---
BinMeth=$(render "${BinMeth}")
WorkDir=$(render "${WorkDir}")
FinalDir=$(render "${FinalDir}")
BinUnsorted=$(render "${BinUnsorted}")
BinFile=$(render "${BinFile}")
BinCount=$(render "${BinCount}")
ProcessR=$(render "${ProcessR}")
ReclustR=$(render "${ReclustR}")
CNVCaller=$(render "${CNVCaller}")
CMD_LINE='for n in [NormalSamples]; do ln -Tsf [NormalDir]/[n] [WorkDir]/[n]; done && ln -Tsf [BamToBedDir]/[SeqID].bed.gz [WorkDir]/[SeqID].bed.gz && ls [WorkDir] | grep ".bed.gz$" > [WorkDir]/list &&
while read -r f; do
  echo ">> Binning: [f]" &&
  SampleName="[f%.bed.gz]" &&
  SampleName="[SampleName%.bed]" &&
  
  # [Normalization] chr 접두사 보정
  if [[ "[f]" =~ \.gz$ ]]; then
    zcat [WorkDir]/[f] | awk '"'"'{if($1!~/^chr/) print "chr"$0; else print $0]'"'"' | gzip > [WorkDir]/[f]_tmp.gz &&
    [BinUnsorted] [BinFile] [BinCount] <(zcat [WorkDir]/[f]_tmp.gz) [SampleName] [WorkDir]/[f]_mapped;
  else
    awk '"'"'{if($1!~/^chr/) print "chr"$0; else print $0]'"'"' [WorkDir]/[f] > [WorkDir]/[f]_tmp &&
    [BinUnsorted] [BinFile] [BinCount] [WorkDir]/[f]_tmp [SampleName] [WorkDir]/[f]_mapped;
  fi;
done < [WorkDir]/list &&
# 데이터 통합 및 분석 실행 paste [WorkDir]/*_mapped > [WorkDir]/data && cd [GinkgoHomeDir] && [ProcessR] [GinkgoHomeDir]/genomes/[Genome]/original [WorkDir] [R_Args_Base] && [ReclustR] [GinkgoHomeDir]/genomes/[Genome]/original [WorkDir] [R_Args_Reclust] && [CNVCaller] [WorkDir]/SegCopy [WorkDir]/CNV1 [WorkDir]/CNV2 &&
# 결과 정리 mkdir -p [FinalDir] && mv [WorkDir]/* [FinalDir]/ && rm -rf [WorkDir]'
CMD=$(render "$CMD_LINE")

echo -e "\n[RUNNING CMD]\n$CMD\n"

# --- Execution ---
mkdir -p "$CWD"
cd "$CWD" && eval "$CMD"
