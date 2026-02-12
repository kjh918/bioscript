#!/usr/bin/env bash
set -euo pipefail

# --- Defaults ---
BamToBedDir=''
BinMeth='variable_[BinSize]000_[ReadLength]_bwa'
BinSize=''
BinUnsorted='[GinkgoHomeDir]/scripts/binUnsorted'
CNVCaller='[GinkgoHomeDir]/scripts/CNVcaller'
ClustMeth=''
DistMeth=''
FinalDir='[ResultBaseDir]/[SeqID]/Binsize_[BinSize]kb'
Genome=''
GenomeDir='[GinkgoHomeDir]/genomes/[Genome]/original'
GinkgoHomeDir=''
NormalDir=''
NormalSamples='Normal_01.bed.gz Normal_02.bed.gz Normal_03.bed.gz'
ProcessR='[GinkgoHomeDir]/scripts/process.R'
ReadLength=''
ReclustR='[GinkgoHomeDir]/scripts/reclust.R'
ResultBaseDir=''
RunID='[SeqID]_[BinSize]kb'
SeqID=''
WorkDir='[GinkgoHomeDir]/uploads/[SeqID]_[BinSize]kb'
file=''
normal=''
CWD="."

# --- CLI Argument Parsing ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    --BamToBedDir) BamToBedDir="$2"; shift 2;;
    --BinMeth) BinMeth="$2"; shift 2;;
    --BinSize) BinSize="$2"; shift 2;;
    --BinUnsorted) BinUnsorted="$2"; shift 2;;
    --CNVCaller) CNVCaller="$2"; shift 2;;
    --ClustMeth) ClustMeth="$2"; shift 2;;
    --DistMeth) DistMeth="$2"; shift 2;;
    --FinalDir) FinalDir="$2"; shift 2;;
    --Genome) Genome="$2"; shift 2;;
    --GenomeDir) GenomeDir="$2"; shift 2;;
    --GinkgoHomeDir) GinkgoHomeDir="$2"; shift 2;;
    --NormalDir) NormalDir="$2"; shift 2;;
    --NormalSamples) NormalSamples="$2"; shift 2;;
    --ProcessR) ProcessR="$2"; shift 2;;
    --ReadLength) ReadLength="$2"; shift 2;;
    --ReclustR) ReclustR="$2"; shift 2;;
    --ResultBaseDir) ResultBaseDir="$2"; shift 2;;
    --RunID) RunID="$2"; shift 2;;
    --SeqID) SeqID="$2"; shift 2;;
    --WorkDir) WorkDir="$2"; shift 2;;
    --file) file="$2"; shift 2;;
    --normal) normal="$2"; shift 2;;
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
    s="${s//\\[BinMeth\\]/${BinMeth}}"
    s="${s//\\[BinSize\\]/${BinSize}}"
    s="${s//\\[BinUnsorted\\]/${BinUnsorted}}"
    s="${s//\\[CNVCaller\\]/${CNVCaller}}"
    s="${s//\\[ClustMeth\\]/${ClustMeth}}"
    s="${s//\\[DistMeth\\]/${DistMeth}}"
    s="${s//\\[FinalDir\\]/${FinalDir}}"
    s="${s//\\[Genome\\]/${Genome}}"
    s="${s//\\[GenomeDir\\]/${GenomeDir}}"
    s="${s//\\[GinkgoHomeDir\\]/${GinkgoHomeDir}}"
    s="${s//\\[NormalDir\\]/${NormalDir}}"
    s="${s//\\[NormalSamples\\]/${NormalSamples}}"
    s="${s//\\[ProcessR\\]/${ProcessR}}"
    s="${s//\\[ReadLength\\]/${ReadLength}}"
    s="${s//\\[ReclustR\\]/${ReclustR}}"
    s="${s//\\[ResultBaseDir\\]/${ResultBaseDir}}"
    s="${s//\\[RunID\\]/${RunID}}"
    s="${s//\\[SeqID\\]/${SeqID}}"
    s="${s//\\[WorkDir\\]/${WorkDir}}"
    s="${s//\\[file\\]/${file}}"
    s="${s//\\[normal\\]/${normal}}"
  done
  # 앞뒤 공백 제거 및 연속 공백 축소
  echo "$s" | xargs
}

# --- Finalize Outputs & Command ---
RunID=$(render "${RunID}")
WorkDir=$(render "${WorkDir}")
FinalDir=$(render "${FinalDir}")
BinUnsorted=$(render "${BinUnsorted}")
ProcessR=$(render "${ProcessR}")
ReclustR=$(render "${ReclustR}")
CNVCaller=$(render "${CNVCaller}")
GenomeDir=$(render "${GenomeDir}")
CMD_LINE='echo "==== [Integrated Ginkgo] Starting [RunID] ====" && mkdir -p [WorkDir] &&
# 1. 데이터 준비 및 링크 for normal in [NormalSamples]; do ln -Tsf [NormalDir]/[normal] [WorkDir]/[normal]; done && ln -Tsf [BamToBedDir]/[SeqID].bed.gz [WorkDir]/[SeqID].bed.gz && ls [WorkDir] | grep ".bed.gz$" > [WorkDir]/list &&
# 2. 전처리 (analyze.sh 로직 재현: chr 접두사 체크 및 Binning) while read file; do \
  echo ">> Binning: [file]"; \
  if [[ "[file]" =~ \.gz$ ]]; then \
    zcat [WorkDir]/[file] | awk '"'"'{if($1!~/^chr/) print "chr"$0; else print $0]'"'"' | gzip > [WorkDir]/[file]_tmp.gz; \
    [BinUnsorted] [GenomeDir]/[BinMeth] $(wc -l < [GenomeDir]/[BinMeth]) <(zcat [WorkDir]/[file]_tmp.gz) [file%.bed.gz] [WorkDir]/[file]_mapped; \
  else \
    awk '"'"'{if($1!~/^chr/) print "chr"$0; else print $0]'"'"' [WorkDir]/[file] > [WorkDir]/[file]_tmp; \
    [BinUnsorted] [GenomeDir]/[BinMeth] $(wc -l < [GenomeDir]/[BinMeth]) [WorkDir]/[file]_tmp [file%.bed] [WorkDir]/[file]_mapped; \
  fi; \
done < [WorkDir]/list &&
# 3. 데이터 통합 (Merge binned reads) paste [WorkDir]/*_mapped > [WorkDir]/data &&
# 4. R 분석 실행 (Ginkgo Core) cd [GinkgoHomeDir] && \ [ProcessR] [GenomeDir] [WorkDir] status.xml data 2 [BinMeth] [ClustMeth] [DistMeth] 1 refDummy.bed_mapped 0 ploidyDummy.txt 0 1 && \ [ReclustR] [GenomeDir] [WorkDir] status.xml [BinMeth] [ClustMeth] [DistMeth] 0 ploidyDummy.txt 0 &&
# 5. CNV Calling (analyze.sh 후반부 로직) [CNVCaller] [WorkDir]/SegCopy [WorkDir]/CNV1 [WorkDir]/CNV2 &&
# 6. 결과 이동 및 청소 mkdir -p [FinalDir] && \ mv [WorkDir]/* [FinalDir]/ && \ rm -rf [WorkDir] && \ echo "==== [Integrated Ginkgo] Completed [RunID] ===="'
CMD=$(render "$CMD_LINE")

echo -e "\n[RUNNING CMD]\n$CMD\n"

# --- Execution ---
mkdir -p "$CWD"
cd "$CWD" && eval "$CMD"
