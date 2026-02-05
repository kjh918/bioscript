#!/usr/bin/env bash
set -euo pipefail

# --- Defaults ---
BinSize=''
BinSizeList='100 150 200 250'
FinalDir=''
FinalResultRoot='[ResultBaseDir]/[SeqID]'
GinkgoScript=/storage/home/kangsm/runScripts/cbNIPT/GCX.cbNIPT_Run.Ginkgo.Analysis.sh
ReadLength=''
ResultBaseDir=''
SeqID=''
threads=1
CWD="."

# --- CLI Argument Parsing ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    --BinSize) BinSize="$2"; shift 2;;
    --BinSizeList) BinSizeList="$2"; shift 2;;
    --FinalDir) FinalDir="$2"; shift 2;;
    --FinalResultRoot) FinalResultRoot="$2"; shift 2;;
    --GinkgoScript) GinkgoScript="$2"; shift 2;;
    --ReadLength) ReadLength="$2"; shift 2;;
    --ResultBaseDir) ResultBaseDir="$2"; shift 2;;
    --SeqID) SeqID="$2"; shift 2;;
    --threads) threads="$2"; shift 2;;
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
    s="${s//\\[BinSize\\]/${BinSize}}"
    s="${s//\\[BinSizeList\\]/${BinSizeList}}"
    s="${s//\\[FinalDir\\]/${FinalDir}}"
    s="${s//\\[FinalResultRoot\\]/${FinalResultRoot}}"
    s="${s//\\[GinkgoScript\\]/${GinkgoScript}}"
    s="${s//\\[ReadLength\\]/${ReadLength}}"
    s="${s//\\[ResultBaseDir\\]/${ResultBaseDir}}"
    s="${s//\\[SeqID\\]/${SeqID}}"
    s="${s//\\[threads\\]/${threads}}"
  done
  # 앞뒤 공백 제거 및 연속 공백 축소
  echo "$s" | xargs
}

# --- Finalize Outputs & Command ---
FinalResultRoot=$(render "${FinalResultRoot}")
CMD_LINE='for BinSize in [BinSizeList]; do
  echo ">> Running Ginkgo for [SeqID] with BinSize: [BinSize]kb" &&
  [GinkgoScript] --SeqID [SeqID] --ReadLength [ReadLength] --BinSize [BinSize] &&
  
  FinalDir=[ResultBaseDir]/[SeqID]/Binsize_[BinSize]kb &&
  mkdir -p [FinalDir] &&
  
  if [ -d "[ResultBaseDir]/[SeqID]_[BinSize]kb" ]; then
    mv [ResultBaseDir]/[SeqID]_[BinSize]kb/* [FinalDir]/ &&
    rm -rf [ResultBaseDir]/[SeqID]_[BinSize]kb
  fi;
done'
CMD=$(render "$CMD_LINE")

echo -e "\n[RUNNING CMD]\n$CMD\n"

# --- Execution ---
mkdir -p "$CWD"
cd "$CWD" && eval "$CMD"
