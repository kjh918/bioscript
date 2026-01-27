#!/usr/bin/env bash
set -euo pipefail

# Runner: GinkgoAnalysis 1.0.0 (cbNIPT_CopyNumber_Ginkgo)

# defaults (override by CLI)
SeqID=''
ReadLength=''
ResultBaseDir=''
threads=1
BinSizeList='100 150 200 250'
GinkgoScript=/storage/home/kangsm/runScripts/cbNIPT/GCX.cbNIPT_Run.Ginkgo.Analysis.sh
FinalResultRoot='[ResultBaseDir]/[SeqID]'

CWD="."
PRINT_CMD=0
PRINT_OUTPUTS=0
EMIT_OUTPUTS=""

OUTPUT_KEYS=(FinalResultRoot)

usage() {
  echo "Usage: $0 --<Key> <val> ... [--cwd <dir>] [--print-cmd] [--print-outputs] [--emit-outputs <file>]" >&2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --SeqID) SeqID="$2"; shift 2;;
    --ReadLength) ReadLength="$2"; shift 2;;
    --ResultBaseDir) ResultBaseDir="$2"; shift 2;;
    --threads) threads="$2"; shift 2;;
    --BinSizeList) BinSizeList="$2"; shift 2;;
    --GinkgoScript) GinkgoScript="$2"; shift 2;;
    --FinalResultRoot) FinalResultRoot="$2"; shift 2;;
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
[[ -n "${ReadLength}" ]] || { echo "Missing required --ReadLength" >&2; exit 2; }

render() {
  local s="$1"
  s="${s//[SeqID]/${SeqID}}"
  s="${s//[ReadLength]/${ReadLength}}"
  s="${s//[ResultBaseDir]/${ResultBaseDir}}"
  s="${s//[threads]/${threads}}"
  s="${s//[BinSizeList]/${BinSizeList}}"
  s="${s//[GinkgoScript]/${GinkgoScript}}"
  s="${s//[FinalResultRoot]/${FinalResultRoot}}"
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
if [[ -z "${FinalResultRoot}" ]]; then
  __tmp='[ResultBaseDir]/[SeqID]'
  FinalResultRoot="$(render "$__tmp")"
fi

CMD_LINE='for BinSize in [BinSizeList]; do
  echo ">> Running Ginkgo for [SeqID] with BinSize: ${BinSize}kb" &&
  [GinkgoScript] --SeqID [SeqID] --ReadLength [ReadLength] --BinSize ${BinSize} &&
  
  FinalDir=[ResultBaseDir]/[SeqID]/Binsize_${BinSize}kb &&
  mkdir -p ${FinalDir} &&
  
  if [ -d "[ResultBaseDir]/[SeqID]_${BinSize}kb" ]; then
    mv [ResultBaseDir]/[SeqID]_${BinSize}kb/* ${FinalDir}/ &&
    rm -rf [ResultBaseDir]/[SeqID]_${BinSize}kb
  fi;
done'
CMD="$(render "$CMD_LINE")"

if [[ "$PRINT_CMD" -eq 1 ]]; then
  echo "$CMD"
  exit 0
fi

mkdir -p "$CWD"
bash -lc "cd \"$CWD\" && $CMD"

emit_outputs
