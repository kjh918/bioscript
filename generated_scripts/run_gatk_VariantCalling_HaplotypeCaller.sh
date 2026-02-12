#!/usr/bin/env bash
set -euo pipefail

# --- Defaults ---
Chromosome=''
DbsnpVcf=''
GatkPath=/storage/apps/gatk-4.4.0.0/gatk
IncludeNonVariant=false
InputBam=''
Memory=32g
OutGvcf='[ResultDir]/[SeqID].[Chromosome].gvcf.gz'
OutVcf='[ResultDir]/[SeqID].[Chromosome].vcf.gz'
Ploidy=2
ReferenceFasta=''
ResultDir=''
SeqID=''
Threads=4
TmpDir='[ResultDir]/tmp'
CWD="."

# --- CLI Argument Parsing ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    --Chromosome) Chromosome="$2"; shift 2;;
    --DbsnpVcf) DbsnpVcf="$2"; shift 2;;
    --GatkPath) GatkPath="$2"; shift 2;;
    --IncludeNonVariant) IncludeNonVariant="$2"; shift 2;;
    --InputBam) InputBam="$2"; shift 2;;
    --Memory) Memory="$2"; shift 2;;
    --OutGvcf) OutGvcf="$2"; shift 2;;
    --OutVcf) OutVcf="$2"; shift 2;;
    --Ploidy) Ploidy="$2"; shift 2;;
    --ReferenceFasta) ReferenceFasta="$2"; shift 2;;
    --ResultDir) ResultDir="$2"; shift 2;;
    --SeqID) SeqID="$2"; shift 2;;
    --Threads) Threads="$2"; shift 2;;
    --TmpDir) TmpDir="$2"; shift 2;;
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
    s="${s//\\[Chromosome\\]/${Chromosome}}"
    s="${s//\\[DbsnpVcf\\]/${DbsnpVcf}}"
    s="${s//\\[GatkPath\\]/${GatkPath}}"
    s="${s//\\[IncludeNonVariant\\]/${IncludeNonVariant}}"
    s="${s//\\[InputBam\\]/${InputBam}}"
    s="${s//\\[Memory\\]/${Memory}}"
    s="${s//\\[OutGvcf\\]/${OutGvcf}}"
    s="${s//\\[OutVcf\\]/${OutVcf}}"
    s="${s//\\[Ploidy\\]/${Ploidy}}"
    s="${s//\\[ReferenceFasta\\]/${ReferenceFasta}}"
    s="${s//\\[ResultDir\\]/${ResultDir}}"
    s="${s//\\[SeqID\\]/${SeqID}}"
    s="${s//\\[Threads\\]/${Threads}}"
    s="${s//\\[TmpDir\\]/${TmpDir}}"
  done
  # 앞뒤 공백 제거 및 연속 공백 축소
  echo "$s" | xargs
}

# --- Finalize Outputs & Command ---
OutVcf=$(render "${OutVcf}")
OutGvcf=$(render "${OutGvcf}")
CMD_LINE='mkdir -p [TmpDir] && 
[GatkPath] --java-options "-XX:ParallelGCThreads=[Threads] -Xmx[Memory] -Djava.io.tmpdir=[TmpDir]"  HaplotypeCaller  -R [ReferenceFasta] -I [InputBam] -L [Chromosome]  -ploidy [Ploidy] -stand-call-conf 30 --dbsnp [DbsnpVcf]  -O [OutGvcf] -ERC GVCF && 
[GatkPath] --java-options "-XX:ParallelGCThreads=[Threads] -Xmx[Memory] -Djava.io.tmpdir=[TmpDir]"  GenotypeGVCFs  --include-non-variant-sites [IncludeNonVariant]  -R [ReferenceFasta] -V [OutGvcf] -O [OutVcf]'
CMD=$(render "$CMD_LINE")

echo -e "\n[RUNNING CMD]\n$CMD\n"

# --- Execution ---
mkdir -p "$CWD"
cd "$CWD" && eval "$CMD"
