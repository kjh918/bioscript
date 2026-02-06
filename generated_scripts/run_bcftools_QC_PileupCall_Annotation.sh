#!/usr/bin/env bash
set -euo pipefail

# --- Defaults ---
AnnVcf='[ResultDir]/[SeqID].[Chromosome].sites.af.vcf.gz'
AnnotationQuery=CHROM,POS,REF,ALT,INFO/AFPOP:=5,INFO/KOVA_AN:=6,INFO/GNOMAD_AF:=7,INFO/GNOMAD_AN:=8
BcftoolsPath=/storage/apps/bcftools-1.3.1/bin/bcftools
BgzipPath=bgzip
Chromosome=''
InputBamDir=''
MinBQ=20
MinMQ=30
OutCountsTsvGz='[ResultDir]/[SeqID].[Chromosome].counts.tsv.gz'
PopAfAnnotTsv=''
PopAfHeaderHdr=''
RawVcf='[ResultDir]/[SeqID].[Chromosome].sites.raw.vcf.gz'
ReferenceFasta=''
ResultDir=''
SeqID=''
SitesVcfGz=''
Threads=4
TmpDir='[ResultDir]/tmp'
VcfQuery='%CHROM|%POS|%REF|%ALT|%INFO/AFPOP|%INFO/KOVA_AN|%INFO/GNOMAD_AF|%INFO/GNOMAD_AN|[%GT]|[%DP]|[%AD]\n'
CWD="."

# --- CLI Argument Parsing ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    --AnnVcf) AnnVcf="$2"; shift 2;;
    --AnnotationQuery) AnnotationQuery="$2"; shift 2;;
    --BcftoolsPath) BcftoolsPath="$2"; shift 2;;
    --BgzipPath) BgzipPath="$2"; shift 2;;
    --Chromosome) Chromosome="$2"; shift 2;;
    --InputBamDir) InputBamDir="$2"; shift 2;;
    --MinBQ) MinBQ="$2"; shift 2;;
    --MinMQ) MinMQ="$2"; shift 2;;
    --OutCountsTsvGz) OutCountsTsvGz="$2"; shift 2;;
    --PopAfAnnotTsv) PopAfAnnotTsv="$2"; shift 2;;
    --PopAfHeaderHdr) PopAfHeaderHdr="$2"; shift 2;;
    --RawVcf) RawVcf="$2"; shift 2;;
    --ReferenceFasta) ReferenceFasta="$2"; shift 2;;
    --ResultDir) ResultDir="$2"; shift 2;;
    --SeqID) SeqID="$2"; shift 2;;
    --SitesVcfGz) SitesVcfGz="$2"; shift 2;;
    --Threads) Threads="$2"; shift 2;;
    --TmpDir) TmpDir="$2"; shift 2;;
    --VcfQuery) VcfQuery="$2"; shift 2;;
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
    s="${s//\\[AnnVcf\\]/${AnnVcf}}"
    s="${s//\\[AnnotationQuery\\]/${AnnotationQuery}}"
    s="${s//\\[BcftoolsPath\\]/${BcftoolsPath}}"
    s="${s//\\[BgzipPath\\]/${BgzipPath}}"
    s="${s//\\[Chromosome\\]/${Chromosome}}"
    s="${s//\\[InputBamDir\\]/${InputBamDir}}"
    s="${s//\\[MinBQ\\]/${MinBQ}}"
    s="${s//\\[MinMQ\\]/${MinMQ}}"
    s="${s//\\[OutCountsTsvGz\\]/${OutCountsTsvGz}}"
    s="${s//\\[PopAfAnnotTsv\\]/${PopAfAnnotTsv}}"
    s="${s//\\[PopAfHeaderHdr\\]/${PopAfHeaderHdr}}"
    s="${s//\\[RawVcf\\]/${RawVcf}}"
    s="${s//\\[ReferenceFasta\\]/${ReferenceFasta}}"
    s="${s//\\[ResultDir\\]/${ResultDir}}"
    s="${s//\\[SeqID\\]/${SeqID}}"
    s="${s//\\[SitesVcfGz\\]/${SitesVcfGz}}"
    s="${s//\\[Threads\\]/${Threads}}"
    s="${s//\\[TmpDir\\]/${TmpDir}}"
    s="${s//\\[VcfQuery\\]/${VcfQuery}}"
  done
  # 앞뒤 공백 제거 및 연속 공백 축소
  echo "$s" | xargs
}

# --- Finalize Outputs & Command ---
RawVcf=$(render "${RawVcf}")
AnnVcf=$(render "${AnnVcf}")
OutCountsTsvGz=$(render "${OutCountsTsvGz}")
CMD_LINE='[BcftoolsPath] mpileup -f [ReferenceFasta] -T [SitesVcfGz] -q [MinMQ] -Q [MinBQ] -a FORMAT/AD,FORMAT/DP -Ou [InputBamDir]/[SeqID].analysisReady.bam | [BcftoolsPath] call -Aim -Oz -o [RawVcf] && [BcftoolsPath] index -f [RawVcf] && [BcftoolsPath] annotate -a [PopAfAnnotTsv] -c [AnnotationQuery] -h [PopAfHeaderHdr] -Oz -o [AnnVcf] [RawVcf] && [BcftoolsPath] index -f [AnnVcf]'
CMD=$(render "$CMD_LINE")

echo -e "\n[RUNNING CMD]\n$CMD\n"

# --- Execution ---
mkdir -p "$CWD"
cd "$CWD" && eval "$CMD"
