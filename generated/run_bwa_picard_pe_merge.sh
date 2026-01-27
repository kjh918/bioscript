#!/usr/bin/env bash
set -euo pipefail

# Runner: bwa_picard bwa_0.7.17-picard_3.1.0 (pe_merge)

# defaults (override by CLI)
SeqID=''
TrimFastqDir=''
BamDir=''
TmpDir=''
ReferenceFasta=''
ReadGroupID=''
ReadGroupPlatform=''
ReadGroupLibrary=''
ReadGroupCenter=''
java_bin=java
picard_jar=/storage/apps/bin/picard.jar
gc_threads=14
xmx_mb=16384
bind=/storage,/data
bwa_sif=/storage/images/bwa-0.7.17.sif
Threads=8
bwa_args='-M -Y -L 50,50'
mba_args='--CREATE_INDEX true --MAX_INSERTIONS_OR_DELETIONS -1 --CLIP_ADAPTERS false --PRIMARY_ALIGNMENT_STRATEGY MostDistant --ATTRIBUTES_TO_RETAIN XS --EXPECTED_ORIENTATIONS FR --EXPECTED_ORIENTATIONS RF
'
unmapped_bam='[BamDir]/[SeqID].fastqtosam.bam'
aligned_sam='[BamDir]/[SeqID].bwa.mem.sam'
primary_bam='[BamDir]/[SeqID].primary.bam'
primary_bai='[BamDir]/[SeqID].primary.bai'

CWD="."
PRINT_CMD=0
PRINT_OUTPUTS=0
EMIT_OUTPUTS=""

OUTPUT_KEYS=(unmapped_bam aligned_sam primary_bam primary_bai)

usage() {
  echo "Usage: $0 --<Key> <val> ... [--cwd <dir>] [--print-cmd] [--print-outputs] [--emit-outputs <file>]" >&2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --SeqID) SeqID="$2"; shift 2;;
    --TrimFastqDir) TrimFastqDir="$2"; shift 2;;
    --BamDir) BamDir="$2"; shift 2;;
    --TmpDir) TmpDir="$2"; shift 2;;
    --ReferenceFasta) ReferenceFasta="$2"; shift 2;;
    --ReadGroupID) ReadGroupID="$2"; shift 2;;
    --ReadGroupPlatform) ReadGroupPlatform="$2"; shift 2;;
    --ReadGroupLibrary) ReadGroupLibrary="$2"; shift 2;;
    --ReadGroupCenter) ReadGroupCenter="$2"; shift 2;;
    --java_bin) java_bin="$2"; shift 2;;
    --picard_jar) picard_jar="$2"; shift 2;;
    --gc_threads) gc_threads="$2"; shift 2;;
    --xmx_mb) xmx_mb="$2"; shift 2;;
    --bind) bind="$2"; shift 2;;
    --bwa_sif) bwa_sif="$2"; shift 2;;
    --Threads) Threads="$2"; shift 2;;
    --bwa_args) bwa_args="$2"; shift 2;;
    --mba_args) mba_args="$2"; shift 2;;
    --unmapped_bam) unmapped_bam="$2"; shift 2;;
    --aligned_sam) aligned_sam="$2"; shift 2;;
    --primary_bam) primary_bam="$2"; shift 2;;
    --primary_bai) primary_bai="$2"; shift 2;;
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
[[ -n "${TrimFastqDir}" ]] || { echo "Missing required --TrimFastqDir" >&2; exit 2; }
[[ -n "${BamDir}" ]] || { echo "Missing required --BamDir" >&2; exit 2; }
[[ -n "${TmpDir}" ]] || { echo "Missing required --TmpDir" >&2; exit 2; }
[[ -n "${ReferenceFasta}" ]] || { echo "Missing required --ReferenceFasta" >&2; exit 2; }
[[ -n "${ReadGroupID}" ]] || { echo "Missing required --ReadGroupID" >&2; exit 2; }
[[ -n "${ReadGroupPlatform}" ]] || { echo "Missing required --ReadGroupPlatform" >&2; exit 2; }
[[ -n "${ReadGroupLibrary}" ]] || { echo "Missing required --ReadGroupLibrary" >&2; exit 2; }
[[ -n "${ReadGroupCenter}" ]] || { echo "Missing required --ReadGroupCenter" >&2; exit 2; }

render() {
  local s="$1"
  s="${s//[SeqID]/${SeqID}}"
  s="${s//[TrimFastqDir]/${TrimFastqDir}}"
  s="${s//[BamDir]/${BamDir}}"
  s="${s//[TmpDir]/${TmpDir}}"
  s="${s//[ReferenceFasta]/${ReferenceFasta}}"
  s="${s//[ReadGroupID]/${ReadGroupID}}"
  s="${s//[ReadGroupPlatform]/${ReadGroupPlatform}}"
  s="${s//[ReadGroupLibrary]/${ReadGroupLibrary}}"
  s="${s//[ReadGroupCenter]/${ReadGroupCenter}}"
  s="${s//[java_bin]/${java_bin}}"
  s="${s//[picard_jar]/${picard_jar}}"
  s="${s//[gc_threads]/${gc_threads}}"
  s="${s//[xmx_mb]/${xmx_mb}}"
  s="${s//[bind]/${bind}}"
  s="${s//[bwa_sif]/${bwa_sif}}"
  s="${s//[Threads]/${Threads}}"
  s="${s//[bwa_args]/${bwa_args}}"
  s="${s//[mba_args]/${mba_args}}"
  s="${s//[unmapped_bam]/${unmapped_bam}}"
  s="${s//[aligned_sam]/${aligned_sam}}"
  s="${s//[primary_bam]/${primary_bam}}"
  s="${s//[primary_bai]/${primary_bai}}"
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
if [[ -z "${unmapped_bam}" ]]; then
  __tmp='[BamDir]/[SeqID].fastqtosam.bam'
  unmapped_bam="$(render "$__tmp")"
fi
if [[ -z "${aligned_sam}" ]]; then
  __tmp='[BamDir]/[SeqID].bwa.mem.sam'
  aligned_sam="$(render "$__tmp")"
fi
if [[ -z "${primary_bam}" ]]; then
  __tmp='[BamDir]/[SeqID].primary.bam'
  primary_bam="$(render "$__tmp")"
fi
if [[ -z "${primary_bai}" ]]; then
  __tmp='[BamDir]/[SeqID].primary.bai'
  primary_bai="$(render "$__tmp")"
fi

CMD_LINE='[java_bin] -XX:ParallelGCThreads=[gc_threads] -Xmx[xmx_mb]m -jar [picard_jar] FastqToSam --FASTQ [TrimFastqDir]/[SeqID].trimmed_R1.fastq.gz --FASTQ2 [TrimFastqDir]/[SeqID].trimmed_R2.fastq.gz --SAMPLE_NAME [SeqID] --OUTPUT [BamDir]/[SeqID].fastqtosam.bam --READ_GROUP_NAME [ReadGroupID] --PLATFORM [ReadGroupPlatform] --LIBRARY_NAME [ReadGroupLibrary] --SEQUENCING_CENTER [ReadGroupCenter] --TMP_DIR [TmpDir] && singularity exec -B [bind] [bwa_sif] bwa mem [bwa_args] -t [Threads] -R "@RG\tID:[ReadGroupID]\tPL:[ReadGroupPlatform]\tLB:[ReadGroupLibrary]\tSM:[SeqID]\tCN:[ReadGroupCenter]" [ReferenceFasta] [TrimFastqDir]/[SeqID].trimmed_R1.fastq.gz [TrimFastqDir]/[SeqID].trimmed_R2.fastq.gz > [BamDir]/[SeqID].bwa.mem.sam && [java_bin] -XX:ParallelGCThreads=[gc_threads] -Xmx[xmx_mb]m -jar [picard_jar] MergeBamAlignment --UNMAPPED_BAM [BamDir]/[SeqID].fastqtosam.bam --ALIGNED_BAM [BamDir]/[SeqID].bwa.mem.sam --REFERENCE_SEQUENCE [ReferenceFasta] --OUTPUT [BamDir]/[SeqID].primary.bam [mba_args]'
CMD="$(render "$CMD_LINE")"

if [[ "$PRINT_CMD" -eq 1 ]]; then
  echo "$CMD"
  exit 0
fi

mkdir -p "$CWD"
bash -lc "cd \"$CWD\" && $CMD"

emit_outputs
