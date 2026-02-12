#!/usr/bin/env bash
set -euo pipefail

# --- Defaults ---
BamDir=''
ReadGroupCenter=''
ReadGroupID=''
ReadGroupLibrary=''
ReadGroupPlatform=''
ReferenceFasta=''
SeqID=''
Threads=8
TmpDir=''
TrimFastqDir=''
aligned_sam='[BamDir]/[SeqID].bwa.mem.sam'
bind=/storage,/data
bwa_args='-M -Y -L 50,50'
bwa_sif=/storage/images/bwa-0.7.17.sif
java_bin=java
mba_args='--CREATE_INDEX true --MAX_INSERTIONS_OR_DELETIONS -1 --CLIP_ADAPTERS false --PRIMARY_ALIGNMENT_STRATEGY MostDistant --ATTRIBUTES_TO_RETAIN XS --EXPECTED_ORIENTATIONS FR --EXPECTED_ORIENTATIONS RF
'
picard_jar=/storage/apps/bin/picard.jar
primary_bai='[BamDir]/[SeqID].primary.bai'
primary_bam='[BamDir]/[SeqID].primary.bam'
unmapped_bam='[BamDir]/[SeqID].fastqtosam.bam'
xmx_mb=16384
CWD="."

# --- CLI Argument Parsing ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    --BamDir) BamDir="$2"; shift 2;;
    --ReadGroupCenter) ReadGroupCenter="$2"; shift 2;;
    --ReadGroupID) ReadGroupID="$2"; shift 2;;
    --ReadGroupLibrary) ReadGroupLibrary="$2"; shift 2;;
    --ReadGroupPlatform) ReadGroupPlatform="$2"; shift 2;;
    --ReferenceFasta) ReferenceFasta="$2"; shift 2;;
    --SeqID) SeqID="$2"; shift 2;;
    --Threads) Threads="$2"; shift 2;;
    --TmpDir) TmpDir="$2"; shift 2;;
    --TrimFastqDir) TrimFastqDir="$2"; shift 2;;
    --aligned_sam) aligned_sam="$2"; shift 2;;
    --bind) bind="$2"; shift 2;;
    --bwa_args) bwa_args="$2"; shift 2;;
    --bwa_sif) bwa_sif="$2"; shift 2;;
    --java_bin) java_bin="$2"; shift 2;;
    --mba_args) mba_args="$2"; shift 2;;
    --picard_jar) picard_jar="$2"; shift 2;;
    --primary_bai) primary_bai="$2"; shift 2;;
    --primary_bam) primary_bam="$2"; shift 2;;
    --unmapped_bam) unmapped_bam="$2"; shift 2;;
    --xmx_mb) xmx_mb="$2"; shift 2;;
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
    s="${s//\\[BamDir\\]/${BamDir}}"
    s="${s//\\[ReadGroupCenter\\]/${ReadGroupCenter}}"
    s="${s//\\[ReadGroupID\\]/${ReadGroupID}}"
    s="${s//\\[ReadGroupLibrary\\]/${ReadGroupLibrary}}"
    s="${s//\\[ReadGroupPlatform\\]/${ReadGroupPlatform}}"
    s="${s//\\[ReferenceFasta\\]/${ReferenceFasta}}"
    s="${s//\\[SeqID\\]/${SeqID}}"
    s="${s//\\[Threads\\]/${Threads}}"
    s="${s//\\[TmpDir\\]/${TmpDir}}"
    s="${s//\\[TrimFastqDir\\]/${TrimFastqDir}}"
    s="${s//\\[aligned_sam\\]/${aligned_sam}}"
    s="${s//\\[bind\\]/${bind}}"
    s="${s//\\[bwa_args\\]/${bwa_args}}"
    s="${s//\\[bwa_sif\\]/${bwa_sif}}"
    s="${s//\\[java_bin\\]/${java_bin}}"
    s="${s//\\[mba_args\\]/${mba_args}}"
    s="${s//\\[picard_jar\\]/${picard_jar}}"
    s="${s//\\[primary_bai\\]/${primary_bai}}"
    s="${s//\\[primary_bam\\]/${primary_bam}}"
    s="${s//\\[unmapped_bam\\]/${unmapped_bam}}"
    s="${s//\\[xmx_mb\\]/${xmx_mb}}"
  done
  # 앞뒤 공백 제거 및 연속 공백 축소
  echo "$s" | xargs
}

# --- Finalize Outputs & Command ---
unmapped_bam=$(render "${unmapped_bam}")
aligned_sam=$(render "${aligned_sam}")
primary_bam=$(render "${primary_bam}")
primary_bai=$(render "${primary_bai}")
CMD_LINE='[java_bin] -XX:ParallelGCThreads=[Threads] -Xmx[xmx_mb]m -jar [picard_jar] FastqToSam --FASTQ [TrimFastqDir]/[SeqID].trimmed_R1.fastq.gz --FASTQ2 [TrimFastqDir]/[SeqID].trimmed_R2.fastq.gz --SAMPLE_NAME [SeqID] --OUTPUT [BamDir]/[SeqID].fastqtosam.bam --READ_GROUP_NAME [ReadGroupID] --PLATFORM [ReadGroupPlatform] --LIBRARY_NAME [ReadGroupLibrary] --SEQUENCING_CENTER [ReadGroupCenter] --TMP_DIR [TmpDir] && singularity exec -B [bind] [bwa_sif] bwa mem [bwa_args] -t [Threads] -R "@RG\tID:[ReadGroupID]\tPL:[ReadGroupPlatform]\tLB:[ReadGroupLibrary]\tSM:[SeqID]\tCN:[ReadGroupCenter]" [ReferenceFasta] [TrimFastqDir]/[SeqID].trimmed_R1.fastq.gz [TrimFastqDir]/[SeqID].trimmed_R2.fastq.gz > [BamDir]/[SeqID].bwa.mem.sam && [java_bin] -XX:ParallelGCThreads=[Threads] -Xmx[xmx_mb]m -jar [picard_jar] MergeBamAlignment --UNMAPPED_BAM [BamDir]/[SeqID].fastqtosam.bam --ALIGNED_BAM [BamDir]/[SeqID].bwa.mem.sam --REFERENCE_SEQUENCE [ReferenceFasta] --OUTPUT [BamDir]/[SeqID].primary.bam [mba_args]'
CMD=$(render "$CMD_LINE")

echo -e "\n[RUNNING CMD]\n$CMD\n"

# --- Execution ---
mkdir -p "$CWD"
cd "$CWD" && eval "$CMD"
