#!/bin/bash
# [METADATA]
# TOOL_NAME = bwa_picard
# VERSION = bwa_0.7.17-picard_3.1.0
# THREADS = 1

# Tool Info: bwa_picard (bwa_0.7.17-picard_3.1.0)
# Profile: pe_align_merge

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           Sequence identifier for sample tracking"
    echo "  --TrimFastqDir    Directory containing input FASTQ files"
    echo "  --BamDir          Directory to store intermediate and final BAM files"
    echo "  --TmpDir          Temporary directory for Picard operations"
    echo "  --ReferenceFasta  Path to the reference genome FASTA file"
    echo "  --ReadGroupID     Read Group ID (ID tag)"
    echo "  --ReadGroupPlatform Sequencing platform (PL tag)"
    echo "  --ReadGroupLibrary Library identifier (LB tag)"
    echo "  --ReadGroupCenter Sequencing center (CN tag)"
    echo ""
    echo "Optional Parameters:"
    echo "  --InputSuffix     Suffix of input FASTQs (e.g., trimmed, raw) (Default: trimmed)"
    echo "  --OutputSuffix    Suffix for the final merged result (Default: primary)"
    echo "  --java_bin        Path to Java executable (Default: java)"
    echo "  --picard_jar      Path to Picard.jar (Default: /storage/apps/bin/picard.jar)"
    echo "  --singularity_bin Path to Singularity (Default: singularity)"
    echo "  --bwa_sif         BWA Singularity image (Default: /storage/images/bwa-0.7.17.sif)"
    echo "  --bind            Mount paths for Singularity (Default: /storage,/data)"
    echo "  --xmx_mb          Max Java heap memory (MB) (Default: 16384)"
    echo "  --Threads         Threads for BWA and Java ParallelGC (Default: 8)"
    echo "  --bwa_bin         No description (Default: bwa)"
    echo "  --mark_short_split No description (Default: -M)"
    echo "  --soft_clipping   No description (Default: -Y)"
    echo "  --clipping_penalty No description (Default: -L 50,50)"
    echo "  --other_args      No description (Default: )"
    echo "  --mba_strategy    No description (Default: MostDistant)"
    echo "  --mba_attributes  No description (Default: XS)"
    echo "  --mba_orientations No description (Default: FR --EXPECTED_ORIENTATIONS RF)"
    echo "  --mba_other_flags No description (Default: --CREATE_INDEX true --MAX_INSERTIONS_OR_DELETIONS -1 --CLIP_ADAPTERS false)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
TrimFastqDir=""
BamDir=""
TmpDir=""
ReferenceFasta=""
ReadGroupID=""
ReadGroupPlatform=""
ReadGroupLibrary=""
ReadGroupCenter=""
InputSuffix="trimmed"
OutputSuffix="primary"
java_bin="java"
picard_jar="/storage/apps/bin/picard.jar"
singularity_bin="singularity"
bwa_sif="/storage/images/bwa-0.7.17.sif"
bind="/storage,/data"
xmx_mb="16384"
Threads="8"
bwa_bin="bwa"
mark_short_split="-M"
soft_clipping="-Y"
clipping_penalty="-L 50,50"
other_args=""
mba_strategy="MostDistant"
mba_attributes="XS"
mba_orientations="FR --EXPECTED_ORIENTATIONS RF"
mba_other_flags="--CREATE_INDEX true --MAX_INSERTIONS_OR_DELETIONS -1 --CLIP_ADAPTERS false"

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --TrimFastqDir) TrimFastqDir="$2"; shift 2 ;;
        --BamDir) BamDir="$2"; shift 2 ;;
        --TmpDir) TmpDir="$2"; shift 2 ;;
        --ReferenceFasta) ReferenceFasta="$2"; shift 2 ;;
        --ReadGroupID) ReadGroupID="$2"; shift 2 ;;
        --ReadGroupPlatform) ReadGroupPlatform="$2"; shift 2 ;;
        --ReadGroupLibrary) ReadGroupLibrary="$2"; shift 2 ;;
        --ReadGroupCenter) ReadGroupCenter="$2"; shift 2 ;;
        --InputSuffix) InputSuffix="$2"; shift 2 ;;
        --OutputSuffix) OutputSuffix="$2"; shift 2 ;;
        --java_bin) java_bin="$2"; shift 2 ;;
        --picard_jar) picard_jar="$2"; shift 2 ;;
        --singularity_bin) singularity_bin="$2"; shift 2 ;;
        --bwa_sif) bwa_sif="$2"; shift 2 ;;
        --bind) bind="$2"; shift 2 ;;
        --xmx_mb) xmx_mb="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        --bwa_bin) bwa_bin="$2"; shift 2 ;;
        --mark_short_split) mark_short_split="$2"; shift 2 ;;
        --soft_clipping) soft_clipping="$2"; shift 2 ;;
        --clipping_penalty) clipping_penalty="$2"; shift 2 ;;
        --other_args) other_args="$2"; shift 2 ;;
        --mba_strategy) mba_strategy="$2"; shift 2 ;;
        --mba_attributes) mba_attributes="$2"; shift 2 ;;
        --mba_orientations) mba_orientations="$2"; shift 2 ;;
        --mba_other_flags) mba_other_flags="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${TrimFastqDir:-}" ]]; then echo "Error: --TrimFastqDir is required"; usage; fi
if [[ -z "${BamDir:-}" ]]; then echo "Error: --BamDir is required"; usage; fi
if [[ -z "${TmpDir:-}" ]]; then echo "Error: --TmpDir is required"; usage; fi
if [[ -z "${ReferenceFasta:-}" ]]; then echo "Error: --ReferenceFasta is required"; usage; fi
if [[ -z "${ReadGroupID:-}" ]]; then echo "Error: --ReadGroupID is required"; usage; fi
if [[ -z "${ReadGroupPlatform:-}" ]]; then echo "Error: --ReadGroupPlatform is required"; usage; fi
if [[ -z "${ReadGroupLibrary:-}" ]]; then echo "Error: --ReadGroupLibrary is required"; usage; fi
if [[ -z "${ReadGroupCenter:-}" ]]; then echo "Error: --ReadGroupCenter is required"; usage; fi

# --- [Output Paths] ---
unmapped_bam="${BamDir}/${SeqID}.${InputSuffix}.u.bam"
aligned_sam="${BamDir}/${SeqID}.${InputSuffix}.bwa_raw.sam"
primary_bam="${BamDir}/${SeqID}.${OutputSuffix}.bam"
primary_bai="${BamDir}/${SeqID}.${OutputSuffix}.bai"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${java_bin} -XX:ParallelGCThreads=${Threads} -Xmx${xmx_mb}m -jar ${picard_jar} FastqToSam --FASTQ ${TrimFastqDir}/${SeqID}.${InputSuffix}_R1.fastq.gz --FASTQ2 ${TrimFastqDir}/${SeqID}.${InputSuffix}_R2.fastq.gz --SAMPLE_NAME ${SeqID} --OUTPUT ${unmapped_bam} --READ_GROUP_NAME ${ReadGroupID} --PLATFORM ${ReadGroupPlatform} --LIBRARY_NAME ${ReadGroupLibrary} --SEQUENCING_CENTER ${ReadGroupCenter} --TMP_DIR ${TmpDir} && ${singularity_bin} exec -B ${bind} ${bwa_sif} ${bwa_bin} mem  ${mark_short_split} ${soft_clipping} ${clipping_penalty} ${other_args}  -t ${Threads} -R '@RG\tID:${ReadGroupID}\tPL:${ReadGroupPlatform}\tLB:${ReadGroupLibrary}\tSM:${SeqID}\tCN:${ReadGroupCenter}' ${ReferenceFasta} ${TrimFastqDir}/${SeqID}.${InputSuffix}_R1.fastq.gz ${TrimFastqDir}/${SeqID}.${InputSuffix}_R2.fastq.gz > ${aligned_sam} && ${java_bin} -XX:ParallelGCThreads=${Threads} -Xmx${xmx_mb}m -jar ${picard_jar} MergeBamAlignment --UNMAPPED_BAM ${unmapped_bam} --ALIGNED_BAM ${aligned_sam} --REFERENCE_SEQUENCE ${ReferenceFasta} --OUTPUT ${primary_bam} --PRIMARY_ALIGNMENT_STRATEGY ${mba_strategy} --ATTRIBUTES_TO_RETAIN ${mba_attributes} --EXPECTED_ORIENTATIONS ${mba_orientations} ${mba_other_flags}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
mkdir -p "$(dirname "${TrimFastqDir}")" 2>/dev/null || mkdir -p "${TrimFastqDir}"
mkdir -p "$(dirname "${BamDir}")" 2>/dev/null || mkdir -p "${BamDir}"
mkdir -p "$(dirname "${TmpDir}")" 2>/dev/null || mkdir -p "${TmpDir}"

eval "$cmd"