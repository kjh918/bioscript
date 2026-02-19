#!/bin/bash
# [METADATA]
# TOOL_NAME = bwa_mem
# VERSION = 0.7.17
# THREADS = 1

# Tool Info: bwa_mem (0.7.17)
# Profile: align_pe_using_singularity

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           Sequence identifier used for file naming and SM tag"
    echo "  --TrimFastqDir    Directory containing FASTQ files"
    echo "  --BamDir          Output directory for aligned BAM files"
    echo "  --ReferenceFasta  Path to the reference genome FASTA file"
    echo "  --ReadGroupID     Read Group ID (ID tag)"
    echo "  --ReadGroupPlatform Sequencing platform (PL tag: e.g., ILLUMINA)"
    echo "  --ReadGroupLibrary Library identifier (LB tag)"
    echo "  --ReadGroupCenter Sequencing center name (CN tag)"
    echo ""
    echo "Optional Parameters:"
    echo "  --InputSuffix     Suffix of input FASTQs (e.g., trimmed, raw, filtered) (Default: trimmed)"
    echo "  --OutputSuffix    Suffix for output BAM (e.g., primary, initial) (Default: primary)"
    echo "  --singularity_bin Path to singularity executable (Default: singularity)"
    echo "  --bwa_bin         BWA binary path inside container (Default: bwa)"
    echo "  --samtools_bin    SAMtools binary path (host or inside if mounted) (Default: samtools)"
    echo "  --sif             BWA singularity image (.sif) (Default: /storage/images/bwa-0.7.17.sif)"
    echo "  --bind            Mount paths for singularity (Default: /storage,/data)"
    echo "  --Threads         Number of threads for BWA and SAMtools (Default: 8)"
    echo "  --mark_short_split Mark shorter split hits as secondary (for Picard compatibility) (Default: -M)"
    echo "  --soft_clipping   Use soft clipping for supplementary alignments (Default: -Y)"
    echo "  --clipping_penalty Penalty for 5'- and 3'-end clipping (Default: -L 50,50)"
    echo "  --other_args      Other additional BWA MEM arguments (Default: )"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
TrimFastqDir=""
BamDir=""
ReferenceFasta=""
ReadGroupID=""
ReadGroupPlatform=""
ReadGroupLibrary=""
ReadGroupCenter=""
InputSuffix="trimmed"
OutputSuffix="primary"
singularity_bin="singularity"
bwa_bin="bwa"
samtools_bin="samtools"
sif="/storage/images/bwa-0.7.17.sif"
bind="/storage,/data"
Threads="8"
mark_short_split="-M"
soft_clipping="-Y"
clipping_penalty="-L 50,50"
other_args=""

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --TrimFastqDir) TrimFastqDir="$2"; shift 2 ;;
        --BamDir) BamDir="$2"; shift 2 ;;
        --ReferenceFasta) ReferenceFasta="$2"; shift 2 ;;
        --ReadGroupID) ReadGroupID="$2"; shift 2 ;;
        --ReadGroupPlatform) ReadGroupPlatform="$2"; shift 2 ;;
        --ReadGroupLibrary) ReadGroupLibrary="$2"; shift 2 ;;
        --ReadGroupCenter) ReadGroupCenter="$2"; shift 2 ;;
        --InputSuffix) InputSuffix="$2"; shift 2 ;;
        --OutputSuffix) OutputSuffix="$2"; shift 2 ;;
        --singularity_bin) singularity_bin="$2"; shift 2 ;;
        --bwa_bin) bwa_bin="$2"; shift 2 ;;
        --samtools_bin) samtools_bin="$2"; shift 2 ;;
        --sif) sif="$2"; shift 2 ;;
        --bind) bind="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        --mark_short_split) mark_short_split="$2"; shift 2 ;;
        --soft_clipping) soft_clipping="$2"; shift 2 ;;
        --clipping_penalty) clipping_penalty="$2"; shift 2 ;;
        --other_args) other_args="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${TrimFastqDir:-}" ]]; then echo "Error: --TrimFastqDir is required"; usage; fi
if [[ -z "${BamDir:-}" ]]; then echo "Error: --BamDir is required"; usage; fi
if [[ -z "${ReferenceFasta:-}" ]]; then echo "Error: --ReferenceFasta is required"; usage; fi
if [[ -z "${ReadGroupID:-}" ]]; then echo "Error: --ReadGroupID is required"; usage; fi
if [[ -z "${ReadGroupPlatform:-}" ]]; then echo "Error: --ReadGroupPlatform is required"; usage; fi
if [[ -z "${ReadGroupLibrary:-}" ]]; then echo "Error: --ReadGroupLibrary is required"; usage; fi
if [[ -z "${ReadGroupCenter:-}" ]]; then echo "Error: --ReadGroupCenter is required"; usage; fi

# --- [Output Paths] ---
primary_bam="${BamDir}/${SeqID}.${OutputSuffix}.bam"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${singularity_bin} exec -B ${bind} ${sif} ${bwa_bin} mem  ${mark_short_split} ${soft_clipping} ${clipping_penalty} ${other_args}  -t ${Threads} -R '@RG\tID:${ReadGroupID}\tPL:${ReadGroupPlatform}\tLB:${ReadGroupLibrary}\tSM:${SeqID}\tCN:${ReadGroupCenter}' ${ReferenceFasta} ${TrimFastqDir}/${SeqID}.${InputSuffix}_R1.fastq.gz ${TrimFastqDir}/${SeqID}.${InputSuffix}_R2.fastq.gz | ${samtools_bin} view -@ ${Threads} -bS -o ${primary_bam} -"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
mkdir -p "$(dirname "${TrimFastqDir}")" 2>/dev/null || mkdir -p "${TrimFastqDir}"
mkdir -p "$(dirname "${BamDir}")" 2>/dev/null || mkdir -p "${BamDir}"

eval "$cmd"