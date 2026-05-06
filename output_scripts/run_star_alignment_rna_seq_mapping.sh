#!/bin/bash
# [METADATA]
# TOOL_NAME = star_alignment
# VERSION = 2.7.11a
# THREADS = 1

# Tool Info: star_alignment (2.7.11a)
# Profile: rna_seq_mapping

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           No description"
    echo "  --FastqDir        No description"
    echo "  --BamDir          No description"
    echo "  --StarIndex       Path to STAR index directory"
    echo ""
    echo "Optional Parameters:"
    echo "  --Threads         No description (Default: 12)"
    echo "  --InputSuffix     No description (Default: .trimmed)"
    echo "  --outSAMtype      No description (Default: BAM SortedByCoordinate)"
    echo "  --quantMode       No description (Default: TranscriptomeSAM)"
    echo "  --outSAMattributes No description (Default: NH HI AS nM XS NM)"
    echo "  --chimSegmentMin  No description (Default: 10)"
    echo "  --twopassMode     No description (Default: Basic)"
    echo "  --outFilterMismatchNmax No description (Default: 10)"
    echo "  --outSAMunmapped  No description (Default: Within)"
    echo "  --singularity_bin No description (Default: singularity)"
    echo "  --star_sif        No description (Default: /storage/images/star-2.7.11.sif)"
    echo "  --bind            No description (Default: /storage,/data)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
FastqDir=""
BamDir=""
StarIndex=""
Threads="12"
InputSuffix=".trimmed"
outSAMtype="BAM SortedByCoordinate"
quantMode="TranscriptomeSAM"
outSAMattributes="NH HI AS nM XS NM"
chimSegmentMin="10"
twopassMode="Basic"
outFilterMismatchNmax="10"
outSAMunmapped="Within"
singularity_bin="singularity"
star_sif="/storage/images/star-2.7.11.sif"
bind="/storage,/data"

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --FastqDir) FastqDir="$2"; shift 2 ;;
        --BamDir) BamDir="$2"; shift 2 ;;
        --StarIndex) StarIndex="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        --InputSuffix) InputSuffix="$2"; shift 2 ;;
        --outSAMtype) outSAMtype="$2"; shift 2 ;;
        --quantMode) quantMode="$2"; shift 2 ;;
        --outSAMattributes) outSAMattributes="$2"; shift 2 ;;
        --chimSegmentMin) chimSegmentMin="$2"; shift 2 ;;
        --twopassMode) twopassMode="$2"; shift 2 ;;
        --outFilterMismatchNmax) outFilterMismatchNmax="$2"; shift 2 ;;
        --outSAMunmapped) outSAMunmapped="$2"; shift 2 ;;
        --singularity_bin) singularity_bin="$2"; shift 2 ;;
        --star_sif) star_sif="$2"; shift 2 ;;
        --bind) bind="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${FastqDir:-}" ]]; then echo "Error: --FastqDir is required"; usage; fi
if [[ -z "${BamDir:-}" ]]; then echo "Error: --BamDir is required"; usage; fi
if [[ -z "${StarIndex:-}" ]]; then echo "Error: --StarIndex is required"; usage; fi

# --- [Output Paths] ---
aligned_bam="${BamDir}/${SeqID}.Aligned.sortedByCoord.out.bam"
aligned_transcriptome_bam="${BamDir}/${SeqID}.Aligned.toTranscriptome.out.bam"
gene_counts="${BamDir}/${SeqID}.ReadsPerGene.out.tab"
mapping_log="${BamDir}/${SeqID}.Log.final.out"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${singularity_bin} exec -B ${bind} ${star_sif} STAR --runThreadN ${Threads} --genomeDir ${StarIndex} --readFilesIn ${FastqDir}/${SeqID}${InputSuffix}_R1.fastq.gz ${FastqDir}/${SeqID}${InputSuffix}_R2.fastq.gz --readFilesCommand zcat --quantMode ${quantMode} --twopassMode ${twopassMode} --chimSegmentMin ${chimSegmentMin} --outFilterMismatchNmax ${outFilterMismatchNmax} --outFileNamePrefix ${BamDir}/${SeqID}. --outSAMtype ${outSAMtype} --outSAMattributes ${outSAMattributes} --outSAMunmapped ${outSAMunmapped}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
mkdir -p "$(dirname "${FastqDir}")" 2>/dev/null || mkdir -p "${FastqDir}"
mkdir -p "$(dirname "${BamDir}")" 2>/dev/null || mkdir -p "${BamDir}"

eval "$cmd"