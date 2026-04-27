#!/bin/bash
# [METADATA]
# TOOL_NAME = bcftools
# VERSION = 1.23
# THREADS = 1

# Tool Info: bcftools (1.23)
# Profile: QC_PileupCall_Annotation

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           Sequence identifier for the sample"
    echo "  --Chromosome      Target chromosome or region"
    echo "  --BamDir          Directory containing the input BAM file"
    echo "  --ResultDir       Output directory for VCF results"
    echo "  --ReferenceFasta  Path to reference genome FASTA"
    echo "  --SitesVcfGz      VCF file containing specific sites to genotype"
    echo "  --PopAfAnnotVcf   VCF/Tabix file for population AF annotation"
    echo "  --PopAfHeaderHdr  Header file describing the annotation columns"
    echo ""
    echo "Optional Parameters:"
    echo "  --InputSuffix     Suffix of input BAM (e.g., analysisReady, recal) (Default: analysisReady)"
    echo "  --bcftools_bin    Path to bcftools binary (Default: /storage/home/jhkim/Apps/bcftools/bcftools)"
    echo "  --Threads         Number of threads for compression/processing (Default: 4)"
    echo "  --MinBQ           Minimum base quality for mpileup (Default: 20)"
    echo "  --MinMQ           Minimum mapping quality for mpileup (Default: 30)"
    echo "  --AnnotationQuery Columns to extract from annotation VCF (Default: CHROM,POS,REF,ALT,INFO/KOVA_AF,INFO/KOVA_AN,INFO/GNOMAD_AF,INFO/GNOMAD_AN)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
Chromosome=""
BamDir=""
ResultDir=""
ReferenceFasta=""
SitesVcfGz=""
PopAfAnnotVcf=""
PopAfHeaderHdr=""
InputSuffix="analysisReady"
bcftools_bin="/storage/home/jhkim/Apps/bcftools/bcftools"
Threads="4"
MinBQ="20"
MinMQ="30"
AnnotationQuery="CHROM,POS,REF,ALT,INFO/KOVA_AF,INFO/KOVA_AN,INFO/GNOMAD_AF,INFO/GNOMAD_AN"

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --Chromosome) Chromosome="$2"; shift 2 ;;
        --BamDir) BamDir="$2"; shift 2 ;;
        --ResultDir) ResultDir="$2"; shift 2 ;;
        --ReferenceFasta) ReferenceFasta="$2"; shift 2 ;;
        --SitesVcfGz) SitesVcfGz="$2"; shift 2 ;;
        --PopAfAnnotVcf) PopAfAnnotVcf="$2"; shift 2 ;;
        --PopAfHeaderHdr) PopAfHeaderHdr="$2"; shift 2 ;;
        --InputSuffix) InputSuffix="$2"; shift 2 ;;
        --bcftools_bin) bcftools_bin="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        --MinBQ) MinBQ="$2"; shift 2 ;;
        --MinMQ) MinMQ="$2"; shift 2 ;;
        --AnnotationQuery) AnnotationQuery="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${Chromosome:-}" ]]; then echo "Error: --Chromosome is required"; usage; fi
if [[ -z "${BamDir:-}" ]]; then echo "Error: --BamDir is required"; usage; fi
if [[ -z "${ResultDir:-}" ]]; then echo "Error: --ResultDir is required"; usage; fi
if [[ -z "${ReferenceFasta:-}" ]]; then echo "Error: --ReferenceFasta is required"; usage; fi
if [[ -z "${SitesVcfGz:-}" ]]; then echo "Error: --SitesVcfGz is required"; usage; fi
if [[ -z "${PopAfAnnotVcf:-}" ]]; then echo "Error: --PopAfAnnotVcf is required"; usage; fi
if [[ -z "${PopAfHeaderHdr:-}" ]]; then echo "Error: --PopAfHeaderHdr is required"; usage; fi

# --- [Output Paths] ---
raw_vcf="${ResultDir}/${SeqID}.${Chromosome}.${InputSuffix}.raw.vcf.gz"
ann_vcf="${ResultDir}/${SeqID}.${Chromosome}.${InputSuffix}.ann.vcf.gz"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${bcftools_bin} mpileup -f ${ReferenceFasta} -T ${SitesVcfGz} -r ${Chromosome} -q ${MinMQ} -Q ${MinBQ} -a FORMAT/AD,FORMAT/DP -Ou --threads ${Threads} ${BamDir}/${SeqID}.${InputSuffix}.bam | ${bcftools_bin} call -Am -Oz -o ${raw_vcf} --threads ${Threads} && ${bcftools_bin} index -f ${raw_vcf} --threads ${Threads} && ${bcftools_bin} annotate -a ${PopAfAnnotVcf} -c ${AnnotationQuery} -h ${PopAfHeaderHdr} -Oz -o ${ann_vcf} ${raw_vcf} --threads ${Threads} && ${bcftools_bin} index -f ${ann_vcf} --threads ${Threads}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
mkdir -p "$(dirname "${BamDir}")" 2>/dev/null || mkdir -p "${BamDir}"
mkdir -p "$(dirname "${ResultDir}")" 2>/dev/null || mkdir -p "${ResultDir}"

eval "$cmd"