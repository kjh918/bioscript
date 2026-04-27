#!/bin/bash
# [METADATA]
# TOOL_NAME = vcf2maf
# VERSION = 1.6.21
# THREADS = 1

# Tool Info: vcf2maf (1.6.21)
# Profile: maf_conversion

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           No description"
    echo "  --vcfDir          No description"
    echo "  --VcfTag          Input VCF tag (e.g. mutect2.filtered.vep.refseq)"
    echo "  --NormalID        Normal sample ID"
    echo "  --genomeFasta     No description"
    echo ""
    echo "Optional Parameters:"
    echo "  --retain_info     MAF에 컬럼으로 남길 VCF INFO 필드 리스트 (Default: DP,GERMQ,MBQ,MFRL,MMQ,MPOS,POPAF,PON,ROQ,TLOD,FILTER,pArtifact,artiStatus)"
    echo "  --retain_ann      MAF에 컬럼으로 남길 VEP ANN(CSQ) 필드 리스트 (Default: HGVSg,am_class,am_pathogenicity,CADD_PHRED,CADD_RAW,ClinVar_CLNSIG,ClinVar_CLNDN,ClinVar_ORIGIN,gnomADe_AF,gnomADe_AFR_AF,gnomADe_AMR_AF,gnomADe_ASJ_AF,gnomADe_EAS_AF,gnomADe_FIN_AF,gnomADe_NFE_AF,gnomADe_OTH_AF,gnomADe_SAS_AF,MAX_AF,MAX_AF_POPS,TRANSCRIPTION_FACTORS,COSMIC,COSMIC_HGVSC,COSMIC_LEGACY_ID,COSMIC_TIER)"
    echo "  --singularity_bin No description (Default: singularity)"
    echo "  --vcf2maf_sif     No description (Default: /storage/images/vcf2maf.sif)"
    echo "  --vcf2maf_bin     No description (Default: vcf2maf.pl)"
    echo "  --bind            No description (Default: /storage,/data)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
vcfDir=""
VcfTag=""
NormalID="NORMAL"
genomeFasta=""
retain_info="DP,GERMQ,MBQ,MFRL,MMQ,MPOS,POPAF,PON,ROQ,TLOD,FILTER,pArtifact,artiStatus"
retain_ann="HGVSg,am_class,am_pathogenicity,CADD_PHRED,CADD_RAW,ClinVar_CLNSIG,ClinVar_CLNDN,ClinVar_ORIGIN,gnomADe_AF,gnomADe_AFR_AF,gnomADe_AMR_AF,gnomADe_ASJ_AF,gnomADe_EAS_AF,gnomADe_FIN_AF,gnomADe_NFE_AF,gnomADe_OTH_AF,gnomADe_SAS_AF,MAX_AF,MAX_AF_POPS,TRANSCRIPTION_FACTORS,COSMIC,COSMIC_HGVSC,COSMIC_LEGACY_ID,COSMIC_TIER"
singularity_bin="singularity"
vcf2maf_sif="/storage/images/vcf2maf.sif"
vcf2maf_bin="vcf2maf.pl"
bind="/storage,/data"

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --vcfDir) vcfDir="$2"; shift 2 ;;
        --VcfTag) VcfTag="$2"; shift 2 ;;
        --NormalID) NormalID="$2"; shift 2 ;;
        --genomeFasta) genomeFasta="$2"; shift 2 ;;
        --retain_info) retain_info="$2"; shift 2 ;;
        --retain_ann) retain_ann="$2"; shift 2 ;;
        --singularity_bin) singularity_bin="$2"; shift 2 ;;
        --vcf2maf_sif) vcf2maf_sif="$2"; shift 2 ;;
        --vcf2maf_bin) vcf2maf_bin="$2"; shift 2 ;;
        --bind) bind="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${vcfDir:-}" ]]; then echo "Error: --vcfDir is required"; usage; fi
if [[ -z "${VcfTag:-}" ]]; then echo "Error: --VcfTag is required"; usage; fi
if [[ -z "${NormalID:-}" ]]; then echo "Error: --NormalID is required"; usage; fi
if [[ -z "${genomeFasta:-}" ]]; then echo "Error: --genomeFasta is required"; usage; fi

# --- [Output Paths] ---
output_maf="${vcfDir}/${SeqID}.${VcfTag}.maf"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${singularity_bin} exec -B ${bind} ${vcf2maf_sif} ${vcf2maf_bin} --inhibit-vep --input-vcf ${vcfDir}/${SeqID}.${VcfTag}.vcf --output-maf ${output_maf} --tumor-id ${SeqID} --normal-id ${NormalID} --ref-fasta ${genomeFasta} --retain-info ${retain_info} --retain-ann ${retain_ann}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
mkdir -p "$(dirname "${vcfDir}")" 2>/dev/null || mkdir -p "${vcfDir}"

eval "$cmd"