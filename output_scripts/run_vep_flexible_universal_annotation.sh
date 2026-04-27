#!/bin/bash
# [METADATA]
# TOOL_NAME = vep_flexible
# VERSION = 110
# THREADS = 1

# Tool Info: vep_flexible (110)
# Profile: universal_annotation

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           No description"
    echo "  --vcfDir          No description"
    echo "  --VcfTag          Input VCF tag (e.g. mutect2.bias.filtered)"
    echo "  --vepCacheDir     No description"
    echo "  --vepPluginDir    No description"
    echo "  --genomeFasta     No description"
    echo "  --clinvarData     ClinVar VCF path (Assembly-specific)"
    echo "  --alphaMissense   No description"
    echo "  --caddSnp         No description"
    echo "  --caddIndel       No description"
    echo ""
    echo "Optional Parameters:"
    echo "  --cosmicData      COSMIC VCF path (Assembly-specific, optional) (Default: )"
    echo "  --assembly        GRCh38 or GRCh37 (Default: GRCh38)"
    echo "  --species         No description (Default: homo_sapiens)"
    echo "  --Threads         No description (Default: 8)"
    echo "  --buffer_size     No description (Default: 50000)"
    echo "  --cosmic_args     Custom flag for COSMIC (e.g. --custom file=...) (Default: )"
    echo "  --clinvar_args    No description (Default: --custom file=[clinvarData],short_name=ClinVar,format=vcf,fields=CLNSIG%CLNDN%ORIGIN)"
    echo "  --singularity_bin No description (Default: singularity)"
    echo "  --vep_sif         No description (Default: /storage/images/vep-110.sif)"
    echo "  --bind            No description (Default: /storage,/data)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
vcfDir=""
VcfTag=""
vepCacheDir=""
vepPluginDir=""
genomeFasta=""
clinvarData=""
cosmicData=""
alphaMissense=""
caddSnp=""
caddIndel=""
assembly="GRCh38"
species="homo_sapiens"
Threads="8"
buffer_size="50000"
cosmic_args=""
clinvar_args="--custom file=[clinvarData],short_name=ClinVar,format=vcf,fields=CLNSIG%CLNDN%ORIGIN"
singularity_bin="singularity"
vep_sif="/storage/images/vep-110.sif"
bind="/storage,/data"

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --vcfDir) vcfDir="$2"; shift 2 ;;
        --VcfTag) VcfTag="$2"; shift 2 ;;
        --vepCacheDir) vepCacheDir="$2"; shift 2 ;;
        --vepPluginDir) vepPluginDir="$2"; shift 2 ;;
        --genomeFasta) genomeFasta="$2"; shift 2 ;;
        --clinvarData) clinvarData="$2"; shift 2 ;;
        --cosmicData) cosmicData="$2"; shift 2 ;;
        --alphaMissense) alphaMissense="$2"; shift 2 ;;
        --caddSnp) caddSnp="$2"; shift 2 ;;
        --caddIndel) caddIndel="$2"; shift 2 ;;
        --assembly) assembly="$2"; shift 2 ;;
        --species) species="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        --buffer_size) buffer_size="$2"; shift 2 ;;
        --cosmic_args) cosmic_args="$2"; shift 2 ;;
        --clinvar_args) clinvar_args="$2"; shift 2 ;;
        --singularity_bin) singularity_bin="$2"; shift 2 ;;
        --vep_sif) vep_sif="$2"; shift 2 ;;
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
if [[ -z "${vepCacheDir:-}" ]]; then echo "Error: --vepCacheDir is required"; usage; fi
if [[ -z "${vepPluginDir:-}" ]]; then echo "Error: --vepPluginDir is required"; usage; fi
if [[ -z "${genomeFasta:-}" ]]; then echo "Error: --genomeFasta is required"; usage; fi
if [[ -z "${clinvarData:-}" ]]; then echo "Error: --clinvarData is required"; usage; fi
if [[ -z "${alphaMissense:-}" ]]; then echo "Error: --alphaMissense is required"; usage; fi
if [[ -z "${caddSnp:-}" ]]; then echo "Error: --caddSnp is required"; usage; fi
if [[ -z "${caddIndel:-}" ]]; then echo "Error: --caddIndel is required"; usage; fi

# --- [Output Paths] ---
vep_vcf="${vcfDir}/${SeqID}.${VcfTag}.vep.refseq.vcf"
vep_stats="${vcfDir}/${SeqID}.${VcfTag}.vep.refseq.vcf_summary.txt"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${singularity_bin} exec -B ${bind} ${vep_sif} vep --force_overwrite --offline --cache --dir_cache ${vepCacheDir} --dir_plugins ${vepPluginDir} --fasta ${genomeFasta} --species ${species} --assembly ${assembly} --input_file ${vcfDir}/${SeqID}.${VcfTag}.vcf --output_file ${vep_vcf} --vcf --stats_text --refseq --show_ref_allele --uploaded_allele --use_transcript_ref --variant_class --sift b --polyphen b --gene_phenotype --numbers --hgvs --hgvsg --symbol --canonical --biotype --regulatory --mirna --check_existing --max_af --af_1kg --af_gnomade --exclude_predicted --pick --flag_pick --fork ${Threads} --buffer_size ${buffer_size} ${clinvar_args} ${cosmic_args} --plugin AlphaMissense,file=${alphaMissense} --plugin CADD,snv=${caddSnp},indels=${caddIndel}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
mkdir -p "$(dirname "${vcfDir}")" 2>/dev/null || mkdir -p "${vcfDir}"
mkdir -p "$(dirname "${vepCacheDir}")" 2>/dev/null || mkdir -p "${vepCacheDir}"
mkdir -p "$(dirname "${vepPluginDir}")" 2>/dev/null || mkdir -p "${vepPluginDir}"

eval "$cmd"