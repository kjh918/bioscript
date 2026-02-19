#!/bin/bash
# [METADATA]
# TOOL_NAME = gatk4
# VERSION = 4.4.0.0
# THREADS = 1

# Tool Info: gatk4 (4.4.0.0)
# Profile: HaplotypeCaller_based_on_chrom

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           Sequence identifier for sample tracking"
    echo "  --Chromosome      Target chromosome or interval (L option)"
    echo "  --BamDir          Directory containing the input BAM file"
    echo "  --ResultDir       Output directory for VCF and GVCF files"
    echo ""
    echo "Optional Parameters:"
    echo "  --ReferenceFasta  Reference genome FASTA (Default: /storage/references_and_index/hg38/fasta/cbNIPT/hg38.fa)"
    echo "  --DbsnpVcf        Known variation database (dbSNP) (Default: /storage/references_and_index/hg38/vcf/Homo_sapiens_assembly38.dbsnp138.vcf.gz)"
    echo "  --InputSuffix     Suffix of input BAM (e.g., analysisReady, recal, dedup) (Default: analysisReady)"
    echo "  --singularity_bin Path to singularity executable (Default: singularity)"
    echo "  --gatk_bin        GATK binary path inside container (Default: gatk)"
    echo "  --sif             GATK singularity image file (Default: /storage/images/gatk-4.4.0.0.sif)"
    echo "  --bind            Mount paths for singularity (Default: /storage,/data)"
    echo "  --Threads         Threads for ParallelGC (Default: 4)"
    echo "  --Memory          Java Xmx memory setting (e.g., 32g) (Default: 32g)"
    echo "  --Ploidy          Sample ploidy (Default: 2) (Default: 2)"
    echo "  --TmpDir          Temporary directory for GATK (Default: [ResultDir]/tmp/[SeqID]_[Chromosome])"
    echo "  --IncludeNonVariant Whether to include non-variant sites in GenotypeGVCFs (Default: false)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
Chromosome=""
BamDir=""
ResultDir=""
ReferenceFasta="/storage/references_and_index/hg38/fasta/cbNIPT/hg38.fa"
DbsnpVcf="/storage/references_and_index/hg38/vcf/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
InputSuffix="analysisReady"
singularity_bin="singularity"
gatk_bin="gatk"
sif="/storage/images/gatk-4.4.0.0.sif"
bind="/storage,/data"
Threads="4"
Memory="32g"
Ploidy="2"
TmpDir="[ResultDir]/tmp/[SeqID]_[Chromosome]"
IncludeNonVariant="false"

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --Chromosome) Chromosome="$2"; shift 2 ;;
        --BamDir) BamDir="$2"; shift 2 ;;
        --ResultDir) ResultDir="$2"; shift 2 ;;
        --ReferenceFasta) ReferenceFasta="$2"; shift 2 ;;
        --DbsnpVcf) DbsnpVcf="$2"; shift 2 ;;
        --InputSuffix) InputSuffix="$2"; shift 2 ;;
        --singularity_bin) singularity_bin="$2"; shift 2 ;;
        --gatk_bin) gatk_bin="$2"; shift 2 ;;
        --sif) sif="$2"; shift 2 ;;
        --bind) bind="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        --Memory) Memory="$2"; shift 2 ;;
        --Ploidy) Ploidy="$2"; shift 2 ;;
        --TmpDir) TmpDir="$2"; shift 2 ;;
        --IncludeNonVariant) IncludeNonVariant="$2"; shift 2 ;;
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

# --- [Output Paths] ---
OutGvcf="${ResultDir}/${SeqID}.${Chromosome}.gvcf.gz"
OutVcf="${ResultDir}/${SeqID}.${Chromosome}.vcf.gz"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${singularity_bin} exec -B ${bind} ${sif} ${gatk_bin} --java-options '-XX:ParallelGCThreads=${Threads} -Xmx${Memory} -Djava.io.tmpdir=${TmpDir}'  HaplotypeCaller  -R ${ReferenceFasta}  -I ${BamDir}/${SeqID}.${InputSuffix}.bam  -L ${Chromosome}  -ploidy ${Ploidy}  -stand-call-conf 30  --dbsnp ${DbsnpVcf}  -O ${OutGvcf}  -ERC GVCF && 
${singularity_bin} exec -B ${bind} ${sif} ${gatk_bin} --java-options '-XX:ParallelGCThreads=${Threads} -Xmx${Memory} -Djava.io.tmpdir=${TmpDir}'  GenotypeGVCFs  --include-non-variant-sites ${IncludeNonVariant}  -R ${ReferenceFasta}  -V ${OutGvcf}  -O ${OutVcf}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
mkdir -p "$(dirname "${BamDir}")" 2>/dev/null || mkdir -p "${BamDir}"
mkdir -p "$(dirname "${ResultDir}")" 2>/dev/null || mkdir -p "${ResultDir}"
mkdir -p "$(dirname "${TmpDir}")" 2>/dev/null || mkdir -p "${TmpDir}"

eval "$cmd"