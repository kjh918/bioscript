#!/bin/bash
# [METADATA]
# TOOL_NAME = sob_detector
# VERSION = 1.0.4
# THREADS = 1

# Tool Info: sob_detector (1.0.4)
# Profile: somatic_artifact_cleanup

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           No description"
    echo "  --vcfDir          필터링할 VCF가 위치한 경로"
    echo "  --BamDir          분석용 BAM(analysisReady) 경로"
    echo ""
    echo "Optional Parameters:"
    echo "  --bias_filtered_vcf SOBDetector에 의해 아티팩트가 제거된 최종 VCF (Default: [vcfDir]/[SeqID].mutect2.bias.filtered.vcf)"
    echo "  --InputVcfSuffix  No description (Default: filtered)"
    echo "  --InputBamSuffix  No description (Default: analysisReady)"
    echo "  --singularity_bin No description (Default: singularity)"
    echo "  --sob_sif         No description (Default: /storage/images/sobdetector.sif)"
    echo "  --sob_bin         No description (Default: SOBDetector)"
    echo "  --bind            No description (Default: /storage,/data)"
    echo "  --minBaseQ        No description (Default: 20)"
    echo "  --minMapQ         No description (Default: 20)"
    echo "  --only_passed     PASS 필터된 변이만 처리할지 여부 (Default: false)"
    echo "  --extra_args      No description (Default: )"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
vcfDir=""
BamDir=""
bias_filtered_vcf="[vcfDir]/[SeqID].mutect2.bias.filtered.vcf"
InputVcfSuffix="filtered"
InputBamSuffix="analysisReady"
singularity_bin="singularity"
sob_sif="/storage/images/sobdetector.sif"
sob_bin="SOBDetector"
bind="/storage,/data"
minBaseQ="20"
minMapQ="20"
only_passed="false"
extra_args=""

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --vcfDir) vcfDir="$2"; shift 2 ;;
        --BamDir) BamDir="$2"; shift 2 ;;
        --bias_filtered_vcf) bias_filtered_vcf="$2"; shift 2 ;;
        --InputVcfSuffix) InputVcfSuffix="$2"; shift 2 ;;
        --InputBamSuffix) InputBamSuffix="$2"; shift 2 ;;
        --singularity_bin) singularity_bin="$2"; shift 2 ;;
        --sob_sif) sob_sif="$2"; shift 2 ;;
        --sob_bin) sob_bin="$2"; shift 2 ;;
        --bind) bind="$2"; shift 2 ;;
        --minBaseQ) minBaseQ="$2"; shift 2 ;;
        --minMapQ) minMapQ="$2"; shift 2 ;;
        --only_passed) only_passed="$2"; shift 2 ;;
        --extra_args) extra_args="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${vcfDir:-}" ]]; then echo "Error: --vcfDir is required"; usage; fi
if [[ -z "${BamDir:-}" ]]; then echo "Error: --BamDir is required"; usage; fi

# --- [Output Paths] ---
bias_filtered_vcf="${vcfDir}/${SeqID}.mutect2.bias.filtered.vcf"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${singularity_bin} exec -B ${bind} ${sob_sif} ${sob_bin}  --input-type VCF --input-variants ${vcfDir}/${SeqID}.mutect2.${InputVcfSuffix}.vcf --input-bam ${BamDir}/${SeqID}.${InputBamSuffix}.bam --output-variants ${bias_filtered_vcf} --minBaseQuality ${minBaseQ} --minMappingQuality ${minMapQ} --only-passed ${only_passed} ${extra_args}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
if [[ -n "${bias_filtered_vcf:-}" ]]; then
  if [[ "${bias_filtered_vcf}" == *.* ]]; then mkdir -p "$(dirname "${bias_filtered_vcf}")"; else mkdir -p "${bias_filtered_vcf}"; fi
fi
if [[ -n "${vcfDir:-}" ]]; then
  if [[ "${vcfDir}" == *.* ]]; then mkdir -p "$(dirname "${vcfDir}")"; else mkdir -p "${vcfDir}"; fi
fi
if [[ -n "${BamDir:-}" ]]; then
  if [[ "${BamDir}" == *.* ]]; then mkdir -p "$(dirname "${BamDir}")"; else mkdir -p "${BamDir}"; fi
fi

eval "$cmd"