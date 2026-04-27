#!/bin/bash
# [METADATA]
# TOOL_NAME = gatk4_filter_mutect_calls
# VERSION = 4.4.0.0
# THREADS = 1

# Tool Info: gatk4_filter_mutect_calls (4.4.0.0)
# Profile: somatic_variant_filtering

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           No description"
    echo "  --vcfDir          입력 및 출력 VCF가 위치한 경로"
    echo "  --qcResDir        Contamination 및 OB-prior 파일 경로"
    echo "  --ReferenceFasta  No description"
    echo "  --TargetInterval  No description"
    echo "  --ContaminationTable [SeqID].contamination.table"
    echo ""
    echo "Optional Parameters:"
    echo "  --Suffix          File name suffix (e.g., 'keep.germline.') (Default: )"
    echo "  --singularity_bin No description (Default: singularity)"
    echo "  --gatk4_sif       No description (Default: /storage/images/gatk-4.4.0.0.sif)"
    echo "  --bind            No description (Default: /storage,/data)"
    echo "  --xmx_mb          16GB memory from original script (Default: 16384)"
    echo "  --Threads         No description (Default: 14)"
    echo "  --extra_args      No description (Default: )"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
vcfDir=""
qcResDir=""
ReferenceFasta=""
TargetInterval=""
ContaminationTable=""
Suffix=""
singularity_bin="singularity"
gatk4_sif="/storage/images/gatk-4.4.0.0.sif"
bind="/storage,/data"
xmx_mb="16384"
Threads="14"
extra_args=""

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --vcfDir) vcfDir="$2"; shift 2 ;;
        --qcResDir) qcResDir="$2"; shift 2 ;;
        --ReferenceFasta) ReferenceFasta="$2"; shift 2 ;;
        --TargetInterval) TargetInterval="$2"; shift 2 ;;
        --ContaminationTable) ContaminationTable="$2"; shift 2 ;;
        --Suffix) Suffix="$2"; shift 2 ;;
        --singularity_bin) singularity_bin="$2"; shift 2 ;;
        --gatk4_sif) gatk4_sif="$2"; shift 2 ;;
        --bind) bind="$2"; shift 2 ;;
        --xmx_mb) xmx_mb="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        --extra_args) extra_args="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${vcfDir:-}" ]]; then echo "Error: --vcfDir is required"; usage; fi
if [[ -z "${qcResDir:-}" ]]; then echo "Error: --qcResDir is required"; usage; fi
if [[ -z "${ReferenceFasta:-}" ]]; then echo "Error: --ReferenceFasta is required"; usage; fi
if [[ -z "${TargetInterval:-}" ]]; then echo "Error: --TargetInterval is required"; usage; fi
if [[ -z "${ContaminationTable:-}" ]]; then echo "Error: --ContaminationTable is required"; usage; fi

# --- [Output Paths] ---
filtered_vcf="${vcfDir}/${SeqID}.mutect2.${Suffix}filtered.vcf"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${singularity_bin} exec -B ${bind} ${gatk4_sif} gatk FilterMutectCalls --java-options '-XX:ParallelGCThreads=${Threads} -Xmx${xmx_mb}m' -V ${vcfDir}/${SeqID}.mutect2.${Suffix}vcf -L ${TargetInterval} --reference ${ReferenceFasta} --contamination-table ${ContaminationTable} --ob-priors ${qcResDir}/${SeqID}.${Suffix}read-orientation-model.tar.gz -O ${filtered_vcf} ${extra_args}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
mkdir -p "$(dirname "${vcfDir}")" 2>/dev/null || mkdir -p "${vcfDir}"
mkdir -p "$(dirname "${qcResDir}")" 2>/dev/null || mkdir -p "${qcResDir}"

eval "$cmd"