#!/bin/bash
# [METADATA]
# TOOL_NAME = gatk4_contamination_estimation
# VERSION = 4.4.0.0
# THREADS = 1

# Tool Info: gatk4_contamination_estimation (4.4.0.0)
# Profile: contamination_calculation

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           No description"
    echo "  --BamDir          No description"
    echo "  --qcResDir        결과 테이블(.table) 저장 경로"
    echo "  --ReferenceFasta  No description"
    echo "  --TargetInterval  분석 대상 영역 (Interval_list)"
    echo "  --VcfGnomad       gnomAD germline resource VCF (Required for pileup)"
    echo ""
    echo "Optional Parameters:"
    echo "  --InputSuffix     No description (Default: analysisReady)"
    echo "  --singularity_bin No description (Default: singularity)"
    echo "  --gatk4_sif       No description (Default: /storage/images/gatk-4.4.0.0.sif)"
    echo "  --bind            No description (Default: /storage,/data)"
    echo "  --Threads         No description (Default: 14)"
    echo "  --pileup_xmx      Xmx for GetPileupSummaries (8G) (Default: 8192)"
    echo "  --contam_xmx      Xmx for CalculateContamination (16G) (Default: 16384)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
BamDir=""
qcResDir=""
ReferenceFasta=""
TargetInterval=""
VcfGnomad=""
InputSuffix="analysisReady"
singularity_bin="singularity"
gatk4_sif="/storage/images/gatk-4.4.0.0.sif"
bind="/storage,/data"
Threads="14"
pileup_xmx="8192"
contam_xmx="16384"

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --BamDir) BamDir="$2"; shift 2 ;;
        --qcResDir) qcResDir="$2"; shift 2 ;;
        --ReferenceFasta) ReferenceFasta="$2"; shift 2 ;;
        --TargetInterval) TargetInterval="$2"; shift 2 ;;
        --VcfGnomad) VcfGnomad="$2"; shift 2 ;;
        --InputSuffix) InputSuffix="$2"; shift 2 ;;
        --singularity_bin) singularity_bin="$2"; shift 2 ;;
        --gatk4_sif) gatk4_sif="$2"; shift 2 ;;
        --bind) bind="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        --pileup_xmx) pileup_xmx="$2"; shift 2 ;;
        --contam_xmx) contam_xmx="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${BamDir:-}" ]]; then echo "Error: --BamDir is required"; usage; fi
if [[ -z "${qcResDir:-}" ]]; then echo "Error: --qcResDir is required"; usage; fi
if [[ -z "${ReferenceFasta:-}" ]]; then echo "Error: --ReferenceFasta is required"; usage; fi
if [[ -z "${TargetInterval:-}" ]]; then echo "Error: --TargetInterval is required"; usage; fi
if [[ -z "${VcfGnomad:-}" ]]; then echo "Error: --VcfGnomad is required"; usage; fi

# --- [Output Paths] ---
pileup_table="${qcResDir}/${SeqID}.targeted_sequencing.table"
contamination_table="${qcResDir}/${SeqID}.contamination.table"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${singularity_bin} exec -B ${bind} ${gatk4_sif} gatk GetPileupSummaries --java-options '-XX:ParallelGCThreads=${Threads} -Xmx${pileup_xmx}m' -I ${BamDir}/${SeqID}.${InputSuffix}.bam -V ${VcfGnomad} -L ${TargetInterval} -R ${ReferenceFasta} -O ${pileup_table} && ${singularity_bin} exec -B ${bind} ${gatk4_sif} gatk CalculateContamination --java-options '-XX:ParallelGCThreads=${Threads} -Xmx${contam_xmx}m' -I ${pileup_table} -O ${contamination_table}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
mkdir -p "$(dirname "${BamDir}")" 2>/dev/null || mkdir -p "${BamDir}"
mkdir -p "$(dirname "${qcResDir}")" 2>/dev/null || mkdir -p "${qcResDir}"

eval "$cmd"