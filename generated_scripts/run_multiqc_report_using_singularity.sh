#!/bin/bash
# [METADATA]
# TOOL_NAME = multiqc
# VERSION = 1.16
# THREADS = 1

# Tool Info: multiqc (1.16)
# Profile: report_using_singularity

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           Sequence identifier used for the report filename"
    echo "  --qcResDir        Directory containing all QC data and target output directory"
    echo ""
    echo "Optional Parameters:"
    echo "  --singularity_bin Path to singularity executable (Default: singularity)"
    echo "  --multiqc_bin     MultiQC binary name inside container (Default: multiqc)"
    echo "  --sif             MultiQC SIF image path (Default: /storage/images/multiqc-1.16.sif)"
    echo "  --bind            Mount paths for singularity (Default: /storage,/data)"
    echo "  --mqc_config      Path to custom MultiQC config file (Default: /storage/home/kangsm/runScripts/NGS_config.MultiQC_Custom.yaml)"
    echo "  --mqc_filename    Base name for the output report file (Default: [SeqID].QC.Results)"
    echo "  --mqc_args        Additional MultiQC flags (Default: --force --data-dir)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
qcResDir=""
singularity_bin="singularity"
multiqc_bin="multiqc"
sif="/storage/images/multiqc-1.16.sif"
bind="/storage,/data"
mqc_config="/storage/home/kangsm/runScripts/NGS_config.MultiQC_Custom.yaml"
mqc_filename="[SeqID].QC.Results"
mqc_args="--force --data-dir"

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --qcResDir) qcResDir="$2"; shift 2 ;;
        --singularity_bin) singularity_bin="$2"; shift 2 ;;
        --multiqc_bin) multiqc_bin="$2"; shift 2 ;;
        --sif) sif="$2"; shift 2 ;;
        --bind) bind="$2"; shift 2 ;;
        --mqc_config) mqc_config="$2"; shift 2 ;;
        --mqc_filename) mqc_filename="$2"; shift 2 ;;
        --mqc_args) mqc_args="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${qcResDir:-}" ]]; then echo "Error: --qcResDir is required"; usage; fi

# --- [Output Paths] ---
report_html="${qcResDir}/${SeqID}.QC.Results.html"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${singularity_bin} exec -B ${bind} ${sif} ${multiqc_bin} ${mqc_args} --filename ${mqc_filename} --outdir ${qcResDir} --config ${mqc_config} ${qcResDir}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
mkdir -p "$(dirname "${qcResDir}")" 2>/dev/null || mkdir -p "${qcResDir}"

eval "$cmd"