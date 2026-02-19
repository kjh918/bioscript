#!/bin/bash
# [METADATA]
# TOOL_NAME = fastqc
# VERSION = 0.12.1
# THREADS = 8

# Tool Info: fastqc (0.12.1)
# Profile: singularity_pe_qc_extract

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           Sequence identifier used to locate FASTQ files"
    echo "  --RawFastqDir     Directory containing raw or trimmed FASTQ files"
    echo "  --qcResDir        Directory where FastQC reports will be saved"
    echo ""
    echo "Optional Parameters:"
    echo "  --singularity_bin Path to singularity executable (Default: singularity)"
    echo "  --fastqc_bin      FastQC binary name inside container (Default: fastqc)"
    echo "  --sif             Path to FastQC SIF image (Default: /storage/images/fastqc-0.12.1.sif)"
    echo "  --bind            Mount paths for singularity exec (Default: /storage,/data)"
    echo "  --threads         Number of threads (FastQC processes files in parallel) (Default: 8)"
    echo "  --fastqc_args     Additional FastQC flags (e.g., --extract, --nogroup) (Default: --extract)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
RawFastqDir=""
qcResDir=""
singularity_bin="singularity"
fastqc_bin="fastqc"
sif="/storage/images/fastqc-0.12.1.sif"
bind="/storage,/data"
threads="8"
fastqc_args="--extract"

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --RawFastqDir) RawFastqDir="$2"; shift 2 ;;
        --qcResDir) qcResDir="$2"; shift 2 ;;
        --singularity_bin) singularity_bin="$2"; shift 2 ;;
        --fastqc_bin) fastqc_bin="$2"; shift 2 ;;
        --sif) sif="$2"; shift 2 ;;
        --bind) bind="$2"; shift 2 ;;
        --threads) threads="$2"; shift 2 ;;
        --fastqc_args) fastqc_args="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${RawFastqDir:-}" ]]; then echo "Error: --RawFastqDir is required"; usage; fi
if [[ -z "${qcResDir:-}" ]]; then echo "Error: --qcResDir is required"; usage; fi

# --- [Output Paths] ---
r1_html="${qcResDir}/${SeqID}_R1_fastqc.html"
r2_html="${qcResDir}/${SeqID}_R2_fastqc.html"
r1_zip="${qcResDir}/${SeqID}_R1_fastqc.zip"
r2_zip="${qcResDir}/${SeqID}_R2_fastqc.zip"
r1_dir="${qcResDir}/${SeqID}_R1_fastqc"
r2_dir="${qcResDir}/${SeqID}_R2_fastqc"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${singularity_bin} exec -B ${bind} ${sif} ${fastqc_bin} ${fastqc_args} --threads ${threads} --outdir ${qcResDir} ${RawFastqDir}/${SeqID}_R1.fastq.gz ${RawFastqDir}/${SeqID}_R2.fastq.gz"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
mkdir -p "$(dirname "${RawFastqDir}")" 2>/dev/null || mkdir -p "${RawFastqDir}"
mkdir -p "$(dirname "${qcResDir}")" 2>/dev/null || mkdir -p "${qcResDir}"

eval "$cmd"