#!/bin/bash
# [METADATA]
# TOOL_NAME = fastp
# VERSION = 0.23.4
# THREADS = 1

# Tool Info: fastp (0.23.4)
# Profile: pe_trim_using_singularity

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           Sequence identifier used for file naming"
    echo "  --RawFastqDir     Directory containing the raw input FASTQ files. Suffix = [SeqID]_R1.fastq.gz and [SeqID]_R2.fastq.gz"
    echo "  --TrimFastqDir    Directory where trimmed FASTQ files will be stored. Trimmed Read = Suffix = [SeqID].trimmed_R1.fastq.gz and [SeqID].trimmed_R2.fastq.gz"
    echo "  --qcResDir        Directory where fastp JSON and HTML reports will be stored"
    echo ""
    echo "Optional Parameters:"
    echo "  --singularity_bin Path to singularity executable (Default: singularity)"
    echo "  --fastp_bin       fastp binary name or path inside the container (Default: fastp)"
    echo "  --sif             Path to fastp Singularity image file (Default: /storage/images/fastp-0.23.4.sif)"
    echo "  --bind            Mount paths for Singularity execution (Default: /storage,/data)"
    echo "  --Threads         Number of threads to be used by fastp (Default: 8)"
    echo "  --length_required Reads shorter than this length will be discarded (Default: 100)"
    echo "  --average_qual    Reads with average quality score lower than this will be discarded (Default: 10)"
    echo "  --qualified_quality_phred The quality threshold that a base is considered as qualified (Default: 15)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
RawFastqDir=""
TrimFastqDir=""
qcResDir=""
singularity_bin="singularity"
fastp_bin="fastp"
sif="/storage/images/fastp-0.23.4.sif"
bind="/storage,/data"
Threads="8"
length_required="100"
average_qual="10"
qualified_quality_phred="15"

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --RawFastqDir) RawFastqDir="$2"; shift 2 ;;
        --TrimFastqDir) TrimFastqDir="$2"; shift 2 ;;
        --qcResDir) qcResDir="$2"; shift 2 ;;
        --singularity_bin) singularity_bin="$2"; shift 2 ;;
        --fastp_bin) fastp_bin="$2"; shift 2 ;;
        --sif) sif="$2"; shift 2 ;;
        --bind) bind="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        --length_required) length_required="$2"; shift 2 ;;
        --average_qual) average_qual="$2"; shift 2 ;;
        --qualified_quality_phred) qualified_quality_phred="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${RawFastqDir:-}" ]]; then echo "Error: --RawFastqDir is required"; usage; fi
if [[ -z "${TrimFastqDir:-}" ]]; then echo "Error: --TrimFastqDir is required"; usage; fi
if [[ -z "${qcResDir:-}" ]]; then echo "Error: --qcResDir is required"; usage; fi

# --- [Output Paths] ---
out_read1="${TrimFastqDir}/${SeqID}.trimmed_R1.fastq.gz"
out_read2="${TrimFastqDir}/${SeqID}.trimmed_R2.fastq.gz"
json="${qcResDir}/${SeqID}.fastp.json"
html="${qcResDir}/${SeqID}.fastp.html"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${singularity_bin} exec -B ${bind} ${sif} ${fastp_bin} --thread ${Threads} --in1 ${RawFastqDir}/${SeqID}_R1.fastq.gz --in2 ${RawFastqDir}/${SeqID}_R2.fastq.gz --out1 ${out_read1} --out2 ${out_read2} --json ${json} --html ${html} --trim_poly_g --detect_adapter_for_pe --length_required ${length_required} --average_qual ${average_qual} --qualified_quality_phred ${qualified_quality_phred}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
mkdir -p "$(dirname "${RawFastqDir}")" 2>/dev/null || mkdir -p "${RawFastqDir}"
mkdir -p "$(dirname "${TrimFastqDir}")" 2>/dev/null || mkdir -p "${TrimFastqDir}"
mkdir -p "$(dirname "${qcResDir}")" 2>/dev/null || mkdir -p "${qcResDir}"

eval "$cmd"