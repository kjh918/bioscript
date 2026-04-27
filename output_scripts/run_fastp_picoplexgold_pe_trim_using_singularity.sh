#!/bin/bash
# [METADATA]
# TOOL_NAME = fastp
# VERSION = 0.23.4
# THREADS = 1

# Tool Info: fastp (0.23.4)
# Profile: PicoplexGold_pe_trim_using_singularity

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           Sequence identifier used for file naming and sample tracking"
    echo "  --RawFastqDir     Directory containing the raw input FASTQ files"
    echo "  --TrimFastqDir    Target directory for the processed (trimmed) FASTQ files"
    echo "  --qcResDir        Directory for quality control reports (JSON/HTML)"
    echo ""
    echo "Optional Parameters:"
    echo "  --singularity_bin Path to the singularity executable (Default: singularity)"
    echo "  --fastp_bin       fastp executable name or path inside container (Default: fastp)"
    echo "  --sif             Path to fastp Singularity image (.sif) (Default: /storage/images/fastp-0.23.4.sif)"
    echo "  --bind            Directories to mount for Singularity (Default: /storage,/data)"
    echo "  --Threads         Number of threads for parallel processing (Default: 8)"
    echo "  --length_required Minimum read length filter (Default: 100)"
    echo "  --average_qual    Minimum average quality filter (Default: 10)"
    echo "  --qualified_quality_phred Phred quality score threshold (base quality) (Default: 15)"
    echo "  --trim_front1     Number of bases to trim from the start of Read 1 (Default: 14)"
    echo "  --trim_front2     Number of bases to trim from the start of Read 2 (Default: 14)"
    echo "  --adapter_sequence Adapter sequence for Read 1 (Default: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA)"
    echo "  --adapter_sequence_r2 Adapter sequence for Read 2 (Default: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT)"
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
trim_front1="14"
trim_front2="14"
adapter_sequence="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
adapter_sequence_r2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

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
        --trim_front1) trim_front1="$2"; shift 2 ;;
        --trim_front2) trim_front2="$2"; shift 2 ;;
        --adapter_sequence) adapter_sequence="$2"; shift 2 ;;
        --adapter_sequence_r2) adapter_sequence_r2="$2"; shift 2 ;;
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
cmd="${singularity_bin} exec -B ${bind} ${sif} ${fastp_bin} --thread ${Threads} --in1 ${RawFastqDir}/${SeqID}_R1.fastq.gz --in2 ${RawFastqDir}/${SeqID}_R2.fastq.gz --out1 ${out_read1} --out2 ${out_read2} --json ${json} --html ${html} --trim_poly_g --detect_adapter_for_pe --adapter_sequence ${adapter_sequence} --adapter_sequence_r2 ${adapter_sequence_r2} --length_required ${length_required} --average_qual ${average_qual} --qualified_quality_phred ${qualified_quality_phred} --trim_front1 ${trim_front1} --trim_front2 ${trim_front2}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
mkdir -p "$(dirname "${RawFastqDir}")" 2>/dev/null || mkdir -p "${RawFastqDir}"
mkdir -p "$(dirname "${TrimFastqDir}")" 2>/dev/null || mkdir -p "${TrimFastqDir}"
mkdir -p "$(dirname "${qcResDir}")" 2>/dev/null || mkdir -p "${qcResDir}"

eval "$cmd"