#!/bin/bash
# [METADATA]
# TOOL_NAME = gatk3
# VERSION = 3.8
# THREADS = 1

# Tool Info: gatk3 (3.8)
# Profile: indel_realign_using_singularity

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           Sequence identifier used for file naming"
    echo "  --BamDir          Directory containing the input BAM and target output BAM"
    echo "  --qcResDir        Directory for the target intervals file output"
    echo "  --ReferenceFasta  Path to the reference genome FASTA file"
    echo "  --KnownIndel1     Path to known Indels VCF (e.g., Mills and 1000G)"
    echo "  --KnownIndel2     Path to additional known Indels VCF"
    echo ""
    echo "Optional Parameters:"
    echo "  --InputSuffix     Suffix of input BAM (e.g., dedup, sorted, primary) (Default: dedup)"
    echo "  --OutputSuffix    Suffix for output BAM (e.g., realign, ir) (Default: realign)"
    echo "  --singularity_bin Path to singularity executable (Default: singularity)"
    echo "  --java_bin        Java executable (inside container) (Default: java)"
    echo "  --gatk_jar        GATK 3.8 jar path inside container (Default: /usr/GenomeAnalysisTK.jar)"
    echo "  --sif             GATK 3.8 singularity image file (Default: /storage/images/gatk-3.8-1.sif)"
    echo "  --bind            Mount paths for singularity (Default: /storage,/data)"
    echo "  --Threads         Number of threads for RealignerTargetCreator (-nt) (Default: 8)"
    echo "  --xmx_mb          Maximum Java heap memory (MB) (Default: 16384)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
BamDir=""
qcResDir=""
ReferenceFasta=""
KnownIndel1=""
KnownIndel2=""
InputSuffix="dedup"
OutputSuffix="realign"
singularity_bin="singularity"
java_bin="java"
gatk_jar="/usr/GenomeAnalysisTK.jar"
sif="/storage/images/gatk-3.8-1.sif"
bind="/storage,/data"
Threads="8"
xmx_mb="16384"

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --BamDir) BamDir="$2"; shift 2 ;;
        --qcResDir) qcResDir="$2"; shift 2 ;;
        --ReferenceFasta) ReferenceFasta="$2"; shift 2 ;;
        --KnownIndel1) KnownIndel1="$2"; shift 2 ;;
        --KnownIndel2) KnownIndel2="$2"; shift 2 ;;
        --InputSuffix) InputSuffix="$2"; shift 2 ;;
        --OutputSuffix) OutputSuffix="$2"; shift 2 ;;
        --singularity_bin) singularity_bin="$2"; shift 2 ;;
        --java_bin) java_bin="$2"; shift 2 ;;
        --gatk_jar) gatk_jar="$2"; shift 2 ;;
        --sif) sif="$2"; shift 2 ;;
        --bind) bind="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        --xmx_mb) xmx_mb="$2"; shift 2 ;;
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
if [[ -z "${KnownIndel1:-}" ]]; then echo "Error: --KnownIndel1 is required"; usage; fi
if [[ -z "${KnownIndel2:-}" ]]; then echo "Error: --KnownIndel2 is required"; usage; fi

# --- [Output Paths] ---
target_intervals="${qcResDir}/${SeqID}.${InputSuffix}.realign_target.intervals"
realigned_bam="${BamDir}/${SeqID}.${OutputSuffix}.bam"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${singularity_bin} exec -B ${bind} ${sif} ${java_bin} -Xmx${xmx_mb}m -jar ${gatk_jar} -T RealignerTargetCreator -R ${ReferenceFasta} -I ${BamDir}/${SeqID}.${InputSuffix}.bam -o ${target_intervals} -known ${KnownIndel1} -known ${KnownIndel2} -nt ${Threads} && ${singularity_bin} exec -B ${bind} ${sif} ${java_bin} -Xmx${xmx_mb}m -jar ${gatk_jar} -T IndelRealigner -R ${ReferenceFasta} -targetIntervals ${target_intervals} -known ${KnownIndel1} -known ${KnownIndel2} -I ${BamDir}/${SeqID}.${InputSuffix}.bam -o ${realigned_bam}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
mkdir -p "$(dirname "${BamDir}")" 2>/dev/null || mkdir -p "${BamDir}"
mkdir -p "$(dirname "${qcResDir}")" 2>/dev/null || mkdir -p "${qcResDir}"

eval "$cmd"