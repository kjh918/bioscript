#!/bin/bash
# [METADATA]
# TOOL_NAME = alfred
# VERSION = 0.2.6
# THREADS = 1

# Tool Info: alfred (0.2.6)
# Profile: bam_qc_extraction

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           No description"
    echo "  --BamDir          No description"
    echo "  --qcResDir        No description"
    echo "  --ReferenceFasta  No description"
    echo "  --TargetBed       Target BED file for WES metrics"
    echo "  --InputSuffix     No description"
    echo ""
    echo "Optional Parameters:"
    echo "  --singularity_bin No description (Default: singularity)"
    echo "  --alfred_sif      No description (Default: /storage/images/alfred-0.2.6.sif)"
    echo "  --alfred_bin      No description (Default: /opt/alfred/bin/alfred)"
    echo "  --bind            No description (Default: /storage,/data)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
BamDir=""
qcResDir=""
ReferenceFasta=""
TargetBed=""
InputSuffix="analysisReady"
singularity_bin="singularity"
alfred_sif="/storage/images/alfred-0.2.6.sif"
alfred_bin="/opt/alfred/bin/alfred"
bind="/storage,/data"

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --BamDir) BamDir="$2"; shift 2 ;;
        --qcResDir) qcResDir="$2"; shift 2 ;;
        --ReferenceFasta) ReferenceFasta="$2"; shift 2 ;;
        --TargetBed) TargetBed="$2"; shift 2 ;;
        --InputSuffix) InputSuffix="$2"; shift 2 ;;
        --singularity_bin) singularity_bin="$2"; shift 2 ;;
        --alfred_sif) alfred_sif="$2"; shift 2 ;;
        --alfred_bin) alfred_bin="$2"; shift 2 ;;
        --bind) bind="$2"; shift 2 ;;
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
if [[ -z "${TargetBed:-}" ]]; then echo "Error: --TargetBed is required"; usage; fi
if [[ -z "${InputSuffix:-}" ]]; then echo "Error: --InputSuffix is required"; usage; fi

# --- [Output Paths] ---
alfred_raw_tsv="${qcResDir}/${SeqID}.alfred.qc.tsv.gz"
chr_map_stats="${qcResDir}/${SeqID}.alfred.chr.map.stats.txt"
target_coverage="${qcResDir}/${SeqID}.alfred.target.coverage.txt"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${singularity_bin} exec -B ${bind} ${alfred_sif} ${alfred_bin} qc --reference ${ReferenceFasta} --bed ${TargetBed} --outfile ${alfred_raw_tsv} ${BamDir}/${SeqID}.${InputSuffix}.bam && zgrep '^CM' ${alfred_raw_tsv} | cut -f 2- > ${chr_map_stats} && zgrep '^TC' ${alfred_raw_tsv} | cut -f 2- > ${target_coverage}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
mkdir -p "$(dirname "${BamDir}")" 2>/dev/null || mkdir -p "${BamDir}"
mkdir -p "$(dirname "${qcResDir}")" 2>/dev/null || mkdir -p "${qcResDir}"

eval "$cmd"