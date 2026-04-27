#!/bin/bash
# [METADATA]
# TOOL_NAME = CopyKit
# VERSION = 1.0.0
# THREADS = 1

# Tool Info: CopyKit (1.0.0)
# Profile: cbNIPT_CopyNumber_Analysis

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           Original sequence identifier from mapping stage"
    echo "  --NGS_DataBaseDir Directory containing the processed BAM files"
    echo "  --ResultBaseDir   Base directory where all sample results are stored"
    echo ""
    echo "Optional Parameters:"
    echo "  --InputSuffix     Suffix of input BAM (e.g., analysisReady, recal, dedup) (Default: analysisReady)"
    echo "  --mkdir_bin       Path to mkdir executable (Default: mkdir)"
    echo "  --rscript_bin     Path to Rscript executable (Default: Rscript)"
    echo "  --Rscript_path    Path to the CopyKit analysis R script (Default: /storage/home/kangsm/runScripts/cbNIPT/GCX.cbNIPT_Run.CopyKit.Analysis.R)"
    echo "  --Threads         Number of threads for parallel processing (Default: 1)"
    echo "  --BinSize         Bin size for copy number estimation (e.g., 220kb, 500kb) (Default: 220kb)"
    echo "  --GenomeVersion   Reference genome version (hg38/hg19) (Default: hg38)"
    echo "  --SamplePloidy    Expected sample ploidy (Default: 2 for human) (Default: 2)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
NGS_DataBaseDir=""
ResultBaseDir=""
InputSuffix="analysisReady"
mkdir_bin="mkdir"
rscript_bin="Rscript"
Rscript_path="/storage/home/kangsm/runScripts/cbNIPT/GCX.cbNIPT_Run.CopyKit.Analysis.R"
Threads="1"
BinSize="220kb"
GenomeVersion="hg38"
SamplePloidy="2"

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --NGS_DataBaseDir) NGS_DataBaseDir="$2"; shift 2 ;;
        --ResultBaseDir) ResultBaseDir="$2"; shift 2 ;;
        --InputSuffix) InputSuffix="$2"; shift 2 ;;
        --mkdir_bin) mkdir_bin="$2"; shift 2 ;;
        --rscript_bin) rscript_bin="$2"; shift 2 ;;
        --Rscript_path) Rscript_path="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        --BinSize) BinSize="$2"; shift 2 ;;
        --GenomeVersion) GenomeVersion="$2"; shift 2 ;;
        --SamplePloidy) SamplePloidy="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${NGS_DataBaseDir:-}" ]]; then echo "Error: --NGS_DataBaseDir is required"; usage; fi
if [[ -z "${ResultBaseDir:-}" ]]; then echo "Error: --ResultBaseDir is required"; usage; fi

# --- [Output Paths] ---
AnalysisRunDir="${ResultBaseDir}/${SeqID}"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${mkdir_bin} -p ${AnalysisRunDir} && ln -Tsf ${NGS_DataBaseDir}/${SeqID}.${InputSuffix}.bam ${AnalysisRunDir}/${SeqID}.bam && ln -Tsf ${NGS_DataBaseDir}/${SeqID}.${InputSuffix}.bam.bai ${AnalysisRunDir}/${SeqID}.bam.bai && ${rscript_bin} ${Rscript_path}  --SeqID ${SeqID}  --AnalysisRunDir ${AnalysisRunDir}  --BinSize ${BinSize}  --Ploidy ${SamplePloidy}  --GenomeVersion ${GenomeVersion}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
mkdir -p "$(dirname "${NGS_DataBaseDir}")" 2>/dev/null || mkdir -p "${NGS_DataBaseDir}"
mkdir -p "$(dirname "${ResultBaseDir}")" 2>/dev/null || mkdir -p "${ResultBaseDir}"
mkdir -p "$(dirname "${AnalysisRunDir}")" 2>/dev/null || mkdir -p "${AnalysisRunDir}"

eval "$cmd"