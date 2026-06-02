#!/bin/bash
# [METADATA]
# TOOL_NAME = bedtools_genomecov_to_bigwig
# VERSION = v2.27.1
# THREADS = 4

# Tool Info: bedtools_genomecov_to_bigwig (v2.27.1)
# Profile: genome_coverage_and_signal_generation

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           Unique sample identifier used to generate the output file names."
    echo "  --InputFile       Path to the input alignment (BAM) or interval file (BED/GFF/VCF) to compute coverage from."
    echo "  --GenomeSizes     Path to the tab-delimited chromosome sizes reference file."
    echo "  --Extension       Intermediate file extension reflecting bedGraph data structure (must be 'bedgraph' or 'bg')."
    echo "  --OutputDir       Target directory where the generated genome coverage (bedGraph and BigWig) files will be saved."
    echo ""
    echo "Optional Parameters:"
    echo "  --Threads         No description (Default: 4)"
    echo "  --InputFlag       The conditional input flag choice: Use '-ibam' for BAM files or '-i' for interval files. (Default: -ibam)"
    echo "  --GenomeSizesFlag The chromosome layout configuration string for bedtools: Pass '-g [GenomeSizes]' for interval files or '' for BAM files. (Default: )"
    echo "  --ScaleArgs       Pre-calculated scale arguments (e.g., '-scale 0.05 -bg'). Note: '-bg' is required to match bedGraph format specifications. (Default: -bg)"
    echo "  --SortCmd         Downstream shell sort pipeline instruction. Pass '' if sorting is disabled. (Default: | LC_ALL=C sort --buffer-size=8G -k1,1 -k2,2n)"
    echo "  --singularity_bin Absolute path to the singularity/apptainer runtime executable binary. (Default: singularity)"
    echo "  --bedtools_sif    Absolute path to the Bedtools Singularity Image File (.sif). (Default: /storage/images/bedtools-2.27.1.sif)"
    echo "  --bigwig_sif      Absolute path to the Singularity Image File containing the bedGraphToBigWig UCSC utility. (Default: /storage/images/ucsc-bedgraphtobigwig-445.sif)"
    echo "  --bind            Comma-separated host-to-container file system mount paths. (Default: /storage,/data)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
InputFile=""
GenomeSizes=""
Extension=""
OutputDir=""
Threads="4"
InputFlag="-ibam"
GenomeSizesFlag=""
ScaleArgs="-bg"
SortCmd="| LC_ALL=C sort --buffer-size=8G -k1,1 -k2,2n"
singularity_bin="singularity"
bedtools_sif="/storage/images/bedtools-2.27.1.sif"
bigwig_sif="/storage/images/ucsc-bedgraphtobigwig-445.sif"
bind="/storage,/data"

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --InputFile) InputFile="$2"; shift 2 ;;
        --GenomeSizes) GenomeSizes="$2"; shift 2 ;;
        --Extension) Extension="$2"; shift 2 ;;
        --OutputDir) OutputDir="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        --InputFlag) InputFlag="$2"; shift 2 ;;
        --GenomeSizesFlag) GenomeSizesFlag="$2"; shift 2 ;;
        --ScaleArgs) ScaleArgs="$2"; shift 2 ;;
        --SortCmd) SortCmd="$2"; shift 2 ;;
        --singularity_bin) singularity_bin="$2"; shift 2 ;;
        --bedtools_sif) bedtools_sif="$2"; shift 2 ;;
        --bigwig_sif) bigwig_sif="$2"; shift 2 ;;
        --bind) bind="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${InputFile:-}" ]]; then echo "Error: --InputFile is required"; usage; fi
if [[ -z "${GenomeSizes:-}" ]]; then echo "Error: --GenomeSizes is required"; usage; fi
if [[ -z "${Extension:-}" ]]; then echo "Error: --Extension is required"; usage; fi
if [[ -z "${OutputDir:-}" ]]; then echo "Error: --OutputDir is required"; usage; fi

# --- [Output Paths] ---

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${singularity_bin} exec -B ${bind} ${bedtools_sif} bedtools genomecov ${InputFlag} ${InputFile} ${GenomeSizesFlag} ${ScaleArgs} ${SortCmd} --parallel=${Threads} > ${OutputDir}/${SeqID}.${Extension} && ${singularity_bin} exec -B ${bind} ${bigwig_sif} bedGraphToBigWig ${OutputDir}/${SeqID}.${Extension} ${GenomeSizes} ${OutputDir}/${SeqID}.bw"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
if [[ -n "${OutputDir:-}" ]]; then
  if [[ "${OutputDir}" == *.* ]]; then mkdir -p "$(dirname "${OutputDir}")"; else mkdir -p "${OutputDir}"; fi
fi

eval "$cmd"