#!/bin/bash
# [METADATA]
# TOOL_NAME = Ginkgo_Final_Optimized
# VERSION = 1.0.0
# THREADS = 1

# Tool Info: Ginkgo_Final_Optimized (1.0.0)
# Profile: absolute_path_runner

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           No description"
    echo "  --BinSize         No description"
    echo "  --ReadLength      No description"
    echo ""
    echo "Optional Parameters:"
    echo "  --Threads         No description (Default: 1)"
    echo "  --GinkgoHomeDir   No description (Default: /storage/apps/ginkgo)"
    echo "  --BamToBedDir     No description (Default: /data/cbNIPT/bamToBeds)"
    echo "  --NormalDir       No description (Default: /data/cbNIPT/bamToBeds)"
    echo "  --ResultBaseDir   No description (Default: /data/cbNIPT/ginkgo_analysis)"
    echo "  --Genome          No description (Default: hg38)"
    echo "  --NormalSamples   No description (Default: Normal_01.bed.gz Normal_02.bed.gz Normal_03.bed.gz)"
    echo "  --R_Args_Base     No description (Default: status.xml data 2 [BinMeth] ward euclidean 1 refDummy.bed_mapped 0 ploidyDummy.txt 0 1)"
    echo "  --R_Args_Reclust  No description (Default: status.xml [BinMeth] ward euclidean 0 ploidyDummy.txt 0)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
BinSize=""
ReadLength=""
Threads="1"
GinkgoHomeDir="/storage/apps/ginkgo"
BamToBedDir="/data/cbNIPT/bamToBeds"
NormalDir="/data/cbNIPT/bamToBeds"
ResultBaseDir="/data/cbNIPT/ginkgo_analysis"
Genome="hg38"
NormalSamples="Normal_01.bed.gz Normal_02.bed.gz Normal_03.bed.gz"
R_Args_Base="status.xml data 2 [BinMeth] ward euclidean 1 refDummy.bed_mapped 0 ploidyDummy.txt 0 1"
R_Args_Reclust="status.xml [BinMeth] ward euclidean 0 ploidyDummy.txt 0"

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --BinSize) BinSize="$2"; shift 2 ;;
        --ReadLength) ReadLength="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        --GinkgoHomeDir) GinkgoHomeDir="$2"; shift 2 ;;
        --BamToBedDir) BamToBedDir="$2"; shift 2 ;;
        --NormalDir) NormalDir="$2"; shift 2 ;;
        --ResultBaseDir) ResultBaseDir="$2"; shift 2 ;;
        --Genome) Genome="$2"; shift 2 ;;
        --NormalSamples) NormalSamples="$2"; shift 2 ;;
        --R_Args_Base) R_Args_Base="$2"; shift 2 ;;
        --R_Args_Reclust) R_Args_Reclust="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${BinSize:-}" ]]; then echo "Error: --BinSize is required"; usage; fi
if [[ -z "${ReadLength:-}" ]]; then echo "Error: --ReadLength is required"; usage; fi

# --- [Output Paths] ---
BinMeth="variable_${BinSize}000_${ReadLength}_bwa"
WorkDir="${ResultBaseDir}/${SeqID}/${SeqID}_${BinSize}kb"
BinUnsorted="${GinkgoHomeDir}/scripts/binUnsorted"
BinFile="${GinkgoHomeDir}/genomes/${Genome}/original/${BinMeth}"
BinCount="$(wc -l < ${BinFile})"
ProcessR="${GinkgoHomeDir}/scripts/process.R"
ReclustR="${GinkgoHomeDir}/scripts/reclust.R"
CNVCaller="${GinkgoHomeDir}/scripts/CNVcaller"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="mkdir -p ${WorkDir} && ln -sf ${NormalDir}/*.bed.gz ${WorkDir}/ && ln -sf ${BamToBedDir}/${SeqID}.bed.gz ${WorkDir}/ &&
find ${WorkDir} -maxdepth 1 -name '*.gz' | xargs -I {} bash -c 'N=$(basename '{}' | sed 's/\.bed.*//'); ${BinUnsorted} ${BinFile} ${BinCount} <(zcat '{}' | awk '{if(\$1!~/^chr/) print \'chr\'\$0; else print \$0}') $N '{}'_mapped' &&
find ${WorkDir} -maxdepth 1 -name '*.bed' ! -name '*.gz' | xargs -I {} bash -c 'N=$(basename '{}' | sed 's/\.bed.*//'); ${BinUnsorted} ${BinFile} ${BinCount} <(awk '{if(\$1!~/^chr/) print \'chr\'\$0; else print \$0}' '{}') $N '{}'_mapped' &&
paste ${WorkDir}/*_mapped > ${WorkDir}/data &&
${ProcessR} ${GinkgoHomeDir}/genomes/${Genome}/original ${WorkDir} ${R_Args_Base} && ${ReclustR} ${GinkgoHomeDir}/genomes/${Genome}/original ${WorkDir} ${R_Args_Reclust} && ${CNVCaller} ${WorkDir}/SegCopy ${WorkDir}/CNV1 ${WorkDir}/CNV2"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
mkdir -p "$(dirname "${GinkgoHomeDir}")" 2>/dev/null || mkdir -p "${GinkgoHomeDir}"
mkdir -p "$(dirname "${BamToBedDir}")" 2>/dev/null || mkdir -p "${BamToBedDir}"
mkdir -p "$(dirname "${NormalDir}")" 2>/dev/null || mkdir -p "${NormalDir}"
mkdir -p "$(dirname "${ResultBaseDir}")" 2>/dev/null || mkdir -p "${ResultBaseDir}"
mkdir -p "$(dirname "${WorkDir}")" 2>/dev/null || mkdir -p "${WorkDir}"

eval "$cmd"