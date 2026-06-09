#!/bin/bash
# [METADATA]
# TOOL_NAME = cbNIPT_MakeBins
# VERSION = 1.0.0
# THREADS = 1

# Tool Info: cbNIPT_MakeBins (1.0.0)
# Profile: Generate reference genome bins

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --ReferenceFasta  참조 유전체 FASTA 파일 경로"
    echo "  --MappabilityBW   Mappability BigWig 파일"
    echo "  --OutBinFile      생성된 빈 정보를 저장할 경로 (예: hg38_100kb_bins.bed.gz)"
    echo "  --script_path     Path to the CNV analyzer script"
    echo ""
    echo "Optional Parameters:"
    echo "  --BinSize         Bin 크기 (100kb) (Default: 100000)"
    echo "  --MinMappability  최소 Mappability 점수 (Default: 0.9)"
    echo "  --MinGC           최소 GC 함량 필터 (Default: 0.3)"
    echo "  --MaxGC           최대 GC 함량 필터 (Default: 0.7)"
    echo "  --IncludeSexChrom 성염색체 포함 여부 플래그 (제외 시 빈 문자열) (Default: --IncludeSexChrom)"
    echo "  --python_bin      Python executable (Default: python)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
ReferenceFasta=""
MappabilityBW=""
OutBinFile=""
BinSize="100000"
MinMappability="0.9"
MinGC="0.3"
MaxGC="0.7"
IncludeSexChrom="--IncludeSexChrom"
python_bin="python"
script_path=""

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --ReferenceFasta) ReferenceFasta="$2"; shift 2 ;;
        --MappabilityBW) MappabilityBW="$2"; shift 2 ;;
        --OutBinFile) OutBinFile="$2"; shift 2 ;;
        --BinSize) BinSize="$2"; shift 2 ;;
        --MinMappability) MinMappability="$2"; shift 2 ;;
        --MinGC) MinGC="$2"; shift 2 ;;
        --MaxGC) MaxGC="$2"; shift 2 ;;
        --IncludeSexChrom) IncludeSexChrom="$2"; shift 2 ;;
        --python_bin) python_bin="$2"; shift 2 ;;
        --script_path) script_path="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${ReferenceFasta:-}" ]]; then echo "Error: --ReferenceFasta is required"; usage; fi
if [[ -z "${MappabilityBW:-}" ]]; then echo "Error: --MappabilityBW is required"; usage; fi
if [[ -z "${OutBinFile:-}" ]]; then echo "Error: --OutBinFile is required"; usage; fi
if [[ -z "${script_path:-}" ]]; then echo "Error: --script_path is required"; usage; fi

# --- [Output Paths] ---
OutBinFile=""

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${python_bin} ${script_path} make-bins --ReferenceFasta ${ReferenceFasta} --MappabilityBW ${MappabilityBW} --OutBinFile ${OutBinFile} --BinSize ${BinSize} --MinMappability ${MinMappability} --MinGC ${MinGC} --MaxGC ${MaxGC} ${IncludeSexChrom}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
if [[ -n "${OutBinFile:-}" ]]; then
  if [[ "${OutBinFile}" == *.* ]]; then mkdir -p "$(dirname "${OutBinFile}")"; else mkdir -p "${OutBinFile}"; fi
fi

eval "$cmd"