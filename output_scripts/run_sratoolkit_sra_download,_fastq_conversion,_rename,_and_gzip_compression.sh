#!/bin/bash
# [METADATA]
# TOOL_NAME = SRAToolkit
# VERSION = 3.0.10
# THREADS = 8

# Tool Info: SRAToolkit (3.0.10)
# Profile: SRA Download, FASTQ Conversion, Rename, and Gzip Compression

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SraAccession    다운로드할 SRA Run accession"
    echo "  --OutDir          최종 FASTQ.GZ 저장 디렉토리"
    echo ""
    echo "Optional Parameters:"
    echo "  --TmpDir          임시 작업 디렉토리 (Default: [OutDir]/tmp)"
    echo "  --sra_dir         prefetch 다운로드 디렉토리 (Default: [TmpDir]/sra/[SraAccession])"
    echo "  --r1_fastq        Read 1 FASTQ (Default: [OutDir]/[SraAccession]_R1.fastq)"
    echo "  --r2_fastq        Read 2 FASTQ (Default: [OutDir]/[SraAccession]_R2.fastq)"
    echo "  --r1_fastq_gz     압축된 Read 1 FASTQ (Default: [OutDir]/[SraAccession]_R1.fastq.gz)"
    echo "  --r2_fastq_gz     압축된 Read 2 FASTQ (Default: [OutDir]/[SraAccession]_R2.fastq.gz)"
    echo "  --prefetch_bin    prefetch 실행 파일 (Default: /storage/apps/sratoolkit.3.0.10/bin/prefetch)"
    echo "  --fasterq_dump_bin fasterq-dump 실행 파일 (Default: /storage/apps/sratoolkit.3.0.10/bin/fasterq-dump)"
    echo "  --Threads         사용할 스레드 수 (Default: 8)"
    echo "  --MaxSize         prefetch 최대 다운로드 크기 (Default: 100G)"
    echo "  --SplitMode       paired-end FASTQ 분리 옵션 (Default: --split-files)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SraAccession=""
OutDir=""
TmpDir="[OutDir]/tmp"
sra_dir="[TmpDir]/sra/[SraAccession]"
r1_fastq="[OutDir]/[SraAccession]_R1.fastq"
r2_fastq="[OutDir]/[SraAccession]_R2.fastq"
r1_fastq_gz="[OutDir]/[SraAccession]_R1.fastq.gz"
r2_fastq_gz="[OutDir]/[SraAccession]_R2.fastq.gz"
prefetch_bin="/storage/apps/sratoolkit.3.0.10/bin/prefetch"
fasterq_dump_bin="/storage/apps/sratoolkit.3.0.10/bin/fasterq-dump"
Threads="8"
MaxSize="100G"
SplitMode="--split-files"

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SraAccession) SraAccession="$2"; shift 2 ;;
        --OutDir) OutDir="$2"; shift 2 ;;
        --TmpDir) TmpDir="$2"; shift 2 ;;
        --sra_dir) sra_dir="$2"; shift 2 ;;
        --r1_fastq) r1_fastq="$2"; shift 2 ;;
        --r2_fastq) r2_fastq="$2"; shift 2 ;;
        --r1_fastq_gz) r1_fastq_gz="$2"; shift 2 ;;
        --r2_fastq_gz) r2_fastq_gz="$2"; shift 2 ;;
        --prefetch_bin) prefetch_bin="$2"; shift 2 ;;
        --fasterq_dump_bin) fasterq_dump_bin="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        --MaxSize) MaxSize="$2"; shift 2 ;;
        --SplitMode) SplitMode="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SraAccession:-}" ]]; then echo "Error: --SraAccession is required"; usage; fi
if [[ -z "${OutDir:-}" ]]; then echo "Error: --OutDir is required"; usage; fi

# --- [Output Paths] ---
sra_dir="${TmpDir}/sra/${SraAccession}"
r1_fastq="${OutDir}/${SraAccession}_R1.fastq"
r2_fastq="${OutDir}/${SraAccession}_R2.fastq"
r1_fastq_gz="${OutDir}/${SraAccession}_R1.fastq.gz"
r2_fastq_gz="${OutDir}/${SraAccession}_R2.fastq.gz"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="mkdir -p ${OutDir} ${TmpDir}/sra ${TmpDir}/fasterq && ${prefetch_bin} --max-size ${MaxSize} --output-directory ${TmpDir}/sra ${SraAccession} && ${fasterq_dump_bin} ${SplitMode} --threads ${Threads} --temp ${TmpDir}/fasterq --outdir ${OutDir} ${sra_dir} && mv ${OutDir}/${SraAccession}_1.fastq ${r1_fastq} && mv ${OutDir}/${SraAccession}_2.fastq ${r2_fastq} && gzip -f ${r1_fastq} && gzip -f ${r2_fastq}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
if [[ -n "${r1_fastq:-}" ]]; then
  if [[ "${r1_fastq}" == *.* ]]; then mkdir -p "$(dirname "${r1_fastq}")"; else mkdir -p "${r1_fastq}"; fi
fi
if [[ -n "${OutDir:-}" ]]; then
  if [[ "${OutDir}" == *.* ]]; then mkdir -p "$(dirname "${OutDir}")"; else mkdir -p "${OutDir}"; fi
fi
if [[ -n "${sra_dir:-}" ]]; then
  if [[ "${sra_dir}" == *.* ]]; then mkdir -p "$(dirname "${sra_dir}")"; else mkdir -p "${sra_dir}"; fi
fi
if [[ -n "${TmpDir:-}" ]]; then
  if [[ "${TmpDir}" == *.* ]]; then mkdir -p "$(dirname "${TmpDir}")"; else mkdir -p "${TmpDir}"; fi
fi
if [[ -n "${r2_fastq_gz:-}" ]]; then
  if [[ "${r2_fastq_gz}" == *.* ]]; then mkdir -p "$(dirname "${r2_fastq_gz}")"; else mkdir -p "${r2_fastq_gz}"; fi
fi
if [[ -n "${r1_fastq_gz:-}" ]]; then
  if [[ "${r1_fastq_gz}" == *.* ]]; then mkdir -p "$(dirname "${r1_fastq_gz}")"; else mkdir -p "${r1_fastq_gz}"; fi
fi
if [[ -n "${r2_fastq:-}" ]]; then
  if [[ "${r2_fastq}" == *.* ]]; then mkdir -p "$(dirname "${r2_fastq}")"; else mkdir -p "${r2_fastq}"; fi
fi

eval "$cmd"