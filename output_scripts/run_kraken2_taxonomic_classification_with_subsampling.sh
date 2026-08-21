#!/bin/bash
# [METADATA]
# TOOL_NAME = Kraken2
# VERSION = 2.17.1
# THREADS = 8

# Tool Info: Kraken2 (2.17.1)
# Profile: Taxonomic Classification with Subsampling

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           분석 고유 ID (Sample ID)"
    echo "  --TrimFastqDir    Trimmed FASTQ 파일이 위치한 디렉토리 경로"
    echo "  --OutDir          결과물이 저장될 디렉토리 경로"
    echo ""
    echo "Optional Parameters:"
    echo "  --KrakenDB        Kraken2-db 명령어 (Default: /storage/home/jhkim/Apps/kraken2/bin/kraken2_db)"
    echo "  --sub_r1          Subsampled Read 1 (Default: [OutDir]/[SeqID].trimmed_sub_R1.fastq.gz)"
    echo "  --sub_r2          Subsampled Read 2 (Default: [OutDir]/[SeqID].trimmed_sub_R2.fastq.gz)"
    echo "  --kraken_out      Kraken2 분류 결과 파일 (Default: [OutDir]/[SeqID].kraken)"
    echo "  --report_out      Kraken2 계층적 요약 리포트 (Default: [OutDir]/[SeqID].kreport)"
    echo "  --plot_out        생성된 Krona 대화형 HTML 플롯 (Default: [OutDir]/[SeqID]_krona.html)"
    echo "  --seqkit_bin      Seqkit 실행 절대 경로 (Default: /storage/apps/seqkit-2.8.2/seqkit)"
    echo "  --SampleSeed      난수 시드 (PE 동기화를 위해 R1/R2 동일 적용) (Default: 100)"
    echo "  --SampleNum       추출할 Read 개수 (-n) (Default: 100000)"
    echo "  --R1_suffix       Read 1 입력 파일의 Suffix (Default: .trimmed_R1.fastq.gz)"
    echo "  --R2_suffix       Read 2 입력 파일의 Suffix (Default: .trimmed_R2.fastq.gz)"
    echo "  --kraken2_bin     Kraken2 명령어 (Default: /storage/home/jhkim/Apps/kraken2/bin/kraken2)"
    echo "  --krona_bin       Krona 명령어 (Default: /storage/home/jhkim/Apps/Krona/KronaTools/bin/ktImportTaxonomy)"
    echo "  --Threads         병렬 처리에 사용할 스레드 수 (Default: 8)"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
TrimFastqDir=""
KrakenDB="/storage/home/jhkim/Apps/kraken2/bin/kraken2_db"
OutDir=""
sub_r1="[OutDir]/[SeqID].trimmed_sub_R1.fastq.gz"
sub_r2="[OutDir]/[SeqID].trimmed_sub_R2.fastq.gz"
kraken_out="[OutDir]/[SeqID].kraken"
report_out="[OutDir]/[SeqID].kreport"
plot_out="[OutDir]/[SeqID]_krona.html"
seqkit_bin="/storage/apps/seqkit-2.8.2/seqkit"
SampleSeed="100"
SampleNum="100000"
R1_suffix=".trimmed_R1.fastq.gz"
R2_suffix=".trimmed_R2.fastq.gz"
kraken2_bin="/storage/home/jhkim/Apps/kraken2/bin/kraken2"
krona_bin="/storage/home/jhkim/Apps/Krona/KronaTools/bin/ktImportTaxonomy"
Threads="8"

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --TrimFastqDir) TrimFastqDir="$2"; shift 2 ;;
        --KrakenDB) KrakenDB="$2"; shift 2 ;;
        --OutDir) OutDir="$2"; shift 2 ;;
        --sub_r1) sub_r1="$2"; shift 2 ;;
        --sub_r2) sub_r2="$2"; shift 2 ;;
        --kraken_out) kraken_out="$2"; shift 2 ;;
        --report_out) report_out="$2"; shift 2 ;;
        --plot_out) plot_out="$2"; shift 2 ;;
        --seqkit_bin) seqkit_bin="$2"; shift 2 ;;
        --SampleSeed) SampleSeed="$2"; shift 2 ;;
        --SampleNum) SampleNum="$2"; shift 2 ;;
        --R1_suffix) R1_suffix="$2"; shift 2 ;;
        --R2_suffix) R2_suffix="$2"; shift 2 ;;
        --kraken2_bin) kraken2_bin="$2"; shift 2 ;;
        --krona_bin) krona_bin="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${TrimFastqDir:-}" ]]; then echo "Error: --TrimFastqDir is required"; usage; fi
if [[ -z "${OutDir:-}" ]]; then echo "Error: --OutDir is required"; usage; fi

# --- [Output Paths] ---
sub_r1="${OutDir}/${SeqID}.trimmed_sub_R1.fastq.gz"
sub_r2="${OutDir}/${SeqID}.trimmed_sub_R2.fastq.gz"
kraken_out="${OutDir}/${SeqID}.kraken"
report_out="${OutDir}/${SeqID}.kreport"
plot_out="${OutDir}/${SeqID}_krona.html"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${seqkit_bin} sample -s ${SampleSeed} -n ${SampleNum} ${TrimFastqDir}/${SeqID}${R1_suffix} | gzip > ${sub_r1} && ${seqkit_bin} sample -s ${SampleSeed} -n ${SampleNum} ${TrimFastqDir}/${SeqID}${R2_suffix} | gzip > ${sub_r2} && ${kraken2_bin} --db ${KrakenDB} --threads ${Threads} --report ${report_out} --paired ${sub_r1} ${sub_r2} > ${kraken_out} && ${krona_bin} -m 4 -t 3 ${kraken_out} -o ${plot_out}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
if [[ -n "${sub_r1:-}" ]]; then
  if [[ "${sub_r1}" == *.* ]]; then mkdir -p "$(dirname "${sub_r1}")"; else mkdir -p "${sub_r1}"; fi
fi
if [[ -n "${kraken_out:-}" ]]; then
  if [[ "${kraken_out}" == *.* ]]; then mkdir -p "$(dirname "${kraken_out}")"; else mkdir -p "${kraken_out}"; fi
fi
if [[ -n "${report_out:-}" ]]; then
  if [[ "${report_out}" == *.* ]]; then mkdir -p "$(dirname "${report_out}")"; else mkdir -p "${report_out}"; fi
fi
if [[ -n "${plot_out:-}" ]]; then
  if [[ "${plot_out}" == *.* ]]; then mkdir -p "$(dirname "${plot_out}")"; else mkdir -p "${plot_out}"; fi
fi
if [[ -n "${OutDir:-}" ]]; then
  if [[ "${OutDir}" == *.* ]]; then mkdir -p "$(dirname "${OutDir}")"; else mkdir -p "${OutDir}"; fi
fi
if [[ -n "${TrimFastqDir:-}" ]]; then
  if [[ "${TrimFastqDir}" == *.* ]]; then mkdir -p "$(dirname "${TrimFastqDir}")"; else mkdir -p "${TrimFastqDir}"; fi
fi
if [[ -n "${sub_r2:-}" ]]; then
  if [[ "${sub_r2}" == *.* ]]; then mkdir -p "$(dirname "${sub_r2}")"; else mkdir -p "${sub_r2}"; fi
fi

eval "$cmd"