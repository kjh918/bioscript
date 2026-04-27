#!/bin/bash
# [METADATA]
# TOOL_NAME = gatk4_mutect2
# VERSION = 4.4.0.0
# THREADS = 1

# Tool Info: gatk4_mutect2 (4.4.0.0)
# Profile: tumor_only_with_f1r2

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Inputs:"
    echo "  --SeqID           No description"
    echo "  --BamDir          No description"
    echo "  --vcfDir          Raw VCF 결과 저장 경로"
    echo "  --qcResDir        f1r2.tar.gz 및 중간 파일 저장 경로"
    echo "  --ReferenceFasta  No description"
    echo "  --TargetInterval  Target BED 또는 interval_list"
    echo "  --VcfGnomad       gnomAD germline resource VCF"
    echo "  --VcfPon          Panel of Normals VCF"
    echo ""
    echo "Optional Parameters:"
    echo "  --InputSuffix     No description (Default: analysisReady)"
    echo "  --singularity_bin No description (Default: singularity)"
    echo "  --gatk4_sif       No description (Default: /storage/images/gatk-4.4.0.0.sif)"
    echo "  --bind            No description (Default: /storage,/data)"
    echo "  --xmx_mb          Java Max Heap Size (MB) (Default: 32768)"
    echo "  --Threads         Parallel GC 및 OpenMP 스레드 수 (Default: 14)"
    echo "  --pairHMM         No description (Default: AVX_LOGLESS_CACHING_OMP)"
    echo "  --f1r2_max_depth  No description (Default: 2500)"
    echo "  --min_base_q      No description (Default: 20)"
    echo "  --extra_args      No description (Default: )"
    echo "  -h, --help       Show this help message"
    exit 1
}

# --- [Default Variable Declarations] ---
# YAML의 기본값들이 이곳에 Key="Value" 형태로 정의됩니다.
SeqID=""
BamDir=""
vcfDir=""
qcResDir=""
ReferenceFasta=""
TargetInterval=""
VcfGnomad=""
VcfPon=""
InputSuffix="analysisReady"
singularity_bin="singularity"
gatk4_sif="/storage/images/gatk-4.4.0.0.sif"
bind="/storage,/data"
xmx_mb="32768"
Threads="14"
pairHMM="AVX_LOGLESS_CACHING_OMP"
f1r2_max_depth="2500"
min_base_q="20"
extra_args=""

# --- [Argument Parsing] ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --SeqID) SeqID="$2"; shift 2 ;;
        --BamDir) BamDir="$2"; shift 2 ;;
        --vcfDir) vcfDir="$2"; shift 2 ;;
        --qcResDir) qcResDir="$2"; shift 2 ;;
        --ReferenceFasta) ReferenceFasta="$2"; shift 2 ;;
        --TargetInterval) TargetInterval="$2"; shift 2 ;;
        --VcfGnomad) VcfGnomad="$2"; shift 2 ;;
        --VcfPon) VcfPon="$2"; shift 2 ;;
        --InputSuffix) InputSuffix="$2"; shift 2 ;;
        --singularity_bin) singularity_bin="$2"; shift 2 ;;
        --gatk4_sif) gatk4_sif="$2"; shift 2 ;;
        --bind) bind="$2"; shift 2 ;;
        --xmx_mb) xmx_mb="$2"; shift 2 ;;
        --Threads) Threads="$2"; shift 2 ;;
        --pairHMM) pairHMM="$2"; shift 2 ;;
        --f1r2_max_depth) f1r2_max_depth="$2"; shift 2 ;;
        --min_base_q) min_base_q="$2"; shift 2 ;;
        --extra_args) extra_args="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# --- [Validation] ---
# 필수 입력값(required: true)이 비어있는지 체크합니다.
if [[ -z "${SeqID:-}" ]]; then echo "Error: --SeqID is required"; usage; fi
if [[ -z "${BamDir:-}" ]]; then echo "Error: --BamDir is required"; usage; fi
if [[ -z "${vcfDir:-}" ]]; then echo "Error: --vcfDir is required"; usage; fi
if [[ -z "${qcResDir:-}" ]]; then echo "Error: --qcResDir is required"; usage; fi
if [[ -z "${ReferenceFasta:-}" ]]; then echo "Error: --ReferenceFasta is required"; usage; fi
if [[ -z "${TargetInterval:-}" ]]; then echo "Error: --TargetInterval is required"; usage; fi
if [[ -z "${VcfGnomad:-}" ]]; then echo "Error: --VcfGnomad is required"; usage; fi
if [[ -z "${VcfPon:-}" ]]; then echo "Error: --VcfPon is required"; usage; fi

# --- [Output Paths] ---
raw_vcf="${vcfDir}/${SeqID}.mutect2.vcf"
f1r2_tar_gz="${qcResDir}/${SeqID}.f1r2.tar.gz"
mutect_stats="${vcfDir}/${SeqID}.mutect2.vcf.stats"

# --- [Command Execution] ---
# [Key]가 $Key 형태로 치환된 최종 커맨드입니다.
cmd="${singularity_bin} exec -B ${bind} ${gatk4_sif} gatk Mutect2 --java-options '-XX:+UseParallelGC -XX:ParallelGCThreads=${Threads} -Xmx${xmx_mb}m' -pairHMM ${pairHMM} --reference ${ReferenceFasta} --intervals ${TargetInterval} --input ${BamDir}/${SeqID}.${InputSuffix}.bam --germline-resource ${VcfGnomad} --panel-of-normals ${VcfPon} --f1r2-tar-gz ${f1r2_tar_gz} --f1r2-max-depth ${f1r2_max_depth} --min-base-quality-score ${min_base_q} --base-quality-score-threshold ${min_base_q} --output ${raw_vcf} ${extra_args}"

echo -e "\\n[RUNNING]\\n$cmd\\n"

# 자동 디렉토리 생성
mkdir -p "$(dirname "${BamDir}")" 2>/dev/null || mkdir -p "${BamDir}"
mkdir -p "$(dirname "${vcfDir}")" 2>/dev/null || mkdir -p "${vcfDir}"
mkdir -p "$(dirname "${qcResDir}")" 2>/dev/null || mkdir -p "${qcResDir}"

eval "$cmd"