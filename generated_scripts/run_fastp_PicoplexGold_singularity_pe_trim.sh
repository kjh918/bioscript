#!/usr/bin/env bash
set -euo pipefail

# --- Defaults ---
RawFastqDir=''
SeqID=''
Threads=8
TrimFastqDir=''
adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
average_qual=10
bind=/storage,/data
html='[qcResDir]/[SeqID].fastp.html'
json='[qcResDir]/[SeqID].fastp.json'
length_required=100
out_read1='[TrimFastqDir]/[SeqID].trimmed_R1.fastq.gz'
out_read2='[TrimFastqDir]/[SeqID].trimmed_R2.fastq.gz'
qcResDir=''
qualified_quality_phred=15
sif=/storage/images/fastp-0.23.4.sif
trim_front1=14
trim_front2=14
CWD="."

# --- CLI Argument Parsing ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    --RawFastqDir) RawFastqDir="$2"; shift 2;;
    --SeqID) SeqID="$2"; shift 2;;
    --Threads) Threads="$2"; shift 2;;
    --TrimFastqDir) TrimFastqDir="$2"; shift 2;;
    --adapter_sequence) adapter_sequence="$2"; shift 2;;
    --adapter_sequence_r2) adapter_sequence_r2="$2"; shift 2;;
    --average_qual) average_qual="$2"; shift 2;;
    --bind) bind="$2"; shift 2;;
    --html) html="$2"; shift 2;;
    --json) json="$2"; shift 2;;
    --length_required) length_required="$2"; shift 2;;
    --out_read1) out_read1="$2"; shift 2;;
    --out_read2) out_read2="$2"; shift 2;;
    --qcResDir) qcResDir="$2"; shift 2;;
    --qualified_quality_phred) qualified_quality_phred="$2"; shift 2;;
    --sif) sif="$2"; shift 2;;
    --trim_front1) trim_front1="$2"; shift 2;;
    --trim_front2) trim_front2="$2"; shift 2;;
    --cwd) CWD="$2"; shift 2;;
    -h|--help) echo "Usage: $0 [options]"; exit 0;;
    *) echo "Unknown argument: $1" >&2; exit 1;;
  esac
done

# --- Template Engine ---
render() {
  local s="$1"
  # 중첩된 변수 치환을 위해 3회 반복 (예: [A] -> [B] -> value)
  for i in {1..3}; do
    s="${s//\\[RawFastqDir\\]/${RawFastqDir}}"
    s="${s//\\[SeqID\\]/${SeqID}}"
    s="${s//\\[Threads\\]/${Threads}}"
    s="${s//\\[TrimFastqDir\\]/${TrimFastqDir}}"
    s="${s//\\[adapter_sequence\\]/${adapter_sequence}}"
    s="${s//\\[adapter_sequence_r2\\]/${adapter_sequence_r2}}"
    s="${s//\\[average_qual\\]/${average_qual}}"
    s="${s//\\[bind\\]/${bind}}"
    s="${s//\\[html\\]/${html}}"
    s="${s//\\[json\\]/${json}}"
    s="${s//\\[length_required\\]/${length_required}}"
    s="${s//\\[out_read1\\]/${out_read1}}"
    s="${s//\\[out_read2\\]/${out_read2}}"
    s="${s//\\[qcResDir\\]/${qcResDir}}"
    s="${s//\\[qualified_quality_phred\\]/${qualified_quality_phred}}"
    s="${s//\\[sif\\]/${sif}}"
    s="${s//\\[trim_front1\\]/${trim_front1}}"
    s="${s//\\[trim_front2\\]/${trim_front2}}"
  done
  # 앞뒤 공백 제거 및 연속 공백 축소
  echo "$s" | xargs
}

# --- Finalize Outputs & Command ---
out_read1=$(render "${out_read1}")
out_read2=$(render "${out_read2}")
json=$(render "${json}")
html=$(render "${html}")
CMD_LINE='singularity exec -B [bind] [sif] fastp --thread [Threads] --in1 [RawFastqDir]/[SeqID]_R1.fastq.gz --in2 [RawFastqDir]/[SeqID]_R2.fastq.gz --out1 [out_read1] --out2 [out_read2] --json [json] --html [html] --trim_poly_g --detect_adapter_for_pe --adapter_sequence [adapter_sequence] --adapter_sequence_r2 [adapter_sequence_r2] --length_required [length_required] --average_qual [average_qual] --qualified_quality_phred [qualified_quality_phred] --trim_front1 [trim_front1] --trim_front2 [trim_front2]'
CMD=$(render "$CMD_LINE")

echo -e "\n[RUNNING CMD]\n$CMD\n"

# --- Execution ---
mkdir -p "$CWD"
cd "$CWD" && eval "$CMD"
