// [METADATA]
// TOOL_NAME = fastp
// THREADS = 1

process FASTP {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(RawFastqDir), path(TrimFastqDir), val(qcResDir)

    output:
    path "*", emit: out_read1
    path "*", emit: out_read2
    path "*", emit: json
    path "*", emit: html

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def RawFastqDir = params.RawFastqDir ?: ""
    def TrimFastqDir = params.TrimFastqDir ?: ""
    def qcResDir = params.qcResDir ?: ""
    def singularity_bin = params.singularity_bin ?: "singularity"
    def fastp_bin = params.fastp_bin ?: "fastp"
    def sif = params.sif ?: "/storage/images/fastp-0.23.4.sif"
    def bind = params.bind ?: "/storage,/data"
    def Threads = params.Threads ?: "8"
    def length_required = params.length_required ?: "100"
    def average_qual = params.average_qual ?: "10"
    def qualified_quality_phred = params.qualified_quality_phred ?: "15"
    """
    ${singularity_bin} exec -B ${bind} ${sif} ${fastp_bin} --thread ${Threads} --in1 ${RawFastqDir}/${SeqID}_R1.fastq.gz --in2 ${RawFastqDir}/${SeqID}_R2.fastq.gz --out1 ${out_read1} --out2 ${out_read2} --json ${json} --html ${html} --trim_poly_g --detect_adapter_for_pe --length_required ${length_required} --average_qual ${average_qual} --qualified_quality_phred ${qualified_quality_phred}
    """
}