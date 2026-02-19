// [METADATA]
// TOOL_NAME = fastqc
// THREADS = 8

process FASTQC {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(RawFastqDir), val(qcResDir)

    output:
    path "*", emit: r1_html
    path "*", emit: r2_html
    path "*", emit: r1_zip
    path "*", emit: r2_zip
    path "*", emit: r1_dir
    path "*", emit: r2_dir

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def RawFastqDir = params.RawFastqDir ?: ""
    def qcResDir = params.qcResDir ?: ""
    def singularity_bin = params.singularity_bin ?: "singularity"
    def fastqc_bin = params.fastqc_bin ?: "fastqc"
    def sif = params.sif ?: "/storage/images/fastqc-0.12.1.sif"
    def bind = params.bind ?: "/storage,/data"
    def threads = params.threads ?: "8"
    def fastqc_args = params.fastqc_args ?: "--extract"
    """
    ${singularity_bin} exec -B ${bind} ${sif} ${fastqc_bin} ${fastqc_args} --threads ${threads} --outdir ${qcResDir} ${RawFastqDir}/${SeqID}_R1.fastq.gz ${RawFastqDir}/${SeqID}_R2.fastq.gz
    """
}