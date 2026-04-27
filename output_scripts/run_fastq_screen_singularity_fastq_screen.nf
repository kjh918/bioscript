// [METADATA]
// TOOL_NAME = fastq_screen
// THREADS = 1

process FASTQ_SCREEN {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(RawFastqDir), val(qcResDir)

    output:
    path "*", emit: screen_txt
    path "*", emit: screen_png
    path "*", emit: screen_html

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def RawFastqDir = params.RawFastqDir ?: ""
    def qcResDir = params.qcResDir ?: ""
    def singularity_bin = params.singularity_bin ?: "singularity"
    def fastq_screen_bin = params.fastq_screen_bin ?: "fastq_screen"
    def sif = params.sif ?: "/storage/images/fastqScreen-0.15.3.sif"
    def bind = params.bind ?: "/storage,/data"
    def aligner = params.aligner ?: "bowtie2"
    def config_path = params.config_path ?: "/storage/home/kangsm/runScripts/NGS_config.FastqScreen.conf"
    def Threads = params.Threads ?: "15"
    def extra_args = params.extra_args ?: "--force"
    """
    ${singularity_bin} exec -B ${bind} ${sif} ${fastq_screen_bin} --aligner ${aligner} --conf ${config_path} --outdir ${qcResDir} --threads ${Threads} ${extra_args} ${RawFastqDir}/${SeqID}_R1.fastq.gz
    """
}