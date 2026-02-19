// [METADATA]
// TOOL_NAME = multiqc
// THREADS = 1

process MULTIQC {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), val(qcResDir)

    output:
    path "*", emit: report_html

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def qcResDir = params.qcResDir ?: ""
    def singularity_bin = params.singularity_bin ?: "singularity"
    def multiqc_bin = params.multiqc_bin ?: "multiqc"
    def sif = params.sif ?: "/storage/images/multiqc-1.16.sif"
    def bind = params.bind ?: "/storage,/data"
    def mqc_config = params.mqc_config ?: "/storage/home/kangsm/runScripts/NGS_config.MultiQC_Custom.yaml"
    def mqc_filename = params.mqc_filename ?: "[SeqID].QC.Results"
    def mqc_args = params.mqc_args ?: "--force --data-dir"
    """
    ${singularity_bin} exec -B ${bind} ${sif} ${multiqc_bin} ${mqc_args} --filename ${mqc_filename} --outdir ${qcResDir} --config ${mqc_config} ${qcResDir}
    """
}