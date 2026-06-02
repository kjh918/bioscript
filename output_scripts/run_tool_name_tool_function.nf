// [METADATA]
// TOOL_NAME = TOOL_NAME
// THREADS = 8

process TOOL_NAME {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), val(InputFileDir), val(qcResDir), val(OutputDir)

    output:
    path "*", emit: OutputSummary

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def InputFileDir = params.InputFileDir ?: ""
    def qcResDir = params.qcResDir ?: ""
    def OutputDir = params.OutputDir ?: ""
    def InputSuffix = params.InputSuffix ?: ""
    def OutputSuffix = params.OutputSuffix ?: ""
    def singularity_bin = params.singularity_bin ?: "singularity"
    def bind = params.bind ?: "/storage,/data"
    def xmx_mb = params.xmx_mb ?: "16384"
    def Threads = params.Threads ?: "8"
    def extra_args = params.extra_args ?: ""
    """
    ${singularity_bin} exec -B ${bind}  ${extra_args}
    """
}