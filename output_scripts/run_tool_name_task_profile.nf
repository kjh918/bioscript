// [METADATA]
// TOOL_NAME = TOOL_NAME
// THREADS = 8

process TOOL_NAME {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), val(InputDir), val(OutDir)

    output:
    path "*", emit: output_file

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def InputDir = params.InputDir ?: ""
    def OutDir = params.OutDir ?: ""
    def tool_bin = params.tool_bin ?: "TOOL_COMMAND"
    def InputSuffix = params.InputSuffix ?: ".UseAnalysis.bam"
    def OutputSuffix = params.OutputSuffix ?: ".txt"
    def Threads = params.Threads ?: "8"
    """
    ${tool_bin} INPUT_OPTIONS OUTPUT_OPTIONS
    """
}