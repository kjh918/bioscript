// [METADATA]
// TOOL_NAME = samtools_flagstats
// THREADS = 1

process SAMTOOLS_FLAGSTATS {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(BamDir)

    output:
    path "*", emit: flag_stats_txt

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def BamDir = params.BamDir ?: ""
    def InputSuffix = params.InputSuffix ?: ".sorted"
    def OutputSuffix = params.OutputSuffix ?: "filtered"
    def samtools_bin = params.samtools_bin ?: "samtools"
    """
    ${samtools_bin} flagstats ${BamDir}/[SeqId]${InputSuffix}.bam > ${flag_stats_txt}
    """
}