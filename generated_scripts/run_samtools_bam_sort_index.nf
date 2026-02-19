// [METADATA]
// TOOL_NAME = samtools
// THREADS = 1

process SAMTOOLS {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(BamDir)

    output:
    path "*", emit: sorted_bam
    path "*", emit: sorted_bai

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def BamDir = params.BamDir ?: ""
    def OutputSuffix = params.OutputSuffix ?: "sorted"
    def samtools_bin = params.samtools_bin ?: "samtools"
    def Threads = params.Threads ?: "8"
    """
    ${samtools_bin} sort -@ ${Threads} -o ${sorted_bam} ${BamDir}/${SeqID}.bam && ${samtools_bin} index -b -@ ${Threads} ${sorted_bam}
    """
}