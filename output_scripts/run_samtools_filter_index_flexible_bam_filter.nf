// [METADATA]
// TOOL_NAME = samtools_filter_index
// THREADS = 1

process SAMTOOLS_FILTER_INDEX {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(BamDir)

    output:
    path "*", emit: filtered_bam
    path "*", emit: filtered_bai
    path "*", emit: analysis_ready_bam
    path "*", emit: analysis_ready_bai

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def BamDir = params.BamDir ?: ""
    def InputSuffix = params.InputSuffix ?: "recal"
    def OutputSuffix = params.OutputSuffix ?: "filtered"
    def samtools_bin = params.samtools_bin ?: "samtools"
    def ln_bin = params.ln_bin ?: "ln"
    def Threads = params.Threads ?: "8"
    def min_mapq = params.min_mapq ?: "20"
    def include_flag = params.include_flag ?: "0x2"
    def exclude_flag = params.exclude_flag ?: "0x104"
    def expr_nm = params.expr_nm ?: "[NM] < 12"
    """
    ${samtools_bin} view -b -h -q ${min_mapq} -f ${include_flag} -F ${exclude_flag} -e '${expr_nm}' --threads ${Threads} ${BamDir}/${SeqID}.${InputSuffix}.bam > ${filtered_bam} && ${samtools_bin} index --threads ${Threads} ${filtered_bam} && ${ln_bin} -Tsf ${filtered_bam} ${analysis_ready_bam} && ${ln_bin} -Tsf ${filtered_bai} ${analysis_ready_bai}
    """
}