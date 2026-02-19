// [METADATA]
// TOOL_NAME = bedtools
// THREADS = 1

process BEDTOOLS {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(BamDir), path(BamToBedDir)

    output:
    path "*", emit: bed_gz
    path "*", emit: final_bed

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def BamDir = params.BamDir ?: ""
    def BamToBedDir = params.BamToBedDir ?: ""
    def InputSuffix = params.InputSuffix ?: "analysisReady"
    def bedtools_bin = params.bedtools_bin ?: "bedtools"
    def bgzip_bin = params.bgzip_bin ?: "bgzip"
    def ln_bin = params.ln_bin ?: "ln"
    def cp_bin = params.cp_bin ?: "cp"
    def Threads = params.Threads ?: "8"
    """
    ${bedtools_bin} bamtobed -i ${BamDir}/${SeqID}.${InputSuffix}.bam > ${BamDir}/${SeqID}.${InputSuffix}.bed && ${bgzip_bin} -f -@ ${Threads} ${BamDir}/${SeqID}.${InputSuffix}.bed && ${ln_bin} -Tsf ${bed_gz} ${BamDir}/${SeqID}.bed.gz && ${cp_bin} ${bed_gz} ${final_bed}
    """
}