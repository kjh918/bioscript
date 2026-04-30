// [METADATA]
// TOOL_NAME = salmon_quant
// THREADS = 1

process SALMON_QUANT {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(FastqDir), val(OutputDir), val(SalmonIndex)

    output:
    path "*", emit: quant_sf
    path "*", emit: quant_log
    path "*", emit: lib_format

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def FastqDir = params.FastqDir ?: ""
    def OutputDir = params.OutputDir ?: ""
    def SalmonIndex = params.SalmonIndex ?: ""
    def Threads = params.Threads ?: "12"
    def libType = params.libType ?: "A"
    def InputSuffix = params.InputSuffix ?: "fastq.gz"
    def validateMappings = params.validateMappings ?: "--validateMappings"
    def extra_args = params.extra_args ?: ""
    """
    salmon quant --threads ${Threads} --index ${SalmonIndex} --libType ${libType} -1 ${FastqDir}/${SeqID}_R1.${InputSuffix} -2 ${FastqDir}/${SeqID}_R2.${InputSuffix} ${validateMappings} -o ${OutputDir}/${SeqID}_quant ${extra_args}
    """
}