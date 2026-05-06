// [METADATA]
// TOOL_NAME = salmon
// THREADS = 1

process SALMON {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(BamDir), val(OutputDir), val(SalmonIndex)

    output:
    path "*", emit: quant_sf
    path "*", emit: quant_log
    path "*", emit: lib_format

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def BamDir = params.BamDir ?: ""
    def OutputDir = params.OutputDir ?: ""
    def SalmonIndex = params.SalmonIndex ?: ""
    def Threads = params.Threads ?: "12"
    def libType = params.libType ?: "A"
    def InputSuffix = params.InputSuffix ?: ".Aligned.toTranscriptome.out.bam"
    def extra_args = params.extra_args ?: ""
    def samlon_bin = params.samlon_bin ?: "/storage/apps/salmon-1.10.0/bin/salmon"
    """
    ${samlon_bin} quant -t ${SalmonIndex} --threads ${Threads} --libType ${libType} -a ${BamDir}/${SeqID}${InputSuffix} -o ${OutputDir} ${extra_args}
    """
}