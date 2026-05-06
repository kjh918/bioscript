// [METADATA]
// TOOL_NAME = rsem
// THREADS = 1

process RSEM {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(TranscriptomeBam), val(InputSuffix), val(RsemIndex), val(OutputDir)

    output:

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def TranscriptomeBam = params.TranscriptomeBam ?: ""
    def InputSuffix = params.InputSuffix ?: ".Aligned.toTranscriptome.out.bam"
    def RsemIndex = params.RsemIndex ?: ""
    def OutputDir = params.OutputDir ?: ""
    def Threads = params.Threads ?: "12"
    def PairedFlag = params.PairedFlag ?: "--paired-end"
    def NoBamOutput = params.NoBamOutput ?: "--no-bam-output"
    def EstimateRspd = params.EstimateRspd ?: "--estimate-rspd"
    def StrandOption = params.StrandOption ?: "none"
    def ExtraArgs = params.ExtraArgs ?: "None"
    def singularity_bin = params.singularity_bin ?: "singularity"
    def rsem_sif = params.rsem_sif ?: "/storage/images/rsem-1.3.3.sif"
    def bind = params.bind ?: "/storage,/data"
    """
    ${singularity_bin} exec -B ${bind} ${rsem_sif} rsem-calculate-expression ${PairedFlag} ${NoBamOutput} ${EstimateRspd} ${ExtraArgs} --num-threads ${Threads} --strandedness ${StrandOption} --alignments ${TranscriptomeBam}/${SeqID}${InputSuffix} ${RsemIndex} ${OutputDir}/${SeqID}
    """
}