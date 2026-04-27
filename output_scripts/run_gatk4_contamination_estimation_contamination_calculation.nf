// [METADATA]
// TOOL_NAME = gatk4_contamination_estimation
// THREADS = 1

process GATK4_CONTAMINATION_ESTIMATION {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(BamDir), val(qcResDir), val(ReferenceFasta), val(TargetInterval), path(VcfGnomad)

    output:
    path "*", emit: pileup_table
    path "*", emit: contamination_table

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def BamDir = params.BamDir ?: ""
    def qcResDir = params.qcResDir ?: ""
    def ReferenceFasta = params.ReferenceFasta ?: ""
    def TargetInterval = params.TargetInterval ?: ""
    def VcfGnomad = params.VcfGnomad ?: ""
    def InputSuffix = params.InputSuffix ?: "analysisReady"
    def singularity_bin = params.singularity_bin ?: "singularity"
    def gatk4_sif = params.gatk4_sif ?: "/storage/images/gatk-4.4.0.0.sif"
    def bind = params.bind ?: "/storage,/data"
    def Threads = params.Threads ?: "14"
    def pileup_xmx = params.pileup_xmx ?: "8192"
    def contam_xmx = params.contam_xmx ?: "16384"
    """
    ${singularity_bin} exec -B ${bind} ${gatk4_sif} gatk GetPileupSummaries --java-options '-XX:ParallelGCThreads=${Threads} -Xmx${pileup_xmx}m' -I ${BamDir}/${SeqID}.${InputSuffix}.bam -V ${VcfGnomad} -L ${TargetInterval} -R ${ReferenceFasta} -O ${pileup_table} && ${singularity_bin} exec -B ${bind} ${gatk4_sif} gatk CalculateContamination --java-options '-XX:ParallelGCThreads=${Threads} -Xmx${contam_xmx}m' -I ${pileup_table} -O ${contamination_table}
    """
}