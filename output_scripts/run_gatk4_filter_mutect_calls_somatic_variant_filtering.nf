// [METADATA]
// TOOL_NAME = gatk4_filter_mutect_calls
// THREADS = 1

process GATK4_FILTER_MUTECT_CALLS {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(vcfDir), val(qcResDir), val(ReferenceFasta), val(TargetInterval), val(ContaminationTable)

    output:
    path "*", emit: filtered_vcf

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def vcfDir = params.vcfDir ?: ""
    def qcResDir = params.qcResDir ?: ""
    def ReferenceFasta = params.ReferenceFasta ?: ""
    def TargetInterval = params.TargetInterval ?: ""
    def ContaminationTable = params.ContaminationTable ?: ""
    def Suffix = params.Suffix ?: ""
    def singularity_bin = params.singularity_bin ?: "singularity"
    def gatk4_sif = params.gatk4_sif ?: "/storage/images/gatk-4.4.0.0.sif"
    def bind = params.bind ?: "/storage,/data"
    def xmx_mb = params.xmx_mb ?: "16384"
    def Threads = params.Threads ?: "14"
    def extra_args = params.extra_args ?: ""
    """
    ${singularity_bin} exec -B ${bind} ${gatk4_sif} gatk FilterMutectCalls --java-options '-XX:ParallelGCThreads=${Threads} -Xmx${xmx_mb}m' -V ${vcfDir}/${SeqID}.mutect2.${Suffix}vcf -L ${TargetInterval} --reference ${ReferenceFasta} --contamination-table ${ContaminationTable} --ob-priors ${qcResDir}/${SeqID}.${Suffix}read-orientation-model.tar.gz -O ${filtered_vcf} ${extra_args}
    """
}