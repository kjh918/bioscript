// [METADATA]
// TOOL_NAME = sob_detector
// THREADS = 1

process SOB_DETECTOR {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(vcfDir), path(BamDir)

    output:
    path "*", emit: bias_filtered_vcf

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def vcfDir = params.vcfDir ?: ""
    def BamDir = params.BamDir ?: ""
    def InputVcfSuffix = params.InputVcfSuffix ?: "filtered"
    def InputBamSuffix = params.InputBamSuffix ?: "analysisReady"
    def singularity_bin = params.singularity_bin ?: "singularity"
    def sob_sif = params.sob_sif ?: "/storage/images/sobdetector.sif"
    def sob_bin = params.sob_bin ?: "SOBDetector"
    def bind = params.bind ?: "/storage,/data"
    def minBaseQ = params.minBaseQ ?: "20"
    def minMapQ = params.minMapQ ?: "20"
    def only_passed = params.only_passed ?: "false"
    def extra_args = params.extra_args ?: ""
    """
    ${singularity_bin} exec -B ${bind} ${sob_sif} ${sob_bin}  --input-type VCF --input-variants ${vcfDir}/${SeqID}.mutect2.${InputVcfSuffix}.vcf --input-bam ${BamDir}/${SeqID}.${InputBamSuffix}.bam --output-variants ${bias_filtered_vcf} --minBaseQuality ${minBaseQ} --minMappingQuality ${minMapQ} --only-passed ${only_passed} ${extra_args}
    """
}