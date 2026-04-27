// [METADATA]
// TOOL_NAME = gatk4_mutect2
// THREADS = 1

process GATK4_MUTECT2 {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(BamDir), path(vcfDir), val(qcResDir), val(ReferenceFasta), val(TargetInterval), path(VcfGnomad), path(VcfPon)

    output:
    path "*", emit: raw_vcf
    path "*", emit: f1r2_tar_gz
    path "*", emit: mutect_stats

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def BamDir = params.BamDir ?: ""
    def vcfDir = params.vcfDir ?: ""
    def qcResDir = params.qcResDir ?: ""
    def ReferenceFasta = params.ReferenceFasta ?: ""
    def TargetInterval = params.TargetInterval ?: ""
    def VcfGnomad = params.VcfGnomad ?: ""
    def VcfPon = params.VcfPon ?: ""
    def InputSuffix = params.InputSuffix ?: "analysisReady"
    def singularity_bin = params.singularity_bin ?: "singularity"
    def gatk4_sif = params.gatk4_sif ?: "/storage/images/gatk-4.4.0.0.sif"
    def bind = params.bind ?: "/storage,/data"
    def xmx_mb = params.xmx_mb ?: "32768"
    def Threads = params.Threads ?: "14"
    def pairHMM = params.pairHMM ?: "AVX_LOGLESS_CACHING_OMP"
    def f1r2_max_depth = params.f1r2_max_depth ?: "2500"
    def min_base_q = params.min_base_q ?: "20"
    def extra_args = params.extra_args ?: ""
    """
    ${singularity_bin} exec -B ${bind} ${gatk4_sif} gatk Mutect2 --java-options '-XX:+UseParallelGC -XX:ParallelGCThreads=${Threads} -Xmx${xmx_mb}m' -pairHMM ${pairHMM} --reference ${ReferenceFasta} --intervals ${TargetInterval} --input ${BamDir}/${SeqID}.${InputSuffix}.bam --germline-resource ${VcfGnomad} --panel-of-normals ${VcfPon} --f1r2-tar-gz ${f1r2_tar_gz} --f1r2-max-depth ${f1r2_max_depth} --min-base-quality-score ${min_base_q} --base-quality-score-threshold ${min_base_q} --output ${raw_vcf} ${extra_args}
    """
}