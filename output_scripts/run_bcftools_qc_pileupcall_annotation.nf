// [METADATA]
// TOOL_NAME = bcftools
// THREADS = 1

process BCFTOOLS {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), val(Chromosome), path(BamDir), val(ResultDir), val(ReferenceFasta), path(SitesVcfGz), path(PopAfAnnotVcf), val(PopAfHeaderHdr)

    output:
    path "*", emit: raw_vcf
    path "*", emit: ann_vcf

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def Chromosome = params.Chromosome ?: ""
    def BamDir = params.BamDir ?: ""
    def ResultDir = params.ResultDir ?: ""
    def ReferenceFasta = params.ReferenceFasta ?: ""
    def SitesVcfGz = params.SitesVcfGz ?: ""
    def PopAfAnnotVcf = params.PopAfAnnotVcf ?: ""
    def PopAfHeaderHdr = params.PopAfHeaderHdr ?: ""
    def InputSuffix = params.InputSuffix ?: "analysisReady"
    def bcftools_bin = params.bcftools_bin ?: "/storage/home/jhkim/Apps/bcftools/bcftools"
    def Threads = params.Threads ?: "4"
    def MinBQ = params.MinBQ ?: "20"
    def MinMQ = params.MinMQ ?: "30"
    def AnnotationQuery = params.AnnotationQuery ?: "CHROM,POS,REF,ALT,INFO/KOVA_AF,INFO/KOVA_AN,INFO/GNOMAD_AF,INFO/GNOMAD_AN"
    """
    ${bcftools_bin} mpileup -f ${ReferenceFasta} -T ${SitesVcfGz} -r ${Chromosome} -q ${MinMQ} -Q ${MinBQ} -a FORMAT/AD,FORMAT/DP -Ou --threads ${Threads} ${BamDir}/${SeqID}.${InputSuffix}.bam | ${bcftools_bin} call -Am -Oz -o ${raw_vcf} --threads ${Threads} && ${bcftools_bin} index -f ${raw_vcf} --threads ${Threads} && ${bcftools_bin} annotate -a ${PopAfAnnotVcf} -c ${AnnotationQuery} -h ${PopAfHeaderHdr} -Oz -o ${ann_vcf} ${raw_vcf} --threads ${Threads} && ${bcftools_bin} index -f ${ann_vcf} --threads ${Threads}
    """
}