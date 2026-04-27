// [METADATA]
// TOOL_NAME = vcf2maf
// THREADS = 1

process VCF2MAF {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(vcfDir), path(VcfTag), val(NormalID), val(genomeFasta)

    output:
    path "*", emit: output_maf

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def vcfDir = params.vcfDir ?: ""
    def VcfTag = params.VcfTag ?: ""
    def NormalID = params.NormalID ?: "NORMAL"
    def genomeFasta = params.genomeFasta ?: ""
    def retain_info = params.retain_info ?: "DP,GERMQ,MBQ,MFRL,MMQ,MPOS,POPAF,PON,ROQ,TLOD,FILTER,pArtifact,artiStatus"
    def retain_ann = params.retain_ann ?: "HGVSg,am_class,am_pathogenicity,CADD_PHRED,CADD_RAW,ClinVar_CLNSIG,ClinVar_CLNDN,ClinVar_ORIGIN,gnomADe_AF,gnomADe_AFR_AF,gnomADe_AMR_AF,gnomADe_ASJ_AF,gnomADe_EAS_AF,gnomADe_FIN_AF,gnomADe_NFE_AF,gnomADe_OTH_AF,gnomADe_SAS_AF,MAX_AF,MAX_AF_POPS,TRANSCRIPTION_FACTORS,COSMIC,COSMIC_HGVSC,COSMIC_LEGACY_ID,COSMIC_TIER"
    def singularity_bin = params.singularity_bin ?: "singularity"
    def vcf2maf_sif = params.vcf2maf_sif ?: "/storage/images/vcf2maf.sif"
    def vcf2maf_bin = params.vcf2maf_bin ?: "vcf2maf.pl"
    def bind = params.bind ?: "/storage,/data"
    """
    ${singularity_bin} exec -B ${bind} ${vcf2maf_sif} ${vcf2maf_bin} --inhibit-vep --input-vcf ${vcfDir}/${SeqID}.${VcfTag}.vcf --output-maf ${output_maf} --tumor-id ${SeqID} --normal-id ${NormalID} --ref-fasta ${genomeFasta} --retain-info ${retain_info} --retain-ann ${retain_ann}
    """
}