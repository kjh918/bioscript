# [METADATA]
# TOOL_NAME = vcf2maf
# VERSION = 1.6.21
# THREADS = 1

rule vcf2maf:
    input:
        SeqID = ""
        vcfDir = ""
        VcfTag = ""
        NormalID = "NORMAL"
        genomeFasta = ""
    output:
        output_maf = "[vcfDir]/[SeqID].[VcfTag].maf"
    params:
        retain_info = "DP,GERMQ,MBQ,MFRL,MMQ,MPOS,POPAF,PON,ROQ,TLOD,FILTER,pArtifact,artiStatus"
        retain_ann = "HGVSg,am_class,am_pathogenicity,CADD_PHRED,CADD_RAW,ClinVar_CLNSIG,ClinVar_CLNDN,ClinVar_ORIGIN,gnomADe_AF,gnomADe_AFR_AF,gnomADe_AMR_AF,gnomADe_ASJ_AF,gnomADe_EAS_AF,gnomADe_FIN_AF,gnomADe_NFE_AF,gnomADe_OTH_AF,gnomADe_SAS_AF,MAX_AF,MAX_AF_POPS,TRANSCRIPTION_FACTORS,COSMIC,COSMIC_HGVSC,COSMIC_LEGACY_ID,COSMIC_TIER"
        singularity_bin = "singularity"
        vcf2maf_sif = "/storage/images/vcf2maf.sif"
        vcf2maf_bin = "vcf2maf.pl"
        bind = "/storage,/data"
    threads: 1
    shell:
        """
        {params.singularity_bin} exec -B {params.bind} {params.vcf2maf_sif} {params.vcf2maf_bin} --inhibit-vep --input-vcf {input.vcfDir}/{input.SeqID}.{input.VcfTag}.vcf --output-maf {output.output_maf} --tumor-id {input.SeqID} --normal-id {input.NormalID} --ref-fasta {input.genomeFasta} --retain-info {params.retain_info} --retain-ann {params.retain_ann}
        """