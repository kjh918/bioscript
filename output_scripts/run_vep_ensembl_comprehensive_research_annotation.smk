# [METADATA]
# TOOL_NAME = vep_ensembl
# VERSION = 110
# THREADS = 1

rule vep_ensembl:
    input:
        SeqID = ""
        vcfDir = ""
        VcfTag = ""
        vepCacheDir = ""
        genomeFasta = ""
        vepBam = ""
        clinvarData = ""
        cosmicData = ""
        alphaMissense = ""
        caddSnp = ""
        caddIndel = ""
    output:
        vep_vcf = "[vcfDir]/[SeqID].[VcfTag].vep.ensembl.vcf"
        vep_stats = "[vcfDir]/[SeqID].[VcfTag].vep.ensembl.vcf_summary.txt"
    params:
        assembly = "GRCh38"
        species = "homo_sapiens"
        Threads = "8"
        buffer_size = "50000"
        cosmic_args = ""
        clinvar_args = "--custom file=[clinvarData],short_name=ClinVar,format=vcf,fields=CLNSIG%CLNDN%ORIGIN"
        singularity_bin = "singularity"
        vep_sif = "/storage/images/vep-110.sif"
        bind = "/storage,/data"
    threads: 1
    shell:
        """
        {params.singularity_bin} exec -B {params.bind} {params.vep_sif} vep --force_overwrite --offline --cache --dir_cache {input.vepCacheDir} --fasta {input.genomeFasta} --bam {input.vepBam} --species {params.species} --assembly {params.assembly} --input_file {input.vcfDir}/{input.SeqID}.{input.VcfTag}.vcf --output_file {output.vep_vcf} --vcf --stats_text --hgvs --hgvsg --canonical --exclude_predicted --ccds --uniprot --domains --fork {threads} --buffer_size {params.buffer_size} {params.clinvar_args} {params.cosmic_args} --plugin AlphaMissense,file={input.alphaMissense} --plugin CADD,snv={input.caddSnp},indels={input.caddIndel}
        """