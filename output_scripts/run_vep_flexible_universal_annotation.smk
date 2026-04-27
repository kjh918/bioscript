# [METADATA]
# TOOL_NAME = vep_flexible
# VERSION = 110
# THREADS = 1

rule vep_flexible:
    input:
        SeqID = ""
        vcfDir = ""
        VcfTag = ""
        vepCacheDir = ""
        vepPluginDir = ""
        genomeFasta = ""
        clinvarData = ""
        cosmicData = ""
        alphaMissense = ""
        caddSnp = ""
        caddIndel = ""
    output:
        vep_vcf = "[vcfDir]/[SeqID].[VcfTag].vep.refseq.vcf"
        vep_stats = "[vcfDir]/[SeqID].[VcfTag].vep.refseq.vcf_summary.txt"
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
        {params.singularity_bin} exec -B {params.bind} {params.vep_sif} vep --force_overwrite --offline --cache --dir_cache {input.vepCacheDir} --dir_plugins {input.vepPluginDir} --fasta {input.genomeFasta} --species {params.species} --assembly {params.assembly} --input_file {input.vcfDir}/{input.SeqID}.{input.VcfTag}.vcf --output_file {output.vep_vcf} --vcf --stats_text --refseq --show_ref_allele --uploaded_allele --use_transcript_ref --variant_class --sift b --polyphen b --gene_phenotype --numbers --hgvs --hgvsg --symbol --canonical --biotype --regulatory --mirna --check_existing --max_af --af_1kg --af_gnomade --exclude_predicted --pick --flag_pick --fork {threads} --buffer_size {params.buffer_size} {params.clinvar_args} {params.cosmic_args} --plugin AlphaMissense,file={input.alphaMissense} --plugin CADD,snv={input.caddSnp},indels={input.caddIndel}
        """