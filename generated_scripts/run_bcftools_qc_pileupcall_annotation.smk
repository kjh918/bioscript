# [METADATA]
# TOOL_NAME = bcftools
# VERSION = 1.23
# THREADS = 1

rule bcftools:
    input:
        SeqID = ""
        Chromosome = ""
        BamDir = ""
        ResultDir = ""
        ReferenceFasta = ""
        SitesVcfGz = ""
        PopAfAnnotVcf = ""
        PopAfHeaderHdr = ""
    output:
        raw_vcf = "[ResultDir]/[SeqID].[Chromosome].[InputSuffix].raw.vcf.gz"
        ann_vcf = "[ResultDir]/[SeqID].[Chromosome].[InputSuffix].ann.vcf.gz"
    params:
        InputSuffix = "analysisReady"
        bcftools_bin = "/storage/home/jhkim/Apps/bcftools/bcftools"
        Threads = "4"
        MinBQ = "20"
        MinMQ = "30"
        AnnotationQuery = "CHROM,POS,REF,ALT,INFO/KOVA_AF,INFO/KOVA_AN,INFO/GNOMAD_AF,INFO/GNOMAD_AN"
    threads: 1
    shell:
        """
        {params.bcftools_bin} mpileup -f {input.ReferenceFasta} -T {input.SitesVcfGz} -r {input.Chromosome} -q {params.MinMQ} -Q {params.MinBQ} -a FORMAT/AD,FORMAT/DP -Ou --threads {threads} {input.BamDir}/{input.SeqID}.{params.InputSuffix}.bam | {params.bcftools_bin} call -Am -Oz -o {output.raw_vcf} --threads {threads} && {params.bcftools_bin} index -f {output.raw_vcf} --threads {threads} && {params.bcftools_bin} annotate -a {input.PopAfAnnotVcf} -c {params.AnnotationQuery} -h {input.PopAfHeaderHdr} -Oz -o {output.ann_vcf} {output.raw_vcf} --threads {threads} && {params.bcftools_bin} index -f {output.ann_vcf} --threads {threads}
        """