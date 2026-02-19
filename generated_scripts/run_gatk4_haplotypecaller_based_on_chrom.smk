# [METADATA]
# TOOL_NAME = gatk4
# VERSION = 4.4.0.0
# THREADS = 1

rule gatk4:
    input:
        SeqID = ""
        Chromosome = ""
        BamDir = ""
        ResultDir = ""
        ReferenceFasta = "/storage/references_and_index/hg38/fasta/cbNIPT/hg38.fa"
        DbsnpVcf = "/storage/references_and_index/hg38/vcf/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
    output:
        OutGvcf = "[ResultDir]/[SeqID].[Chromosome].gvcf.gz"
        OutVcf = "[ResultDir]/[SeqID].[Chromosome].vcf.gz"
    params:
        InputSuffix = "analysisReady"
        singularity_bin = "singularity"
        gatk_bin = "gatk"
        sif = "/storage/images/gatk-4.4.0.0.sif"
        bind = "/storage,/data"
        Threads = "4"
        Memory = "32g"
        Ploidy = "2"
        TmpDir = "[ResultDir]/tmp/[SeqID]_[Chromosome]"
        IncludeNonVariant = "false"
    threads: 1
    shell:
        """
        {params.singularity_bin} exec -B {params.bind} {params.sif} {params.gatk_bin} --java-options '-XX:ParallelGCThreads={threads} -Xmx{params.Memory} -Djava.io.tmpdir={params.TmpDir}'  HaplotypeCaller  -R {input.ReferenceFasta}  -I {input.BamDir}/{input.SeqID}.{params.InputSuffix}.bam  -L {input.Chromosome}  -ploidy {params.Ploidy}  -stand-call-conf 30  --dbsnp {input.DbsnpVcf}  -O {output.OutGvcf}  -ERC GVCF &&
        {params.singularity_bin} exec -B {params.bind} {params.sif} {params.gatk_bin} --java-options '-XX:ParallelGCThreads={threads} -Xmx{params.Memory} -Djava.io.tmpdir={params.TmpDir}'  GenotypeGVCFs  --include-non-variant-sites {params.IncludeNonVariant}  -R {input.ReferenceFasta}  -V {output.OutGvcf}  -O {output.OutVcf}
        """