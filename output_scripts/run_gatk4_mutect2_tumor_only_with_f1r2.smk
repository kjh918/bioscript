# [METADATA]
# TOOL_NAME = gatk4_mutect2
# VERSION = 4.4.0.0
# THREADS = 1

rule gatk4_mutect2:
    input:
        SeqID = ""
        BamDir = ""
        vcfDir = ""
        qcResDir = ""
        ReferenceFasta = ""
        TargetInterval = ""
        VcfGnomad = ""
        VcfPon = ""
    output:
        raw_vcf = "[vcfDir]/[SeqID].mutect2.vcf"
        f1r2_tar_gz = "[qcResDir]/[SeqID].f1r2.tar.gz"
        mutect_stats = "[vcfDir]/[SeqID].mutect2.vcf.stats"
    params:
        InputSuffix = "analysisReady"
        singularity_bin = "singularity"
        gatk4_sif = "/storage/images/gatk-4.4.0.0.sif"
        bind = "/storage,/data"
        xmx_mb = "32768"
        Threads = "14"
        pairHMM = "AVX_LOGLESS_CACHING_OMP"
        f1r2_max_depth = "2500"
        min_base_q = "20"
        extra_args = ""
    threads: 1
    shell:
        """
        {params.singularity_bin} exec -B {params.bind} {params.gatk4_sif} gatk Mutect2 --java-options '-XX:+UseParallelGC -XX:ParallelGCThreads={threads} -Xmx{params.xmx_mb}m' -pairHMM {params.pairHMM} --reference {input.ReferenceFasta} --intervals {input.TargetInterval} --input {input.BamDir}/{input.SeqID}.{params.InputSuffix}.bam --germline-resource {input.VcfGnomad} --panel-of-normals {input.VcfPon} --f1r2-tar-gz {output.f1r2_tar_gz} --f1r2-max-depth {params.f1r2_max_depth} --min-base-quality-score {params.min_base_q} --base-quality-score-threshold {params.min_base_q} --output {output.raw_vcf} {params.extra_args}
        """