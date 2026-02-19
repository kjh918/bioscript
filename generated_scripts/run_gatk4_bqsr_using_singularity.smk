# [METADATA]
# TOOL_NAME = gatk4
# VERSION = 4.4.0.0
# THREADS = 1

rule gatk4:
    input:
        SeqID = ""
        BamDir = ""
        qcResDir = ""
        ReferenceFasta = ""
        KnownSnp = ""
        KnownIndel1 = ""
        KnownIndel2 = ""
    output:
        recal_table = "[qcResDir]/[SeqID].[InputSuffix].recal_table.txt"
        recal_bam = "[BamDir]/[SeqID].[OutputSuffix].bam"
        recal_bai = "[BamDir]/[SeqID].[OutputSuffix].bam.bai"
    params:
        InputSuffix = "dedup"
        OutputSuffix = "recal"
        singularity_bin = "singularity"
        gatk_bin = "gatk"
        sif = "/storage/images/gatk-4.4.0.0.sif"
        bind = "/storage,/data"
        Threads = "14"
        xmx_mb = "16384"
    threads: 1
    shell:
        """
        {params.singularity_bin} exec -B {params.bind} {params.sif} {params.gatk_bin} BaseRecalibrator --java-options '-XX:ParallelGCThreads={threads} -Xmx{params.xmx_mb}m' --input {input.BamDir}/{input.SeqID}.{params.InputSuffix}.bam --reference {input.ReferenceFasta} --output {output.recal_table} --known-sites {input.KnownSnp} --known-sites {input.KnownIndel1} --known-sites {input.KnownIndel2} && {params.singularity_bin} exec -B {params.bind} {params.sif} {params.gatk_bin} ApplyBQSR --java-options '-XX:ParallelGCThreads={threads} -Xmx{params.xmx_mb}m' --input {input.BamDir}/{input.SeqID}.{params.InputSuffix}.bam --bqsr-recal-file {output.recal_table} --output {output.recal_bam} --create-output-bam-index true
        """