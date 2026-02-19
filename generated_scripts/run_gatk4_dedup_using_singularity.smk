# [METADATA]
# TOOL_NAME = gatk4
# VERSION = 4.4.0.0
# THREADS = 1

rule gatk4:
    input:
        SeqID = ""
        BamDir = ""
        qcResDir = ""
    output:
        dedup_bam = "[BamDir]/[SeqID].[InputSuffix].[OutputSuffix].bam"
        dedup_bai = "[BamDir]/[SeqID].[InputSuffix].[OutputSuffix].bam.bai"
        metrics_txt = "[qcResDir]/[SeqID].[InputSuffix].[OutputSuffix].metrics.txt"
    params:
        InputSuffix = "sorted"
        OutputSuffix = "dedup"
        singularity_bin = "singularity"
        gatk_bin = "gatk"
        sif = "/storage/images/gatk-4.4.0.0.sif"
        bind = "/storage,/data"
        Threads = "14"
        xmx_mb = "16384"
        remove_seq_dups = "true"
        create_index = "true"
        other_md_args = ""
    threads: 1
    shell:
        """
        {params.singularity_bin} exec -B {params.bind} {params.sif} {params.gatk_bin} MarkDuplicates --java-options '-XX:ParallelGCThreads={threads} -Xmx{params.xmx_mb}m' --INPUT {input.BamDir}/{input.SeqID}.{params.InputSuffix}.bam --OUTPUT {output.dedup_bam} --METRICS_FILE {output.metrics_txt} --CREATE_INDEX {params.create_index} --REMOVE_SEQUENCING_DUPLICATES {params.remove_seq_dups} {params.other_md_args}
        """