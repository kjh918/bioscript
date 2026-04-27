# [METADATA]
# TOOL_NAME = gatk4_contamination_estimation
# VERSION = 4.4.0.0
# THREADS = 1

rule gatk4_contamination_estimation:
    input:
        SeqID = ""
        BamDir = ""
        qcResDir = ""
        ReferenceFasta = ""
        TargetInterval = ""
        VcfGnomad = ""
    output:
        pileup_table = "[qcResDir]/[SeqID].targeted_sequencing.table"
        contamination_table = "[qcResDir]/[SeqID].contamination.table"
    params:
        InputSuffix = "analysisReady"
        singularity_bin = "singularity"
        gatk4_sif = "/storage/images/gatk-4.4.0.0.sif"
        bind = "/storage,/data"
        Threads = "14"
        pileup_xmx = "8192"
        contam_xmx = "16384"
    threads: 1
    shell:
        """
        {params.singularity_bin} exec -B {params.bind} {params.gatk4_sif} gatk GetPileupSummaries --java-options '-XX:ParallelGCThreads={threads} -Xmx{params.pileup_xmx}m' -I {input.BamDir}/{input.SeqID}.{params.InputSuffix}.bam -V {input.VcfGnomad} -L {input.TargetInterval} -R {input.ReferenceFasta} -O {output.pileup_table} && {params.singularity_bin} exec -B {params.bind} {params.gatk4_sif} gatk CalculateContamination --java-options '-XX:ParallelGCThreads={threads} -Xmx{params.contam_xmx}m' -I {output.pileup_table} -O {output.contamination_table}
        """