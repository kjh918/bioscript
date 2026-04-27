# [METADATA]
# TOOL_NAME = gatk4_learn_read_orientation_model
# VERSION = 4.4.0.0
# THREADS = 1

rule gatk4_learn_read_orientation_model:
    input:
        SeqID = ""
        qcResDir = ""
    output:
        read_orientation_model = "[qcResDir]/[SeqID].[OutputSuffix].read-orientation-model.tar.gz"
    params:
        InputSuffix = ""
        OutputSuffix = ""
        singularity_bin = "singularity"
        gatk4_sif = "/storage/images/gatk-4.4.0.0.sif"
        bind = "/storage,/data"
        xmx_mb = "8192"
        Threads = "14"
    threads: 1
    shell:
        """
        {params.singularity_bin} exec -B {params.bind} {params.gatk4_sif} gatk LearnReadOrientationModel --java-options '-XX:ParallelGCThreads={threads} -Xmx{params.xmx_mb}m' -I {input.qcResDir}/{input.SeqID}.{params.InputSuffix}f1r2.tar.gz -O {output.read_orientation_model}
        """