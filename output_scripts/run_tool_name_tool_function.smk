# [METADATA]
# TOOL_NAME = TOOL_NAME
# VERSION = TASK_VERSION
# THREADS = 8

rule tool_name:
    input:
        SeqID = ""
        InputFileDir = ""
        qcResDir = ""
        OutputDir = ""
    output:
        OutputSummary = ""
    params:
        InputSuffix = ""
        OutputSuffix = ""
        singularity_bin = "singularity"
        bind = "/storage,/data"
        xmx_mb = "16384"
        Threads = "8"
        extra_args = ""
    threads: 8
    shell:
        """
        {params.singularity_bin} exec -B {params.bind}  {params.extra_args}
        """