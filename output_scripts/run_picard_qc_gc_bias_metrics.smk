# [METADATA]
# TOOL_NAME = picard
# VERSION = 3.1.0
# THREADS = 1

rule picard:
    input:
        SeqID = ""
        BamDir = ""
        qcResDir = ""
    output:
        out_bam = "[BamDir]/[SeqID].sorted.dedup.bam"
        metrics = "[qcResDir]/[SeqID].mark.duplicates.metrics.txt"
    params:
        java_bin = "java"
        picard_jar = "/storage/apps/bin/picard.jar"
        Threads = "14"
        Memory = "16384m"
        TmpDir = "/tmp"
        create_index = "true"
        remove_duplicates = "true"
    threads: 1
    shell:
        """
        {params.java_bin} -XX:ParallelGCThreads={threads} -Xmx{params.Memory} -jar {params.picard_jar} MarkDuplicates  INPUT={input.BamDir}/{input.SeqID}.sorted.bam  OUTPUT={output.out_bam}  METRICS_FILE={output.metrics}  CREATE_INDEX={params.create_index}  REMOVE_DUPLICATES={params.remove_duplicates}  TMP_DIR={params.TmpDir}
        """