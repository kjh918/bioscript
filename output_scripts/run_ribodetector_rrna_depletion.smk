# [METADATA]
# TOOL_NAME = ribodetector
# VERSION = 0.3.3
# THREADS = 1

rule ribodetector:
    input:
        SeqID = ""
        FastqDir = ""
        OutputDir = ""
    output:
        non_rrna_r1 = "[OutputDir]/[SeqID].nonrrna.1.fq"
        non_rrna_r2 = "[OutputDir]/[SeqID].nonrrna.2.fq"
    params:
        Threads = "10"
        ReadLen = "100"
        ChunkSize = "256"
        ExcludeMode = "rrna"
        InputSuffix = "fastq.gz"
    threads: 1
    shell:
        """
        ribodetector_cpu -t {threads} -l {params.ReadLen} -i {input.FastqDir}/{input.SeqID}_R1.{params.InputSuffix} {input.FastqDir}/{input.SeqID}_R2.{params.InputSuffix} -e {params.ExcludeMode} --chunk_size {params.ChunkSize} -o {output.non_rrna_r1} {output.non_rrna_r2}
        """