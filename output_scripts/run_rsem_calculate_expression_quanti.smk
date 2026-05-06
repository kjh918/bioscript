# [METADATA]
# TOOL_NAME = rsem
# VERSION = 1.3.3
# THREADS = 1

rule rsem:
    input:
        SeqID = ""
        TranscriptomeBam = ""
        InputSuffix = ".Aligned.toTranscriptome.out.bam"
        RsemIndex = ""
        OutputDir = ""
    output:
    params:
        Threads = "12"
        PairedFlag = "--paired-end"
        NoBamOutput = "--no-bam-output"
        EstimateRspd = "--estimate-rspd"
        StrandOption = "none"
        ExtraArgs = "None"
        singularity_bin = "singularity"
        rsem_sif = "/storage/images/rsem-1.3.3.sif"
        bind = "/storage,/data"
    threads: 1
    shell:
        """
        {params.singularity_bin} exec -B {params.bind} {params.rsem_sif} rsem-calculate-expression {params.PairedFlag} {params.NoBamOutput} {params.EstimateRspd} {params.ExtraArgs} --num-threads {threads} --strandedness {params.StrandOption} --alignments {input.TranscriptomeBam}/{input.SeqID}{input.InputSuffix} {input.RsemIndex} {input.OutputDir}/{input.SeqID}
        """