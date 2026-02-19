# [METADATA]
# TOOL_NAME = bwa_mem
# VERSION = 0.7.17
# THREADS = 1

rule bwa_mem:
    input:
        SeqID = ""
        TrimFastqDir = ""
        BamDir = ""
        ReferenceFasta = ""
        ReadGroupID = ""
        ReadGroupPlatform = ""
        ReadGroupLibrary = ""
        ReadGroupCenter = ""
    output:
        primary_bam = "[BamDir]/[SeqID].[OutputSuffix].bam"
    params:
        InputSuffix = "trimmed"
        OutputSuffix = "primary"
        singularity_bin = "singularity"
        bwa_bin = "bwa"
        samtools_bin = "samtools"
        sif = "/storage/images/bwa-0.7.17.sif"
        bind = "/storage,/data"
        Threads = "8"
        mark_short_split = "-M"
        soft_clipping = "-Y"
        clipping_penalty = "-L 50,50"
        other_args = ""
    threads: 1
    shell:
        """
        {params.singularity_bin} exec -B {params.bind} {params.sif} {params.bwa_bin} mem  {params.mark_short_split} {params.soft_clipping} {params.clipping_penalty} {params.other_args}  -t {threads} -R '@RG\tID:{input.ReadGroupID}\tPL:{input.ReadGroupPlatform}\tLB:{input.ReadGroupLibrary}\tSM:{input.SeqID}\tCN:{input.ReadGroupCenter}' {input.ReferenceFasta} {input.TrimFastqDir}/{input.SeqID}.{params.InputSuffix}_R1.fastq.gz {input.TrimFastqDir}/{input.SeqID}.{params.InputSuffix}_R2.fastq.gz | {params.samtools_bin} view -@ {threads} -bS -o {output.primary_bam} -
        """