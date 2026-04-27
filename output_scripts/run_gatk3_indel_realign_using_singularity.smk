# [METADATA]
# TOOL_NAME = gatk3
# VERSION = 3.8
# THREADS = 1

rule gatk3:
    input:
        SeqID = ""
        BamDir = ""
        qcResDir = ""
        ReferenceFasta = ""
        TargetBed = ""
        KnownIndel1 = ""
        KnownIndel2 = ""
    output:
        target_intervals = "[qcResDir]/[SeqID].[OutputSuffix].realignertargetcreator.intervals"
        realigned_bam = "[BamDir]/[SeqID].[OutputSuffix].bam"
    params:
        InputSuffix = "merged.dup.marked"
        OutputSuffix = "merged.dup.marked.realign"
        singularity_bin = "singularity"
        java_bin = "java"
        gatk_jar = "/usr/GenomeAnalysisTK.jar"
        sif = "/storage/images/gatk-3.8-1.sif"
        bind = "/storage,/data"
        Threads = "1"
        xmx_mb = "16384"
    threads: 1
    shell:
        """
        {params.singularity_bin} exec -B {params.bind} {params.sif} {params.java_bin} -Xmx{params.xmx_mb}m -jar {params.gatk_jar} -T RealignerTargetCreator -R {input.ReferenceFasta} -L {input.TargetBed} -I {input.BamDir}/{input.SeqID}.{params.InputSuffix}.bam -o {output.target_intervals} -known {input.KnownIndel1} -known {input.KnownIndel2} -nt {threads} && {params.singularity_bin} exec -B {params.bind} {params.sif} {params.java_bin} -Xmx{params.xmx_mb}m -jar {params.gatk_jar} -T IndelRealigner -R {input.ReferenceFasta} -targetIntervals {output.target_intervals} -known {input.KnownIndel1} -known {input.KnownIndel2} -I {input.BamDir}/{input.SeqID}.{params.InputSuffix}.bam -o {output.realigned_bam}
        """