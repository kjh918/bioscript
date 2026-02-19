# [METADATA]
# TOOL_NAME = bwa_picard
# VERSION = bwa_0.7.17-picard_3.1.0
# THREADS = 1

rule bwa_picard:
    input:
        SeqID = ""
        TrimFastqDir = ""
        BamDir = ""
        TmpDir = ""
        ReferenceFasta = ""
        ReadGroupID = ""
        ReadGroupPlatform = ""
        ReadGroupLibrary = ""
        ReadGroupCenter = ""
    output:
        unmapped_bam = "[BamDir]/[SeqID].[InputSuffix].u.bam"
        aligned_sam = "[BamDir]/[SeqID].[InputSuffix].bwa_raw.sam"
        primary_bam = "[BamDir]/[SeqID].[OutputSuffix].bam"
        primary_bai = "[BamDir]/[SeqID].[OutputSuffix].bai"
    params:
        InputSuffix = "trimmed"
        OutputSuffix = "primary"
        java_bin = "java"
        picard_jar = "/storage/apps/bin/picard.jar"
        singularity_bin = "singularity"
        bwa_sif = "/storage/images/bwa-0.7.17.sif"
        bind = "/storage,/data"
        xmx_mb = "16384"
        Threads = "8"
        bwa_bin = "bwa"
        mark_short_split = "-M"
        soft_clipping = "-Y"
        clipping_penalty = "-L 50,50"
        other_args = ""
        mba_strategy = "MostDistant"
        mba_attributes = "XS"
        mba_orientations = "FR --EXPECTED_ORIENTATIONS RF"
        mba_other_flags = "--CREATE_INDEX true --MAX_INSERTIONS_OR_DELETIONS -1 --CLIP_ADAPTERS false"
    threads: 1
    shell:
        """
        {params.java_bin} -XX:ParallelGCThreads={threads} -Xmx{params.xmx_mb}m -jar {params.picard_jar} FastqToSam --FASTQ {input.TrimFastqDir}/{input.SeqID}.{params.InputSuffix}_R1.fastq.gz --FASTQ2 {input.TrimFastqDir}/{input.SeqID}.{params.InputSuffix}_R2.fastq.gz --SAMPLE_NAME {input.SeqID} --OUTPUT {output.unmapped_bam} --READ_GROUP_NAME {input.ReadGroupID} --PLATFORM {input.ReadGroupPlatform} --LIBRARY_NAME {input.ReadGroupLibrary} --SEQUENCING_CENTER {input.ReadGroupCenter} --TMP_DIR {input.TmpDir} && {params.singularity_bin} exec -B {params.bind} {params.bwa_sif} {params.bwa_bin} mem  {params.mark_short_split} {params.soft_clipping} {params.clipping_penalty} {params.other_args}  -t {threads} -R '@RG\tID:{input.ReadGroupID}\tPL:{input.ReadGroupPlatform}\tLB:{input.ReadGroupLibrary}\tSM:{input.SeqID}\tCN:{input.ReadGroupCenter}' {input.ReferenceFasta} {input.TrimFastqDir}/{input.SeqID}.{params.InputSuffix}_R1.fastq.gz {input.TrimFastqDir}/{input.SeqID}.{params.InputSuffix}_R2.fastq.gz > {output.aligned_sam} && {params.java_bin} -XX:ParallelGCThreads={threads} -Xmx{params.xmx_mb}m -jar {params.picard_jar} MergeBamAlignment --UNMAPPED_BAM {output.unmapped_bam} --ALIGNED_BAM {output.aligned_sam} --REFERENCE_SEQUENCE {input.ReferenceFasta} --OUTPUT {output.primary_bam} --PRIMARY_ALIGNMENT_STRATEGY {params.mba_strategy} --ATTRIBUTES_TO_RETAIN {params.mba_attributes} --EXPECTED_ORIENTATIONS {params.mba_orientations} {params.mba_other_flags}
        """