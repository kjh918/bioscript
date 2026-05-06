# [METADATA]
# TOOL_NAME = star_alignment
# VERSION = 2.7.11a
# THREADS = 1

rule star_alignment:
    input:
        SeqID = ""
        FastqDir = ""
        BamDir = ""
        StarIndex = ""
    output:
        aligned_bam = "[BamDir]/[SeqID].Aligned.sortedByCoord.out.bam"
        aligned_transcriptome_bam = "[BamDir]/[SeqID].Aligned.toTranscriptome.out.bam"
        gene_counts = "[BamDir]/[SeqID].ReadsPerGene.out.tab"
        mapping_log = "[BamDir]/[SeqID].Log.final.out"
    params:
        Threads = "12"
        InputSuffix = ".trimmed"
        outSAMtype = "BAM SortedByCoordinate"
        quantMode = "TranscriptomeSAM"
        outSAMattributes = "NH HI AS nM XS NM"
        chimSegmentMin = "10"
        twopassMode = "Basic"
        outFilterMismatchNmax = "10"
        outSAMunmapped = "Within"
        singularity_bin = "singularity"
        star_sif = "/storage/images/star-2.7.11.sif"
        bind = "/storage,/data"
    threads: 1
    shell:
        """
        {params.singularity_bin} exec -B {params.bind} {params.star_sif} STAR --runThreadN {threads} --genomeDir {input.StarIndex} --readFilesIn {input.FastqDir}/{input.SeqID}{params.InputSuffix}_R1.fastq.gz {input.FastqDir}/{input.SeqID}{params.InputSuffix}_R2.fastq.gz --readFilesCommand zcat --quantMode {params.quantMode} --twopassMode {params.twopassMode} --chimSegmentMin {params.chimSegmentMin} --outFilterMismatchNmax {params.outFilterMismatchNmax} --outFileNamePrefix {input.BamDir}/{input.SeqID}. --outSAMtype {params.outSAMtype} --outSAMattributes {params.outSAMattributes} --outSAMunmapped {params.outSAMunmapped}
        """