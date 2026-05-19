# [METADATA]
# TOOL_NAME = bedtools_genomecov_to_bigwig
# VERSION = v2.27.1
# THREADS = 1

rule bedtools_genomecov_to_bigwig:
    input:
        SeqID = ""
        InputFile = ""
        GenomeSizes = ""
        Extension = ""
        OutputDir = ""
    output:
    params:
        Threads = "4"
        InputFlag = "-ibam"
        GenomeSizesFlag = ""
        ScaleArgs = "-bg"
        SortCmd = "| LC_ALL=C sort --buffer-size=8G -k1,1 -k2,2n"
        singularity_bin = "singularity"
        bedtools_sif = "/storage/images/bedtools-2.27.1.sif"
        bigwig_sif = "/storage/images/ucsc-bedgraphtobigwig-445.sif"
        bind = "/storage,/data"
    threads: 1
    shell:
        """
        {params.singularity_bin} exec -B {params.bind} {params.bedtools_sif} bedtools genomecov {params.InputFlag} {input.InputFile} {params.GenomeSizesFlag} {params.ScaleArgs} {params.SortCmd} --parallel={threads} > {input.OutputDir}/{input.SeqID}.{input.Extension} && {params.singularity_bin} exec -B {params.bind} {params.bigwig_sif} bedGraphToBigWig {input.OutputDir}/{input.SeqID}.{input.Extension} {input.GenomeSizes} {input.OutputDir}/{input.SeqID}.bw
        """