# [METADATA]
# TOOL_NAME = sob_detector
# VERSION = 1.0.4
# THREADS = 1

rule sob_detector:
    input:
        SeqID = ""
        vcfDir = ""
        BamDir = ""
    output:
        bias_filtered_vcf = "[vcfDir]/[SeqID].mutect2.bias.filtered.vcf"
    params:
        InputVcfSuffix = "filtered"
        InputBamSuffix = "analysisReady"
        singularity_bin = "singularity"
        sob_sif = "/storage/images/sobdetector.sif"
        sob_bin = "SOBDetector"
        bind = "/storage,/data"
        minBaseQ = "20"
        minMapQ = "20"
        only_passed = "false"
        extra_args = ""
    threads: 1
    shell:
        """
        {params.singularity_bin} exec -B {params.bind} {params.sob_sif} {params.sob_bin}  --input-type VCF --input-variants {input.vcfDir}/{input.SeqID}.mutect2.{params.InputVcfSuffix}.vcf --input-bam {input.BamDir}/{input.SeqID}.{params.InputBamSuffix}.bam --output-variants {output.bias_filtered_vcf} --minBaseQuality {params.minBaseQ} --minMappingQuality {params.minMapQ} --only-passed {params.only_passed} {params.extra_args}
        """