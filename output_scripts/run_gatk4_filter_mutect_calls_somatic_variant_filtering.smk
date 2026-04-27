# [METADATA]
# TOOL_NAME = gatk4_filter_mutect_calls
# VERSION = 4.4.0.0
# THREADS = 1

rule gatk4_filter_mutect_calls:
    input:
        SeqID = ""
        vcfDir = ""
        qcResDir = ""
        ReferenceFasta = ""
        TargetInterval = ""
        ContaminationTable = ""
    output:
        filtered_vcf = "[vcfDir]/[SeqID].mutect2.[Suffix]filtered.vcf"
    params:
        Suffix = ""
        singularity_bin = "singularity"
        gatk4_sif = "/storage/images/gatk-4.4.0.0.sif"
        bind = "/storage,/data"
        xmx_mb = "16384"
        Threads = "14"
        extra_args = ""
    threads: 1
    shell:
        """
        {params.singularity_bin} exec -B {params.bind} {params.gatk4_sif} gatk FilterMutectCalls --java-options '-XX:ParallelGCThreads={threads} -Xmx{params.xmx_mb}m' -V {input.vcfDir}/{input.SeqID}.mutect2.{params.Suffix}vcf -L {input.TargetInterval} --reference {input.ReferenceFasta} --contamination-table {input.ContaminationTable} --ob-priors {input.qcResDir}/{input.SeqID}.{params.Suffix}read-orientation-model.tar.gz -O {output.filtered_vcf} {params.extra_args}
        """