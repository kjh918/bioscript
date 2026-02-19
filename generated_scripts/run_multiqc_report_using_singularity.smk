# [METADATA]
# TOOL_NAME = multiqc
# VERSION = 1.16
# THREADS = 1

rule multiqc:
    input:
        SeqID = ""
        qcResDir = ""
    output:
        report_html = "[qcResDir]/[SeqID].QC.Results.html"
    params:
        singularity_bin = "singularity"
        multiqc_bin = "multiqc"
        sif = "/storage/images/multiqc-1.16.sif"
        bind = "/storage,/data"
        mqc_config = "/storage/home/kangsm/runScripts/NGS_config.MultiQC_Custom.yaml"
        mqc_filename = "[SeqID].QC.Results"
        mqc_args = "--force --data-dir"
    threads: 1
    shell:
        """
        {params.singularity_bin} exec -B {params.bind} {params.sif} {params.multiqc_bin} {params.mqc_args} --filename {params.mqc_filename} --outdir {input.qcResDir} --config {params.mqc_config} {input.qcResDir}
        """