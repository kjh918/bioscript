# [METADATA]
# TOOL_NAME = fastq_screen
# VERSION = 0.15.3
# THREADS = 1

rule fastq_screen:
    input:
        SeqID = ""
        RawFastqDir = ""
        qcResDir = ""
    output:
        screen_txt = "[qcResDir]/[SeqID]_R1_screen.txt"
        screen_png = "[qcResDir]/[SeqID]_R1_screen.png"
        screen_html = "[qcResDir]/[SeqID]_R1_screen.html"
    params:
        singularity_bin = "singularity"
        fastq_screen_bin = "fastq_screen"
        sif = "/storage/images/fastqScreen-0.15.3.sif"
        bind = "/storage,/data"
        aligner = "bowtie2"
        config_path = "/storage/home/kangsm/runScripts/NGS_config.FastqScreen.conf"
        Threads = "15"
        extra_args = "--force"
    threads: 1
    shell:
        """
        {params.singularity_bin} exec -B {params.bind} {params.sif} {params.fastq_screen_bin} --aligner {params.aligner} --conf {params.config_path} --outdir {input.qcResDir} --threads {threads} {params.extra_args} {input.RawFastqDir}/{input.SeqID}_R1.fastq.gz
        """