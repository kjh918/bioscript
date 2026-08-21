# [METADATA]
# TOOL_NAME = SRAToolkit
# VERSION = 3.0.10
# THREADS = 8

rule sratoolkit:
    input:
        SraAccession = ""
        OutDir = ""
        TmpDir = "[OutDir]/tmp"
    output:
        sra_dir = "[TmpDir]/sra/[SraAccession]"
        r1_fastq = "[OutDir]/[SraAccession]_R1.fastq"
        r2_fastq = "[OutDir]/[SraAccession]_R2.fastq"
        r1_fastq_gz = "[OutDir]/[SraAccession]_R1.fastq.gz"
        r2_fastq_gz = "[OutDir]/[SraAccession]_R2.fastq.gz"
    params:
        prefetch_bin = "/storage/apps/sratoolkit.3.0.10/bin/prefetch"
        fasterq_dump_bin = "/storage/apps/sratoolkit.3.0.10/bin/fasterq-dump"
        Threads = "8"
        MaxSize = "100G"
        SplitMode = "--split-files"
    threads: 8
    shell:
        """
        mkdir -p {input.OutDir} {input.TmpDir}/sra {input.TmpDir}/fasterq && {params.prefetch_bin} --max-size {params.MaxSize} --output-directory {input.TmpDir}/sra {input.SraAccession} && {params.fasterq_dump_bin} {params.SplitMode} --threads {threads} --temp {input.TmpDir}/fasterq --outdir {input.OutDir} {output.sra_dir} && mv {input.OutDir}/{input.SraAccession}_1.fastq {output.r1_fastq} && mv {input.OutDir}/{input.SraAccession}_2.fastq {output.r2_fastq} && gzip -f {output.r1_fastq} && gzip -f {output.r2_fastq}
        """