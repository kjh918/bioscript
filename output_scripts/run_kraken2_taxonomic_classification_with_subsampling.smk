# [METADATA]
# TOOL_NAME = Kraken2
# VERSION = 2.17.1
# THREADS = 8

rule kraken2:
    input:
        SeqID = ""
        TrimFastqDir = ""
        KrakenDB = ""
        OutDir = ""
    output:
        sub_r1 = "[OutDir]/[SeqID].trimmed_sub_R1.fastq.gz"
        sub_r2 = "[OutDir]/[SeqID].trimmed_sub_R2.fastq.gz"
        kraken_out = "[OutDir]/[SeqID].kraken"
        report_out = "[OutDir]/[SeqID].kreport"
        plot_out = "[OutDir]/[SeqID]_krona.html"
    params:
        seqkit_bin = "/storage/apps/seqkit-2.8.2/seqkit"
        SampleSeed = "100"
        SampleNum = "100000"
        R1_suffix = ".trimmed_R1.fastq.gz"
        R2_suffix = ".trimmed_R2.fastq.gz"
        kraken2_bin = "/storage/home/jhkim/Apps/kraken2/bin/kraken2"
        KrakenDB = "/storage/home/jhkim/Apps/kraken2/bin/kraken2_db"
        krona_bin = "/storage/home/jhkim/Apps/Krona/KronaTools/bin/ktImportTaxonomy"
        Threads = "8"
    threads: 8
    shell:
        """
        {params.seqkit_bin} sample -s {params.SampleSeed} -n {params.SampleNum} {input.TrimFastqDir}/{input.SeqID}{params.R1_suffix} | gzip > {output.sub_r1} && {params.seqkit_bin} sample -s {params.SampleSeed} -n {params.SampleNum} {input.TrimFastqDir}/{input.SeqID}{params.R2_suffix} | gzip > {output.sub_r2} && {params.kraken2_bin} --db {input.KrakenDB} --threads {threads} --report {output.report_out} --paired {output.sub_r1} {output.sub_r2} > {output.kraken_out} && {params.krona_bin} -m 4 -t 3 {output.kraken_out} -o {output.plot_out}
        """