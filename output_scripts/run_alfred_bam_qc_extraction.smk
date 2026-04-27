# [METADATA]
# TOOL_NAME = alfred
# VERSION = 0.2.6
# THREADS = 1

rule alfred:
    input:
        SeqID = ""
        BamDir = ""
        qcResDir = ""
        ReferenceFasta = ""
        TargetBed = ""
        InputSuffix = "analysisReady"
    output:
        alfred_raw_tsv = "[qcResDir]/[SeqID].alfred.qc.tsv.gz"
        chr_map_stats = "[qcResDir]/[SeqID].alfred.chr.map.stats.txt"
        target_coverage = "[qcResDir]/[SeqID].alfred.target.coverage.txt"
    params:
        singularity_bin = "singularity"
        alfred_sif = "/storage/images/alfred-0.2.6.sif"
        alfred_bin = "/opt/alfred/bin/alfred"
        bind = "/storage,/data"
    threads: 1
    shell:
        """
        {params.singularity_bin} exec -B {params.bind} {params.alfred_sif} {params.alfred_bin} qc --reference {input.ReferenceFasta} --bed {input.TargetBed} --outfile {output.alfred_raw_tsv} {input.BamDir}/{input.SeqID}.{input.InputSuffix}.bam && zgrep '^CM' {output.alfred_raw_tsv} | cut -f 2- > {output.chr_map_stats} && zgrep '^TC' {output.alfred_raw_tsv} | cut -f 2- > {output.target_coverage}
        """