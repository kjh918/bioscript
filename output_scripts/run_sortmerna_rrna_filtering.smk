# [METADATA]
# TOOL_NAME = sortmerna
# VERSION = 4.3.7
# THREADS = 1

rule sortmerna:
    input:
        SeqID = ""
        FastqDir = ""
        qcResDir = ""
        RefArgs = ""
        IndexDir = ""
    output:
        non_rrna_r1 = "[qcResDir]/[SeqID]_1.non_rRNA.fastq.gz"
        non_rrna_r2 = "[qcResDir]/[SeqID]_2.non_rRNA.fastq.gz"
        sortmerna_log = "[qcResDir]/[SeqID].sortmerna.log"
    params:
        Threads = "8"
        InputSuffix = "fastq.gz"
        paired_cmd = "--paired_in --out2"
        extra_args = ""
    threads: 1
    shell:
        """
        sortmerna {input.RefArgs} --reads {input.FastqDir}/{input.SeqID}_1.{params.InputSuffix} --reads {input.FastqDir}/{input.SeqID}_2.{params.InputSuffix} --threads {threads} --workdir . --aligned {input.qcResDir}/rRNA_reads --fastx --other {input.qcResDir}/non_rRNA_reads {params.paired_cmd} {params.extra_args} && mv {input.qcResDir}/non_rRNA_reads_fwd.f*q.gz {output.non_rrna_r1} && mv {input.qcResDir}/non_rRNA_reads_rev.f*q.gz {output.non_rrna_r2} && mv {input.qcResDir}/rRNA_reads.log {output.sortmerna_log}
        """