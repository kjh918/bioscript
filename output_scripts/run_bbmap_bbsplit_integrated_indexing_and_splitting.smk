# [METADATA]
# TOOL_NAME = bbmap_bbsplit
# VERSION = 39.81
# THREADS = 1

rule bbmap_bbsplit:
    input:
        SeqID = ""
        Read1 = ""
        Read2 = ""
        PrimaryRef = ""
        OtherRefsArgs = ""
        OutputDir = ""
    output:
    params:
        MemoryMB = ""
        Threads = ""
        BBMapPath = ""
        BaseNamePattern = ""
        ExtraArgs = "ambiguous2=all"
        singularity_bin = "singularity"
        bbmap_sif = "/storage/home/jhkim/Apps/bbmap_39.81.sif"
        bind = "/storage,/data"
    threads: 1
    shell:
        """
        {params.singularity_bin} exec -B {params.bind} {params.bbmap_sif} bbsplit.sh -Xmx{params.MemoryMB}M ref_primary={input.PrimaryRef} {input.OtherRefsArgs} in1={input.Read1} in2={input.Read2} basename={params.BaseNamePattern} refstats={input.OutputDir}/{input.SeqID}.refstats.txt threads={threads} {params.ExtraArgs}
        """