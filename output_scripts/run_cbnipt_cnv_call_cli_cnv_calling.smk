# [METADATA]
# TOOL_NAME = cbnipt_cnv_call_cli
# VERSION = 1.0.0
# THREADS = 1

rule cbnipt_cnv_call_cli:
    input:
        SeqID = ""
        BamPath = ""
        ReferenceFasta = ""
        OutDir = ""
    output:
        OutDir = ""
    params:
        PythonBin = "/storage/home/jhkim/Apps/Python-3.11.13/python"
        CliScript = "/storage/home/jhkim/scripts/bioscript/manual/cbnipt_manual_cnv_call/scripts/cli.py"
        Threads = "4"
    threads: 1
    shell:
        """
        {params.PythonBin} {params.CliScript} --SeqID {input.SeqID} --BamPath {input.BamPath} --ReferenceFasta {input.ReferenceFasta} --OutDir {input.OutDir} --Threads {threads}
        """