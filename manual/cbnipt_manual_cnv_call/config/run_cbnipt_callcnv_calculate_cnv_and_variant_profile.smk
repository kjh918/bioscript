# [METADATA]
# TOOL_NAME = cbNIPT_CallCNV
# VERSION = 1.0.0
# THREADS = 4

rule cbnipt_callcnv:
    input:
        SeqID = ""
        BamPath = ""
        AnnotatedBins = ""
    output:
        OutDir = ""
    params:
        MinMapQ = "20"
        LowessFrac = "0.2"
        BaselinePloidy = "2"
        SmoothWindow = "5"
        MinDepth = "2.0"
        MinCoverage = "0.5"
        MaskLowerQ = "0.01"
        MaskUpperQ = "0.99"
        SegPenalty = "10.0"
        Threads = "4"
        python_bin = "python"
        script_path = ""
    threads: 4
    shell:
        """
        {params.python_bin} {params.script_path} call-cnv --SeqID {input.SeqID} --BamPath {input.BamPath} --AnnotatedBins {input.AnnotatedBins} --OutDir {output.OutDir} --MinMapQ {params.MinMapQ} --LowessFrac {params.LowessFrac} --BaselinePloidy {params.BaselinePloidy} --SmoothWindow {params.SmoothWindow} --MinDepth {params.MinDepth} --MinCoverage {params.MinCoverage} --MaskLowerQ {params.MaskLowerQ} --MaskUpperQ {params.MaskUpperQ} --SegPenalty {params.SegPenalty} --Threads {threads}
        """