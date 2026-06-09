# [METADATA]
# TOOL_NAME = cbNIPT_MakeBins
# VERSION = 1.0.0
# THREADS = 1

rule cbnipt_makebins:
    input:
        ReferenceFasta = ""
        MappabilityBW = ""
    output:
        OutBinFile = ""
    params:
        BinSize = "100000"
        MinMappability = "0.9"
        MinGC = "0.3"
        MaxGC = "0.7"
        IncludeSexChrom = "--IncludeSexChrom"
        python_bin = "python"
        script_path = ""
    threads: 1
    shell:
        """
        {params.python_bin} {params.script_path} make-bins --ReferenceFasta {input.ReferenceFasta} --MappabilityBW {input.MappabilityBW} --OutBinFile {output.OutBinFile} --BinSize {params.BinSize} --MinMappability {params.MinMappability} --MinGC {params.MinGC} --MaxGC {params.MaxGC} {params.IncludeSexChrom}
        """