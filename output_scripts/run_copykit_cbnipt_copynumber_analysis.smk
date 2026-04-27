# [METADATA]
# TOOL_NAME = CopyKit
# VERSION = 1.0.0
# THREADS = 1

rule copykit:
    input:
        SeqID = ""
        NGS_DataBaseDir = ""
        ResultBaseDir = ""
    output:
        AnalysisRunDir = "[ResultBaseDir]/[SeqID]"
    params:
        InputSuffix = "analysisReady"
        mkdir_bin = "mkdir"
        rscript_bin = "Rscript"
        Rscript_path = "/storage/home/kangsm/runScripts/cbNIPT/GCX.cbNIPT_Run.CopyKit.Analysis.R"
        Threads = "1"
        BinSize = "220kb"
        GenomeVersion = "hg38"
        SamplePloidy = "2"
    threads: 1
    shell:
        """
        {params.mkdir_bin} -p {output.AnalysisRunDir} && ln -Tsf {input.NGS_DataBaseDir}/{input.SeqID}.{params.InputSuffix}.bam {output.AnalysisRunDir}/{input.SeqID}.bam && ln -Tsf {input.NGS_DataBaseDir}/{input.SeqID}.{params.InputSuffix}.bam.bai {output.AnalysisRunDir}/{input.SeqID}.bam.bai && {params.rscript_bin} {params.Rscript_path}  --SeqID {input.SeqID}  --AnalysisRunDir {output.AnalysisRunDir}  --BinSize {params.BinSize}  --Ploidy {params.SamplePloidy}  --GenomeVersion {params.GenomeVersion}
        """