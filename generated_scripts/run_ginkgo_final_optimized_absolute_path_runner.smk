# [METADATA]
# TOOL_NAME = Ginkgo_Final_Optimized
# VERSION = 1.0.0
# THREADS = 1

rule ginkgo_final_optimized:
    input:
        SeqID = ""
        BinSize = ""
        ReadLength = ""
        Threads = "1"
        GinkgoHomeDir = "/storage/apps/ginkgo"
        BamToBedDir = "/data/cbNIPT/bamToBeds"
        NormalDir = "/data/cbNIPT/bamToBeds"
        ResultBaseDir = "/data/cbNIPT/ginkgo_analysis"
        Genome = "hg38"
    output:
        BinMeth = "variable_[BinSize]000_[ReadLength]_bwa"
        WorkDir = "[ResultBaseDir]/[SeqID]/[SeqID]_[BinSize]kb"
        BinUnsorted = "[GinkgoHomeDir]/scripts/binUnsorted"
        BinFile = "[GinkgoHomeDir]/genomes/[Genome]/original/[BinMeth]"
        BinCount = "$(wc -l < [BinFile])"
        ProcessR = "[GinkgoHomeDir]/scripts/process.R"
        ReclustR = "[GinkgoHomeDir]/scripts/reclust.R"
        CNVCaller = "[GinkgoHomeDir]/scripts/CNVcaller"
    params:
        NormalSamples = "Normal_01.bed.gz Normal_02.bed.gz Normal_03.bed.gz"
        R_Args_Base = "status.xml data 2 [BinMeth] ward euclidean 1 refDummy.bed_mapped 0 ploidyDummy.txt 0 1"
        R_Args_Reclust = "status.xml [BinMeth] ward euclidean 0 ploidyDummy.txt 0"
    threads: 1
    shell:
        """
        mkdir -p {output.WorkDir} && ln -sf {input.NormalDir}/*.bed.gz {output.WorkDir}/ && ln -sf {input.BamToBedDir}/{input.SeqID}.bed.gz {output.WorkDir}/ &&
        find {output.WorkDir} -maxdepth 1 -name '*.gz' | xargs -I {} bash -c 'N=$(basename '{}' | sed 's/\.bed.*//'); {output.BinUnsorted} {output.BinFile} {output.BinCount} <(zcat '{}' | awk '{if(\$1!~/^chr/) print \'chr\'\$0; else print \$0}') $N '{}'_mapped' &&
        find {output.WorkDir} -maxdepth 1 -name '*.bed' ! -name '*.gz' | xargs -I {} bash -c 'N=$(basename '{}' | sed 's/\.bed.*//'); {output.BinUnsorted} {output.BinFile} {output.BinCount} <(awk '{if(\$1!~/^chr/) print \'chr\'\$0; else print \$0}' '{}') $N '{}'_mapped' &&
        paste {output.WorkDir}/*_mapped > {output.WorkDir}/data &&
        {output.ProcessR} {input.GinkgoHomeDir}/genomes/{input.Genome}/original {output.WorkDir} {params.R_Args_Base} && {output.ReclustR} {input.GinkgoHomeDir}/genomes/{input.Genome}/original {output.WorkDir} {params.R_Args_Reclust} && {output.CNVCaller} {output.WorkDir}/SegCopy {output.WorkDir}/CNV1 {output.WorkDir}/CNV2
        """