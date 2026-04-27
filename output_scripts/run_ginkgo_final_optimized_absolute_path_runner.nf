// [METADATA]
// TOOL_NAME = Ginkgo_Final_Optimized
// THREADS = 1

process GINKGO_FINAL_OPTIMIZED {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), val(BinSize), val(ReadLength), val(Threads), val(GinkgoHomeDir), path(BamToBedDir), val(NormalDir), val(ResultBaseDir), val(Genome)

    output:
    path "*", emit: BinMeth
    path "*", emit: WorkDir
    path "*", emit: BinUnsorted
    path "*", emit: BinFile
    path "*", emit: BinCount
    path "*", emit: ProcessR
    path "*", emit: ReclustR
    path "*", emit: CNVCaller

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def BinSize = params.BinSize ?: ""
    def ReadLength = params.ReadLength ?: ""
    def Threads = params.Threads ?: "1"
    def GinkgoHomeDir = params.GinkgoHomeDir ?: "/storage/apps/ginkgo"
    def BamToBedDir = params.BamToBedDir ?: "/data/cbNIPT/bamToBeds"
    def NormalDir = params.NormalDir ?: "/data/cbNIPT/bamToBeds"
    def ResultBaseDir = params.ResultBaseDir ?: "/data/cbNIPT/ginkgo_analysis"
    def Genome = params.Genome ?: "hg38"
    def NormalSamples = params.NormalSamples ?: "Normal_01.bed.gz Normal_02.bed.gz Normal_03.bed.gz"
    def R_Args_Base = params.R_Args_Base ?: "status.xml data 2 [BinMeth] ward euclidean 1 refDummy.bed_mapped 0 ploidyDummy.txt 0 1"
    def R_Args_Reclust = params.R_Args_Reclust ?: "status.xml [BinMeth] ward euclidean 0 ploidyDummy.txt 0"
    """
    mkdir -p ${WorkDir} && ln -sf ${NormalDir}/*.bed.gz ${WorkDir}/ && ln -sf ${BamToBedDir}/${SeqID}.bed.gz ${WorkDir}/ &&
find ${WorkDir} -maxdepth 1 -name '*.gz' | xargs -I {} bash -c 'N=$(basename '{}' | sed 's/\.bed.*//'); ${BinUnsorted} ${BinFile} ${BinCount} <(zcat '{}' | awk '{if(\$1!~/^chr/) print \'chr\'\$0; else print \$0}') $N '{}'_mapped' &&
find ${WorkDir} -maxdepth 1 -name '*.bed' ! -name '*.gz' | xargs -I {} bash -c 'N=$(basename '{}' | sed 's/\.bed.*//'); ${BinUnsorted} ${BinFile} ${BinCount} <(awk '{if(\$1!~/^chr/) print \'chr\'\$0; else print \$0}' '{}') $N '{}'_mapped' &&
paste ${WorkDir}/*_mapped > ${WorkDir}/data &&
${ProcessR} ${GinkgoHomeDir}/genomes/${Genome}/original ${WorkDir} ${R_Args_Base} && ${ReclustR} ${GinkgoHomeDir}/genomes/${Genome}/original ${WorkDir} ${R_Args_Reclust} && ${CNVCaller} ${WorkDir}/SegCopy ${WorkDir}/CNV1 ${WorkDir}/CNV2
    """
}