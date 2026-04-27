// [METADATA]
// TOOL_NAME = gatk4
// THREADS = 1

process GATK4 {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), val(Chromosome), path(BamDir), val(ResultDir), val(ReferenceFasta), path(DbsnpVcf)

    output:
    path "*", emit: OutGvcf
    path "*", emit: OutVcf

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def Chromosome = params.Chromosome ?: ""
    def BamDir = params.BamDir ?: ""
    def ResultDir = params.ResultDir ?: ""
    def ReferenceFasta = params.ReferenceFasta ?: "/storage/references_and_index/hg38/fasta/cbNIPT/hg38.fa"
    def DbsnpVcf = params.DbsnpVcf ?: "/storage/references_and_index/hg38/vcf/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
    def InputSuffix = params.InputSuffix ?: "analysisReady"
    def singularity_bin = params.singularity_bin ?: "singularity"
    def gatk_bin = params.gatk_bin ?: "gatk"
    def sif = params.sif ?: "/storage/images/gatk-4.4.0.0.sif"
    def bind = params.bind ?: "/storage,/data"
    def Threads = params.Threads ?: "4"
    def Memory = params.Memory ?: "32g"
    def Ploidy = params.Ploidy ?: "2"
    def TmpDir = params.TmpDir ?: "[ResultDir]/tmp/[SeqID]_[Chromosome]"
    def IncludeNonVariant = params.IncludeNonVariant ?: "false"
    """
    ${singularity_bin} exec -B ${bind} ${sif} ${gatk_bin} --java-options '-XX:ParallelGCThreads=${Threads} -Xmx${Memory} -Djava.io.tmpdir=${TmpDir}'  HaplotypeCaller  -R ${ReferenceFasta}  -I ${BamDir}/${SeqID}.${InputSuffix}.bam  -L ${Chromosome}  -ploidy ${Ploidy}  -stand-call-conf 30  --dbsnp ${DbsnpVcf}  -O ${OutGvcf}  -ERC GVCF && 
${singularity_bin} exec -B ${bind} ${sif} ${gatk_bin} --java-options '-XX:ParallelGCThreads=${Threads} -Xmx${Memory} -Djava.io.tmpdir=${TmpDir}'  GenotypeGVCFs  --include-non-variant-sites ${IncludeNonVariant}  -R ${ReferenceFasta}  -V ${OutGvcf}  -O ${OutVcf}
    """
}