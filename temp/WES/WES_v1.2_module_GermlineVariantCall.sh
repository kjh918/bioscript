
SEQ_FOLDER="23WES01"
SEQ_ID="23WES01_02"
LOG_NAME="23WES01_02.germline.processing.log"

DATA_DIR=/data/wes/${SEQ_FOLDER}

mkdir -p ${DATA_DIR}/${SEQ_ID}/gvcf

GVCF_DIR=${DATA_DIR}/${SEQ_ID}/gvcf
QC_DIR=${DATA_DIR}/${SEQ_ID}/qcfiles
BAM_DIR=${DATA_DIR}/${SEQ_ID}/bam
TMP_DIR=${DATA_DIR}/${SEQ_ID}/tmp

#---| REFERENCES
    REF_DIR=/storage/references_and_index/hg19
    BWA_INDEX=${REF_DIR}/fasta/human_g1k_v37_decoy.fasta
    GENOME_FASTA=${REF_DIR}/fasta/human_g1k_v37_decoy.fasta
    BAIT_BED=${REF_DIR}/bed/Twist_WES_2.0/probe.bed
    TARGET_BED=${REF_DIR}/bed/Twist_WES_2.0/target.bed
    BAIT_INTERVAL=${REF_DIR}/bed/Twist_WES_2.0/probe.interval_list
    TARGET_INTERVAL=${REF_DIR}/bed/Twist_WES_2.0/target.interval_list
  ##
    VCF_DBSNP=${REF_DIR}/vcf/dbsnp_138.b37.vcf.gz
    VCF_1KG_SNP=${REF_DIR}/vcf/1000G_phase1.snps.high_confidence.b37.vcf.gz
    VCF_GNOMAD=${REF_DIR}/vcf/af-only-gnomad.raw.sites.grch37.vcf.gz
    VCF_INDEL=${REF_DIR}/vcf/Homo_sapiens_assembly19.known_indels.vcf.gz
    VCF_1KG_INDEL=${REF_DIR}/vcf/Mills_and_1000G_gold_standard.indels.b37.vcf.gz
    VCF_PON=${REF_DIR}/vcf/gatk-practice-mutect2-panel-of-normal.vcf.gz
    VCF_HAPMAP=${REF_DIR}/vcf/hapmap_3.3.b37.vcf.gz
#---/

LOGFILE=/data/wes/${SEQ_FOLDER}/${SEQ_ID}/${LOG_NAME}

#----
    IMAGE_GATK=/storage/images/gatk-4.4.0.0.sif

singularity exec -B /storage,/data ${IMAGE_GATK} gatk --java-options "-XX:ConcGCThreads=10 -Xmx32G" HaplotypeCaller \
    -R ${GENOME_FASTA} \
    -I ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.realign.recal.bam \
    -O ${GVCF_DIR}/${SEQ_ID}.gvcf.gz \
    -ERC GVCF \
    -pairHMM AVX_LOGLESS_CACHING_OMP \
    --native-pair-hmm-threads 10 \
    --min-base-quality-score 20 \
    -L ${TARGET_INTERVAL} \
    -D ${VCF_DBSNP} \
    -G StandardAnnotation \
    -G AS_StandardAnnotation

singularity exec -B /storage,/data ${IMAGE_GATK} gatk --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads=10 -Xmx16G" GenotypeGVCFs \
    -R ${GENOME_FASTA} \
    -V ${GVCF_DIR}/${SEQ_ID}.gvcf.gz \
    -O ${GVCF_DIR}/${SEQ_ID}.genotype.vcf.gz


singularity exec -B /storage,/data ${IMAGE_GATK} gatk --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads=10 -Xmx16G" CNNScoreVariants \
    -I ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.realign.recal.bam \
    -V ${GVCF_DIR}/${SEQ_ID}.genotype.vcf.gz \
    -R ${GENOME_FASTA} \
    -O ${GVCF_DIR}/${SEQ_ID}.genotype.cnn.vcf \
    -tensor-type read_tensor


singularity exec -B /storage,/data ${IMAGE_GATK} gatk --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads=10 -Xmx16G" FilterVariantTranches \
    -V ${GVCF_DIR}/${SEQ_ID}.genotype.cnn.vcf \
    --resource ${VCF_HAPMAP} \
    --resource ${VCF_1KG_INDEL} \
    --info-key CNN_2D \
    --snp-tranche 99.95 \
    --indel-tranche 99.5 \
    --invalidate-previous-filters \
    -O ${GVCF_DIR}/${SEQ_ID}.genotype.cnn.filtered.vcf

    