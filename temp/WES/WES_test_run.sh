
#!/bin/bash

#---|
    ARGS=$(getopt -a -o s:h: --long sampleID:,help -- "$@" )
    VALID_ARGS=$?
    if [ "$VALID_ARGS" != "0" ]; then 
        usage >&2 
        exit 2
    fi

    eval set -- "$ARGS"
    while :
    do
    case "$1" in
        -s | --sampleID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  SEQ_ID=$2 ; shift 2 ;;
            esac  ;;
        -h | --help )
        usage >&2 ; exit 2 ;;   
        --) shift ; break ;;
        *)  usage >&2 ; exit 2 ;;
    esac
    done
#---|
    SEQ_FOLDER=/data/wes/23WES02
    FASTQ_DIR=${SEQ_FOLDER}/${SEQ_ID}/fastq
    BAM_DIR=${SEQ_FOLDER}/${SEQ_ID}/bam
    VCF_DIR=${SEQ_FOLDER}/${SEQ_ID}/vcf
    QC_DIR=${SEQ_FOLDER}/${SEQ_ID}/qcfiles
    TMP_DIR=${SEQ_FOLDER}/${SEQ_ID}/tmp
#---| REFERENCES
    REF_DIR=/storage/references_and_index/gcx.hg19
    BWA_INDEX=${REF_DIR}/bwa-index/gcx.hg19
    GENOME_FASTA=${REF_DIR}/gcx.hg19.fasta
    BAIT_BED=${REF_DIR}/TwistExome2.0/gcx.hg19.twist.exome2.probe.bed
    TARGET_BED=${REF_DIR}/TwistExome2.0/gcx.hg19.twist.exome2.target.bed
    BAIT_INTERVAL=${REF_DIR}/TwistExome2.0/gcx.hg19.twist.exome2.probe.interval_list
    TARGET_INTERVAL=${REF_DIR}/TwistExome2.0/gcx.hg19.twist.exome2.target.interval_list
  ##
    VCF_DBSNP=${REF_DIR}/vcf/gcx.hg19.dbsnp138.vcf.gz
    VCF_1KG_SNP=${REF_DIR}/vcf/gcx.hg19.1000G_phase1.snps.high_confidence.vcf.gz
    VCF_GNOMAD=${REF_DIR}/vcf/gcx.hg19.af-only-gnomad.raw.sites.vcf.gz
    VCF_INDEL=${REF_DIR}/vcf/gcx.hg19.known_indels.vcf.gz
    VCF_1KG_INDEL=${REF_DIR}/vcf/gcx.hg19.Mills_and_1000G_gold_standard.indels.sites.vcf.gz
    VCF_PON=${REF_DIR}/vcf/gcx.hg19.mutect2-exome-panel.vcf.gz
#---/

    /usr/lib/jvm/java-1.8.0-openjdk-amd64/bin/java -Xmx8G -jar /storage/apps/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
        -T RealignerTargetCreator \
        -R ${GENOME_FASTA} \
        -L ${TARGET_BED} \
        -known ${VCF_INDEL} \
        -known ${VCF_1KG_INDEL} \
        -I ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.recal.bam \
        -o ${QC_DIR}/${SEQ_ID}.merged.dup.marked.recal.bam.realignertargetcreator.intervals
  ##
    /usr/lib/jvm/java-1.8.0-openjdk-amd64/bin/java -Xmx16G -jar /storage/apps/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
        -T IndelRealigner \
        -R ${GENOME_FASTA} \
        -targetIntervals ${QC_DIR}/${SEQ_ID}.merged.dup.marked.recal.bam.realignertargetcreator.intervals \
        -known ${VCF_INDEL} \
        -known ${VCF_1KG_INDEL} \
        -I ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.recal.bam \
        -o ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.recal.realign.bam  

    java -XX:+UseParallelGC -XX:ParallelGCThreads=10 -Xmx16g -jar /storage/apps/bin/gatk-local.jar CollectSequencingArtifactMetrics \
        --INPUT ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.recal.realign.bam  \
        --OUTPUT ${QC_DIR}/${SEQ_ID}.realign --FILE_EXTENSION .txt \
        --REFERENCE_SEQUENCE ${GENOME_FASTA}

    java -XX:+UseParallelGC -XX:ParallelGCThreads=10 -Xmx16g -jar /storage/apps/bin/gatk-local.jar CollectAlignmentSummaryMetrics \
        --INPUT ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.recal.realign.bam \
        --OUTPUT ${QC_DIR}/${SEQ_ID}.merged.dup.marked.recal.realign.bam.AlignSummary.metrics.txt

    for SEQ_ID in 3-100X 3-200X 3-400X 17-00X 17-200X 17-400X tissue-200X
    do
        BAM_DIR=${SEQ_FOLDER}/${SEQ_ID}/bam
        VCF_DIR=${SEQ_FOLDER}/${SEQ_ID}/vcf
        QC_DIR=${SEQ_FOLDER}/${SEQ_ID}/qcfiles
    
        java -XX:+UseParallelGC -XX:ParallelGCThreads=10 -Xmx16g -jar /storage/apps/bin/gatk-local.jar CollectHsMetrics \
            --INPUT ${SEQ_FOLDER}/${SEQ_ID}/bam/${SEQ_ID}.merged.dup.marked.recal.realign.bam \
            --OUTPUT ${SEQ_FOLDER}/${SEQ_ID}/qcfiles/${SEQ_ID}.merged.dup.marked.recal.realign.HS.metrics.txt \
            --BAIT_INTERVALS ${BAIT_INTERVAL} \
            --TARGET_INTERVALS ${TARGET_INTERVAL}
    done


    /storage/apps/mosdepth-0.3.5/mosdepth --threads 20 \
        --by ${TARGET_BED} \
        --mapq 20 \
        ${QC_DIR}/${SEQ_ID}.realign \
        ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.recal.realign.bam 
    
    /storage/apps/alfred-0.2.7/bin/alfred qc \
        --reference ${GENOME_FASTA} \
        --bed ${TARGET_BED} \
        --outfile ${QC_DIR}/${SEQ_ID}.reslign.alfred.qc.tsv.gz \
        ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.recal.realign.bam 
    ##
    zgrep ^CM ${QC_DIR}/${SEQ_ID}.reslign.alfred.qc.tsv.gz | cut -f 2- > ${QC_DIR}/${SEQ_ID}.reslign.alfred.chr.map.stats.txt
    zgrep ^TC ${QC_DIR}/${SEQ_ID}.reslign.alfred.qc.tsv.gz | cut -f 2- > ${QC_DIR}/${SEQ_ID}.reslign.alfred.target.coverage.txt

    java -Xmx32G -jar /storage/apps/bin/gatk-local.jar Mutect2 \
        --reference ${GENOME_FASTA} \
        --intervals ${TARGET_BED} \
        --input ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.recal.realign.bam  \
        --germline-resource ${VCF_GNOMAD} \
        --panel-of-normals ${VCF_PON} \
        --f1r2-tar-gz ${QC_DIR}/${SEQ_ID}.realign.f1r2.tar.gz \
        --output ${VCF_DIR}/${SEQ_ID}.realign.mutect2.vcf

    java -XX:+UseParallelGC -XX:ParallelGCThreads=10 -Xmx8G -jar /storage/apps/bin/gatk-local.jar LearnReadOrientationModel \
        -I ${QC_DIR}/${SEQ_ID}.realign.f1r2.tar.gz \
        -O ${QC_DIR}/${SEQ_ID}.realign.read-orientation-model.tar.gz

    java -XX:+UseParallelGC -XX:ParallelGCThreads=10 -Xmx8G -jar /storage/apps/bin/gatk-local.jar GetPileupSummaries \
        -I ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.recal.realign.bam \
        -O ${QC_DIR}/${SEQ_ID}.realign.targeted_sequencing.table \
        -V ${VCF_GNOMAD} \
        -L ${TARGET_BED} \
        -R ${GENOME_FASTA}

    mv ${QC_DIR}/${SEQ_ID}.relaign.targeted_sequencing.table ${QC_DIR}/${SEQ_ID}.realign.targeted_sequencing.table 

    java -XX:+UseParallelGC -XX:ParallelGCThreads=10 -Xmx16G -jar /storage/apps/bin/gatk-local.jar CalculateContamination \
        -I ${QC_DIR}/${SEQ_ID}.realign.targeted_sequencing.table \
        -O ${QC_DIR}/${SEQ_ID}.realign.targeted_sequencing.contamination.table

    java -XX:+UseParallelGC -XX:ParallelGCThreads=10 -Xmx16G -jar /storage/apps/bin/gatk-local.jar FilterMutectCalls \
        -V ${VCF_DIR}/${SEQ_ID}.realign.mutect2.vcf \
        -L ${TARGET_BED} \
        --reference ${GENOME_FASTA} \
        --contamination-table ${QC_DIR}/${SEQ_ID}.realign.targeted_sequencing.contamination.table \
        --ob-priors ${QC_DIR}/${SEQ_ID}.realign.read-orientation-model.tar.gz \
        -O ${VCF_DIR}/${SEQ_ID}.realign.mutect2.filtered.vcf

    DATA_DIR=/data/wes/23WES02
    CUSTOM_DIR=/storage/DB
    CACHE_DIR=/storage/DB/ensembl-vep-cache/vep-110

    SEQ_IDS="3-100X 3-200X 3-400X 17-100X 17-200X 17-400X tissue-200X"

    parallel -j 7 -k "singularity exec -B $CUSTOM_DIR:/opt/vep/.vep,$DATA_DIR:/vep-input,$CACHE_DIR:/vep-cache /storage/images/vep-110.sif /opt/vep/src/ensembl-vep/vep \
        --cache --dir_cache /vep-cache --force_overwrite --merged --stats_text --species homo_sapiens --assembly GRCh37 \
        --fork 30 --show_ref_allele --uploaded_allele --allele_number --numbers \
        --everything \
        --plugin AlphaMissense,file=/opt/vep/.vep/alpha-missense-prediction/AlphaMissense_b37.tsv.gz \
        --input_file /vep-input/{}/vcf/{}.realign.mutect2.filtered.vcf \
        --vcf --output_file /vep-input/{}/vcf/{}.realign.mutect2.filtered.vep.annotate.vcf " ::: ${SEQ_IDS}


    parallel -j 7 -k "perl /storage/apps/vcf2maf-1.6.21/vcf2maf.pl \
        --inhibit-vep \
        --input-vcf /data/wes/23WES02/{}/vcf/{}.realign.mutect2.filtered.vep.annotate.vcf \
        --output-maf /data/wes/23WES02/{}/vcf/{}.realign.mutect2.filtered.vep.annotate.maf \
        --tumor-id {} \
        --retain-info DP,FILTER \
        --retain-ann am_class,am_pathogenicity" ::: ${SEQ_IDS}


singularity exec -B $CUSTOM_DIR:/opt/vep/.vep,$DATA_DIR:/vep-input,$CACHE_DIR:/vep-cache /storage/images/vep-110.sif /opt/vep/src/ensembl-vep/vep \
        --cache --dir_cache /vep-cache --force_overwrite --merged --stats_text --species homo_sapiens --assembly GRCh37 \
        --fork 35 --show_ref_allele --uploaded_allele --allele_number --numbers \
        --everything \
        --plugin AlphaMissense,file=/opt/vep/.vep/alpha-missense-prediction/AlphaMissense_b37.tsv.gz \
        --input_file /vep-input/${SEQ_ID}/vcf/${SEQ_ID}.realign.mutect2.filtered.vcf