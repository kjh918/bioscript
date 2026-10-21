
### hg19 ###
    BASE_DIR=/storage/home/kangsm/analysis/cnv_seqc2
    SEQ_FOLDER=hg19
    GENOME_ASSEMBLY="hg19"
    WES_LIB_KIT="sureselect.v6"
    PURECN_RESOURCES=/storage/home/kangsm/myDB/cnv_purecn
#--------------------------------------------------------------------------------------------------#

#---| SWTOOLS RUN |--------------------------------------------------------------------------------#   
    RUN_PURECN=/storage/apps/R/4.3.1/lib/R/library/PureCN/extdata
    RUN_SINGULARITY="singularity exec -B /storage,/data /storage/images"
    GATK3="${RUN_SINGULARITY}/gatk-3.8-1.sif java -Xmx16384m -jar /usr/GenomeAnalysisTK.jar"
#--------------------------------------------------------------------------------------------------#

#---| REFERENCES |---------------------------------------------------------------------------------#
    if [ "${GENOME_ASSEMBLY}" = "hg38" ]; then
        GENOME_FASTA=/storage/references_and_index/hg38/fasta/Homo_sapiens_assembly38.fasta
    else
        GENOME_FASTA=/storage/references_and_index/hg19/fasta/human_g1k_v37_decoy.fasta
    fi
    #--------------------------------------------------------------------------#
    SIMPLE_REPEAT=${PURECN_RESOURCES}/snp_blacklist/${GENOME_ASSEMBLY}_simple.repeat.bed
    INTERVALS=${PURECN_RESOURCES}/interval_files/${GENOME_ASSEMBLY}_${WES_LIB_KIT}.baits_intervals.txt
    PON_RDS=${PURECN_RESOURCES}/normalDB/${GENOME_ASSEMBLY}_${WES_LIB_KIT}/normalDB_${GENOME_ASSEMBLY}_${WES_LIB_KIT}_${GENOME_ASSEMBLY}.rds
#--------------------------------------------------------------------------------------------------#

    for TUMOR_SEQ_ID in FD2_N FD2_T FD3_N FD3_T NV1_N NV1_T NV2_N NV2_T NV3_N NV3_T
    do
        DATA_DIR=${BASE_DIR}/${SEQ_FOLDER}
        TUMOR_BAM_DIR=${DATA_DIR}/${TUMOR_SEQ_ID}/bam
        CNV_RES_DIR=${DATA_DIR}/${TUMOR_SEQ_ID}/cnv
        if [ ! -d ${CNV_RES_DIR} ]; then mkdir -p ${CNV_RES_DIR}; fi
        #--------------------------------------------------------------------------#  
        PURECN_RES_DIR=${CNV_RES_DIR}/purecn
        if [ ! -d ${PURECN_RES_DIR} ]; then mkdir -p ${PURECN_RES_DIR}; fi

        Rscript ${RUN_PURECN}/Coverage_rev.R \
            --out-dir ${PURECN_RES_DIR} \
            --bam ${TUMOR_BAM_DIR}/${TUMOR_SEQ_ID}.analysisReady.bam \
            --intervals ${INTERVALS} \
            --force 
    done

#--------------------------------------------------------------------------------------------------#
    INPUT_LIST_FILE=${BASE_DIR}/${SEQ_FOLDER}/meta/NT_match.list

    for SGE_TASK_ID in 2 3 4 5 6
    do
        SEQ_INFO=$(sed -n -e "${SGE_TASK_ID} p" ${INPUT_LIST_FILE})
        #----------------------------------------------------------------------#
        TUMOR_SEQ_ID=$(echo ${SEQ_INFO} | awk '{print $2}')
        NORMAL_SEQ_ID=$(echo ${SEQ_INFO} | awk '{print $1}')
        NORMAL_COVERAGE=${DATA_DIR}/${NORMAL_SEQ_ID}/cnv/purecn/${NORMAL_SEQ_ID}.analysisReady_coverage_loess.txt.gz
        TUMOR_VCF=${DATA_DIR}/${TUMOR_SEQ_ID}/vcf/${TUMOR_SEQ_ID}.mutect2.NT.purecn.filtered.vcf
        PURECN_RES_DIR=${DATA_DIR}/${TUMOR_SEQ_ID}/cnv/purecn
        CALLABLE_CDS_BED=${PURECN_RES_DIR}/${TUMOR_SEQ_ID}.minDP.20.status.CDS.bed
        #----------------------------------------------------------------------#
        Rscript ${RUN_PURECN}/PureCN.R \
            --out ${PURECN_RES_DIR}/NT_mode \
            --tumor ${PURECN_RES_DIR}/${TUMOR_SEQ_ID}.analysisReady_coverage_loess.txt.gz \
            --normal ${NORMAL_COVERAGE} \
            --sampleid ${TUMOR_SEQ_ID} \
            --vcf ${TUMOR_VCF} \
            --intervals ${INTERVALS} \
            --genome ${GENOME_ASSEMBLY} \
            --snp-blacklist ${SIMPLE_REPEAT}
        #----------------------------------------------------------------------#
        Rscript ${RUN_PURECN}/PureCN.R \
            --out ${PURECN_RES_DIR}/NT_mode2 \
            --tumor ${PURECN_RES_DIR}/${TUMOR_SEQ_ID}.analysisReady_coverage_loess.txt.gz \
            --normal ${NORMAL_COVERAGE} \
            --normaldb ${PON_RDS} \
            --fun-segmentation PSCBS \
            --sampleid ${TUMOR_SEQ_ID} \
            --vcf ${TUMOR_VCF} \
            --intervals ${INTERVALS} \
            --genome ${GENOME_ASSEMBLY} \
            --snp-blacklist ${SIMPLE_REPEAT}
        #----------------------------------------------------------------------#
        Rscript ${RUN_PURECN}/PureCN.R \
            --out ${PURECN_RES_DIR}/TumorOnly_mode \
            --tumor ${PURECN_RES_DIR}/${TUMOR_SEQ_ID}.analysisReady_coverage_loess.txt.gz \
            --normaldb ${PON_RDS} \
            --sampleid ${TUMOR_SEQ_ID} \
            --vcf ${TUMOR_VCF} \
            --intervals ${INTERVALS} \
            --genome ${GENOME_ASSEMBLY} \
            --snp-blacklist ${SIMPLE_REPEAT}
        #----------------------------------------------------------------------#
        Rscript ${RUN_PURECN}/Dx.R \
            --out ${PURECN_RES_DIR}/NT_mode/${TUMOR_SEQ_ID} \
            --rds ${PURECN_RES_DIR}/NT_mode/${TUMOR_SEQ_ID}.rds \
            --callable ${CALLABLE_CDS_BED} \
            --exclude ${SIMPLE_REPEAT} \
            --signatures \
            --force    
        #----------------------------------------------------------------------#
        Rscript ${RUN_PURECN}/Dx.R \
            --out ${PURECN_RES_DIR}/TumorOnly_mode/${TUMOR_SEQ_ID} \
            --rds ${PURECN_RES_DIR}/TumorOnly_mode/${TUMOR_SEQ_ID}.rds \
            --callable ${CALLABLE_CDS_BED} \
            --exclude ${SIMPLE_REPEAT} \
            --signatures \
            --force
        #----------------------------------------------------------------------#
        Rscript ${RUN_PURECN}/Dx.R \
            --out ${PURECN_RES_DIR}/NT_mode2/${TUMOR_SEQ_ID} \
            --rds ${PURECN_RES_DIR}/NT_mode2/${TUMOR_SEQ_ID}.rds \
            --callable ${CALLABLE_CDS_BED} \
            --exclude ${SIMPLE_REPEAT} \
            --signatures \
            --force    
        #----------------------------------------------------------------------#
        echo -e ">>>>>\n     ${TUMOR_SEQ_ID} DONE.  \n<<<<<"
        #----------------------------------------------------------------------#
    done



### hg38 ###
    BASE_DIR=/storage/home/kangsm/analysis/cnv_seqc2
    SEQ_FOLDER=hg38
    GENOME_ASSEMBLY="hg38"
    WES_LIB_KIT="sureselect.v6"
    PURECN_RESOURCES=/storage/home/kangsm/myDB/cnv_purecn
    RUN_PURECN=/storage/apps/R/4.3.1/lib/R/library/PureCN/extdata
    RUN_SINGULARITY="singularity exec -B /storage,/data /storage/images"
    GATK3="${RUN_SINGULARITY}/gatk-3.8-1.sif java -Xmx16384m -jar /usr/GenomeAnalysisTK.jar"
    if [ "${GENOME_ASSEMBLY}" = "hg38" ]; then
        GENOME_FASTA=/storage/references_and_index/hg38/fasta/Homo_sapiens_assembly38.fasta
    else
        GENOME_FASTA=/storage/references_and_index/hg19/fasta/human_g1k_v37_decoy.fasta
    fi
    #--------------------------------------------------------------------------#
    SIMPLE_REPEAT=${PURECN_RESOURCES}/snp_blacklist/${GENOME_ASSEMBLY}_simple.repeat.bed
    INTERVALS=${PURECN_RESOURCES}/interval_files/${GENOME_ASSEMBLY}_${WES_LIB_KIT}.baits_intervals.txt
    PON_RDS=${PURECN_RESOURCES}/normalDB/${GENOME_ASSEMBLY}_${WES_LIB_KIT}/normalDB_${GENOME_ASSEMBLY}_${WES_LIB_KIT}_${GENOME_ASSEMBLY}.rds
    
    for TUMOR_SEQ_ID in FD1_N FD1_T
    do
        DATA_DIR=${BASE_DIR}/${SEQ_FOLDER}
        TUMOR_BAM_DIR=${DATA_DIR}/${TUMOR_SEQ_ID}/bam
        CNV_RES_DIR=${DATA_DIR}/${TUMOR_SEQ_ID}/cnv
        if [ ! -d ${CNV_RES_DIR} ]; then mkdir -p ${CNV_RES_DIR}; fi
        #--------------------------------------------------------------------------#  
        PURECN_RES_DIR=${CNV_RES_DIR}/purecn
        if [ ! -d ${PURECN_RES_DIR} ]; then mkdir -p ${PURECN_RES_DIR}; fi

        Rscript ${RUN_PURECN}/Coverage_rev.R \
            --out-dir ${PURECN_RES_DIR} \
            --bam ${TUMOR_BAM_DIR}/${TUMOR_SEQ_ID}.analysisReady.bam \
            --intervals ${INTERVALS} \
            --force 
    done


    INPUT_LIST_FILE=${BASE_DIR}/${SEQ_FOLDER}/meta/NT_match.list

    for SGE_TASK_ID in 2 3 4 5 6
    do
        SEQ_INFO=$(sed -n -e "${SGE_TASK_ID} p" ${INPUT_LIST_FILE})
        #----------------------------------------------------------------------#
        TUMOR_SEQ_ID=$(echo ${SEQ_INFO} | awk '{print $2}')
        NORMAL_SEQ_ID=$(echo ${SEQ_INFO} | awk '{print $1}')
        NORMAL_COVERAGE=${DATA_DIR}/${NORMAL_SEQ_ID}/cnv/purecn/${NORMAL_SEQ_ID}.analysisReady_coverage_loess.txt.gz
        TUMOR_VCF=${DATA_DIR}/${TUMOR_SEQ_ID}/vcf/${TUMOR_SEQ_ID}.mutect2.NT.purecn.filtered.vcf
        PURECN_RES_DIR=${DATA_DIR}/${TUMOR_SEQ_ID}/cnv/purecn
        CALLABLE_CDS_BED=${PURECN_RES_DIR}/${TUMOR_SEQ_ID}.minDP.20.status.CDS.bed
        #----------------------------------------------------------------------#
        Rscript ${RUN_PURECN}/PureCN.R \
            --out ${PURECN_RES_DIR}/NT_mode \
            --tumor ${PURECN_RES_DIR}/${TUMOR_SEQ_ID}.analysisReady_coverage_loess.txt.gz \
            --normal ${NORMAL_COVERAGE} \
            --sampleid ${TUMOR_SEQ_ID} \
            --vcf ${TUMOR_VCF} \
            --intervals ${INTERVALS} \
            --genome ${GENOME_ASSEMBLY} \
            --snp-blacklist ${SIMPLE_REPEAT} \
            --force
        #----------------------------------------------------------------------#
        Rscript ${RUN_PURECN}/PureCN.R \
            --out ${PURECN_RES_DIR}/NT_mode2 \
            --tumor ${PURECN_RES_DIR}/${TUMOR_SEQ_ID}.analysisReady_coverage_loess.txt.gz \
            --normal ${NORMAL_COVERAGE} \
            --normaldb ${PON_RDS} \
            --fun-segmentation PSCBS \
            --sampleid ${TUMOR_SEQ_ID} \
            --vcf ${TUMOR_VCF} \
            --intervals ${INTERVALS} \
            --genome ${GENOME_ASSEMBLY} \
            --snp-blacklist ${SIMPLE_REPEAT} \
            --force
        #----------------------------------------------------------------------#
        Rscript ${RUN_PURECN}/PureCN.R \
            --out ${PURECN_RES_DIR}/TumorOnly_mode \
            --tumor ${PURECN_RES_DIR}/${TUMOR_SEQ_ID}.analysisReady_coverage_loess.txt.gz \
            --normaldb ${PON_RDS} \
            --sampleid ${TUMOR_SEQ_ID} \
            --vcf ${TUMOR_VCF} \
            --intervals ${INTERVALS} \
            --genome ${GENOME_ASSEMBLY} \
            --snp-blacklist ${SIMPLE_REPEAT} \
            --force
        #----------------------------------------------------------------------#
        Rscript ${RUN_PURECN}/Dx.R \
            --out ${PURECN_RES_DIR}/NT_mode/${TUMOR_SEQ_ID} \
            --rds ${PURECN_RES_DIR}/NT_mode/${TUMOR_SEQ_ID}.rds \
            --callable ${CALLABLE_CDS_BED} \
            --exclude ${SIMPLE_REPEAT} \
            --signatures \
            --force    
        #----------------------------------------------------------------------#
        Rscript ${RUN_PURECN}/Dx.R \
            --out ${PURECN_RES_DIR}/TumorOnly_mode/${TUMOR_SEQ_ID} \
            --rds ${PURECN_RES_DIR}/TumorOnly_mode/${TUMOR_SEQ_ID}.rds \
            --callable ${CALLABLE_CDS_BED} \
            --exclude ${SIMPLE_REPEAT} \
            --signatures \
            --force
        #----------------------------------------------------------------------#
        Rscript ${RUN_PURECN}/Dx.R \
            --out ${PURECN_RES_DIR}/NT_mode2/${TUMOR_SEQ_ID} \
            --rds ${PURECN_RES_DIR}/NT_mode2/${TUMOR_SEQ_ID}.rds \
            --callable ${CALLABLE_CDS_BED} \
            --exclude ${SIMPLE_REPEAT} \
            --signatures \
            --force    
        #----------------------------------------------------------------------#
        echo -e ">>>>>\n     ${TUMOR_SEQ_ID} DONE.  \n<<<<<"
        #----------------------------------------------------------------------#
    done







