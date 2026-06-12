

### Lab-Genomic File Treat -------------------------------------------------------------------------
    #Rscript /storage/home/jhkim/Projects/NGS/WES/v1/NGS_service.LabGenomicsFastq.R \
    #    --WES_BatchID WES_26_01 --LabGenomicsDataDir /data/ngs/lab_genomics.ngs_ruo_service/2025-12-30/251226_gencurix
###        
### CREATE PROJECT ---------------------------------------------------------------------------------
    #Rscript /storage/home/jhkim/Projects/NGS/WES/v1/NGS_service.CreateProject.R \
    #    --BASE_DIR "/data/wes" \
    #    --BATCH_ID "WES_26_01" \
    #    --CLIENT_ID 11 \
    #    --SAMPLE_INFO_FILE "WES_26_01_seqid_info.xlsx" \
    #    --ORDER_ID "GCX-C11-260331" \
    #    --DP_RD_FROM_FASTQ TRUE \
    #    --SEQ_TYPE "DNA" \
    #    --THREADS 15 \
    #    --WES_NGS_LIB "twist.exome.2.0" \
    #    --WES_RunGermlineVariant TRUE
###
### Run SGE Sample Processing ----------------------------------------------------------------------
#    qsub -q all.q@ngsnode1 /storage/home/kangsm/myScripts/Analysis/142_GCX-C15-251126_WES_26_01/WES_26_01_wes_SeqID_Processing_TumorOnly_Analysis_500X.sh
#    qsub -q all.q@ngsmaster /storage/home/kangsm/myScripts/Analysis/142_GCX-C15-251126_WES_26_01/WES_26_01_wes_SeqID_Processing_TumorOnly_Analysis_100X.sh
####
#### Individual Matching ----------------------------------------------------------------------------
    #Rscript /storage/home/jhkim/Projects/NGS/WES/v1/WES_run.Individual_Matching_Rev.R \
    #    --baseDir /data/wes \
    #    --seqFolder WES_26_01 \
    #    --assembly hg19 \
    #    --importDB TRUE
###
### Run SGE Matched-Analysis
#    qsub -q all.q@ngsnode1 /data/wes/WES_26_01/meta/WES_26_01_wes_SampleGroup_MatchedNormal_Analysis.sh
####

#### ANALYSIS 
    #cp /storage/home/jhkim/Projects/NGS/WES/v1/DEFAULT_wes_report_config.yaml /data/wes/WES_26_01/meta/wes_report_config.yaml
    #ln -s /data/wes/WES_26_01/meta/wes_report_config.yaml /storage/home/kangsm/myScripts/Analysis/142_GCX-C15-251126_WES_26_01/wes_report_config.yaml
    
    #Rscript /storage/home/jhkim/Projects/NGS/WES/v1/WES_run.ReportAnalysis_OrganoidScience.R \
    #    --BASE_DIR "/data/wes" --SEQ_FOLDER "WES_26_01" --PRESET_RDS "preset.genes_organoid.science.rds" \
    #    --CLIENT_ID 11 --ORDER_ID "GCX-C11-260331" --NGS_LIB "twist" --TONLY_CALL_MODE Tonly
####

#### Cumulative oncoplots
    #Rscript /storage/home/jhkim/Projects/NGS/WES/v1/NGS_service.DrawCumulativeOncoplots.R \
    #    --BASE_DIR /data/wes --SEQ_FOLDER "WES_26_01" --CLIENT_ID 11
####
#### preset gene variants excel file
    #Rscript /storage/home/jhkim/Projects/NGS/WES/v1/NGS_service.ExportPresetGeneVariantsAsExcel.R \
    #    --BASE_DIR /data/wes \
    #    --SEQ_FOLDER "WES_26_01" 
####
#### data export
    /storage/home/jhkim/Projects/NGS/WES/v1/NGS_service.ExportData_v2.sh \
        --baseDir /data/wes \
        --seqFolder "WES_26_01" \
        --checksum md5
####

#    ID_LIST_FILE=/data/wes/WES_26_01/meta/WES_26_01_export_data.list
#    EXPORT_DIR=/data/wes/WES_26_01/export_data
#    DATA_DIR=/data/wes/WES_26_01

#    for IDS in $(seq 1 $(cat $ID_LIST_FILE | wc -l))
#    do
#        TASK_ID=$(sed -n -e "${IDS} p" ${ID_LIST_FILE})
#        SEQ_ID=$(echo ${TASK_ID} | awk '{print $2}')
#        SAMPLE_ID=$(echo ${TASK_ID} | awk '{print $3}')
#        EXPORT_TONLY=$(echo ${TASK_ID} | awk '{print $4}')
#        EXPORT_TS_N=$(echo ${TASK_ID} | awk '{print $5'})
#        EXPORT_ORG_N=$(echo ${TASK_ID} | awk '{print $6}')
#        EXPORT_GERMLINE=$(echo ${TASK_ID} | awk '{print $7}')
#        EXPORT_TONLY_GF=$(echo ${TASK_ID} | awk '{print $8}')
#        TAG_TS_N="mutect2.NT.bias.filtered"

#        echo -e "|--->>> ${SEQ_ID} matched normal tissue call annotated MAF copying..."
#        MAF_FILE=${DATA_DIR}/${SEQ_ID}/vcf/${SEQ_ID}.${TAG_TS_N}.vep.annotated.maf
#        if [ -f ${MAF_FILE} ]; then 
#            cp ${MAF_FILE} ${EXPORT_DIR}/maf/${SAMPLE_ID}.matched.normal.TISSUE.maf
#            sed -i -e "s/${SEQ_ID}/${SAMPLE_ID}/g" ${EXPORT_DIR}/maf/${SAMPLE_ID}.matched.normal.TISSUE.maf
#        fi

#        echo -e "|--->>> ${SEQ_ID} matched normal tissue call somatic variants MAF copying..."
#        SOMATIC_MAF_FILE=${DATA_DIR}/${SEQ_ID}/vcf/${SEQ_ID}.${TAG_TS_N}.somatic.variants.only.maf
#        if [ -f ${SOMATIC_MAF_FILE} ]; then 
#            cp ${SOMATIC_MAF_FILE} ${EXPORT_DIR}/maf/${SAMPLE_ID}.matched.normal.TISSUE.somatic.variants.only.maf
#            sed -i -e "s/${SEQ_ID}/${SAMPLE_ID}/g" ${EXPORT_DIR}/maf/${SAMPLE_ID}.matched.normal.TISSUE.somatic.variants.only.maf
#        fi
#    done
####



