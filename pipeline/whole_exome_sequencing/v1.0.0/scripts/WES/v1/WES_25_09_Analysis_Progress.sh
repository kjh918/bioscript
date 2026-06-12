

### Lab-Genomic File Treat -------------------------------------------------------------------------
    Rscript /storage/home/kangsm/runScripts/NGS_service.LabGenomicsFastq.R \
        --WES_BatchID WES_25_09 --LabGenomicsDataDir /data/ngs/lab_genomics.ngs_ruo_service/2025-12-30/251226_gencurix
###        
### CREATE PROJECT ---------------------------------------------------------------------------------
    Rscript /storage/home/kangsm/runScripts/NGS_service.CreateProject.R \
        --BASE_DIR "/data/wes" \
        --BATCH_ID "WES_25_09" \
        --CLIENT_ID 15 \
        --SAMPLE_INFO_FILE "WES_25_09_seqid_info.xlsx" \
        --ORDER_ID "GCX-C15-251126" \
        --DP_RD_FROM_FASTQ TRUE \
        --SEQ_TYPE "DNA" \
        --THREADS 15 \
        --WES_NGS_LIB "twist.exome.2.0" \
        --WES_RunGermlineVariant TRUE
###
### Run SGE Sample Processing ----------------------------------------------------------------------
    qsub -q all.q@ngsnode1 /storage/home/kangsm/myScripts/Analysis/142_GCX-C15-251126_WES_25_09/WES_25_09_wes_SeqID_Processing_TumorOnly_Analysis_500X.sh
    qsub -q all.q@ngsmaster /storage/home/kangsm/myScripts/Analysis/142_GCX-C15-251126_WES_25_09/WES_25_09_wes_SeqID_Processing_TumorOnly_Analysis_100X.sh
###
### Individual Matching ----------------------------------------------------------------------------
    Rscript /storage/home/kangsm/runScripts/WES_run.Individual_Matching_Rev.R \
        --baseDir /data/wes \
        --seqFolder WES_25_09 \
        --assembly hg19 \
        --importDB TRUE
###
### Run SGE Matched-Analysis
    qsub -q all.q@ngsnode1 /data/wes/WES_25_09/meta/WES_25_09_wes_SampleGroup_MatchedNormal_Analysis.sh
###

### ANALYSIS 
    cp /storage/home/kangsm/myScripts/Default_Scripts/DEFAULT_wes_report_config.yaml /data/wes/WES_25_09/meta/wes_report_config.yaml
    ln -s /data/wes/WES_25_09/meta/wes_report_config.yaml /storage/home/kangsm/myScripts/Analysis/142_GCX-C15-251126_WES_25_09/wes_report_config.yaml
    
    Rscript /storage/home/kangsm/runScripts/WES_run.ReportAnalysis_OrganoidScience.R \
        --BASE_DIR "/data/wes" --SEQ_FOLDER "WES_25_09" --PRESET_RDS "preset.genes_organoid.science.rds" \
        --CLIENT_ID 15 --ORDER_ID "GCX-C15-251126" --NGS_LIB "twist" --TONLY_CALL_MODE Tonly
###

### Cumulative oncoplots
    Rscript /storage/home/kangsm/runScripts/NGS_service.DrawCumulativeOncoplots.R \
        --BASE_DIR /data/wes --SEQ_FOLDER "WES_25_09" --CLIENT_ID 15
###
### preset gene variants excel file
    Rscript /storage/home/kangsm/runScripts/NGS_service.ExportPresetGeneVariantsAsExcel.R \
        --BASE_DIR /data/wes \
        --SEQ_FOLDER "WES_25_09" 
###
### data export
    /storage/home/kangsm/runScripts/NGS_service.ExportData_v2.sh \
        --baseDir /data/wes \
        --seqFolder "WES_25_09" \
        --checksum md5
###

    ID_LIST_FILE=/data/wes/WES_25_09/meta/WES_25_09_export_data.list
    EXPORT_DIR=/data/wes/WES_25_09/export_data
    DATA_DIR=/data/wes/WES_25_09

    for IDS in $(seq 1 $(cat $ID_LIST_FILE | wc -l))
    do
        TASK_ID=$(sed -n -e "${IDS} p" ${ID_LIST_FILE})
        SEQ_ID=$(echo ${TASK_ID} | awk '{print $2}')
        SAMPLE_ID=$(echo ${TASK_ID} | awk '{print $3}')
        EXPORT_TONLY=$(echo ${TASK_ID} | awk '{print $4}')
        EXPORT_TS_N=$(echo ${TASK_ID} | awk '{print $5'})
        EXPORT_ORG_N=$(echo ${TASK_ID} | awk '{print $6}')
        EXPORT_GERMLINE=$(echo ${TASK_ID} | awk '{print $7}')
        EXPORT_TONLY_GF=$(echo ${TASK_ID} | awk '{print $8}')
        TAG_TS_N="mutect2.NT.bias.filtered"

        echo -e "|--->>> ${SEQ_ID} matched normal tissue call annotated MAF copying..."
        MAF_FILE=${DATA_DIR}/${SEQ_ID}/vcf/${SEQ_ID}.${TAG_TS_N}.vep.annotated.maf
        if [ -f ${MAF_FILE} ]; then 
            cp ${MAF_FILE} ${EXPORT_DIR}/maf/${SAMPLE_ID}.matched.normal.TISSUE.maf
            sed -i -e "s/${SEQ_ID}/${SAMPLE_ID}/g" ${EXPORT_DIR}/maf/${SAMPLE_ID}.matched.normal.TISSUE.maf
        fi

        echo -e "|--->>> ${SEQ_ID} matched normal tissue call somatic variants MAF copying..."
        SOMATIC_MAF_FILE=${DATA_DIR}/${SEQ_ID}/vcf/${SEQ_ID}.${TAG_TS_N}.somatic.variants.only.maf
        if [ -f ${SOMATIC_MAF_FILE} ]; then 
            cp ${SOMATIC_MAF_FILE} ${EXPORT_DIR}/maf/${SAMPLE_ID}.matched.normal.TISSUE.somatic.variants.only.maf
            sed -i -e "s/${SEQ_ID}/${SAMPLE_ID}/g" ${EXPORT_DIR}/maf/${SAMPLE_ID}.matched.normal.TISSUE.somatic.variants.only.maf
        fi
    done
###



