

#Rscript WES/v1/NGS_service.CreateProject.R \
Rscript create_project.R \
    --BASE_DIR "/data/wes" \
    --BATCH_ID "WES_26_01" \
    --CLIENT_ID 11 \
    --SAMPLE_INFO_FILE "WES_26_01_seqid_info.xlsx" \
    --ORDER_ID "GCX-C11-260331" \
    --DP_RD_FROM_FASTQ TRUE \
    --SEQ_TYPE "DNA" \
    --THREADS 15 \
    --WES_NGS_LIB "twist.exome.2.0" \
    --WES_RunGermlineVariant TRUE