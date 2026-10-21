

#! /bin/bash

#$ -N HMFTOOLS
#$ -S /bin/bash
#$ -pe mpi 24
#$ -t 1-1
#$ -tc 1
#$ -e /data/wgs/seqc2_wgs/meta/analysis_log/hmftools.run.test.err
#$ -o /data/wgs/seqc2_wgs/meta/analysis_log/hmftools.run.test.out

/storage/home/kangsm/runScripts/WGS/WGS_M05-1_RUN-HMFTOOLS-ANALYSIS.sh \
    --BASE_DIR "/data/wgs" \
    --BATCH_ID "seqc2_wgs" \
    --TUMOR_SEQ_ID "HCC1395" \
    --NORMAL_SEQ_ID "HCC1395BL" \
    --THREADS 24 \
    --GENOME_ASSEMBLY "hg38"    
