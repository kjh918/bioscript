
#! /bin/bash

#$ -N seqc2_hg38_single
#$ -S /bin/bash
#$ -pe mpi 10
#$ -t 1-13
#$ -tc 8
#$ -e /storage/home/kangsm/analysis/cnv_seqc2/sge_log/cnv_re_run_Tonly.20240417.err
#$ -o /storage/home/kangsm/analysis/cnv_seqc2/sge_log/cnv_re_run_Tonly.20240417.out

#--------------------------------------------------------------------------------------------------#
BASE_DIR=/storage/home/kangsm/analysis/cnv_seqc2
SEQ_FOLDER=hg38
INPUT_LIST_FILE=${BASE_DIR}/${SEQ_FOLDER}/meta/sample2.list
#--------------------------------------------------------------------------------------------------#
SEQ_INFO=$(sed -n -e "${SGE_TASK_ID} p" ${INPUT_LIST_FILE})
SEQ_ID=$(echo ${SEQ_INFO} | awk '{print $1}')
SEQ_DEPTH=$(echo ${SEQ_INFO} | awk '{print $2}')
#FLOWCELL_ID=$(echo ${SEQ_INFO} | awk '{print $3}')
#READ_LENGTH=$(echo ${SEQ_INFO} | awk '{print $4}')
#--------------------------------------------------------------------------------------------------#
THREADS=20
GENOME_ASSEMBLY="hg38"
WES_LIB_KIT="sureselect.v6"
#--------------------------------------------------------------------------------------------------#
RGLB=${WES_LIB_KIT}
RUN_SCRIPT=/storage/home/kangsm/runScripts
BED_PREFIX=/storage/references_and_index/${GENOME_ASSEMBLY}/bed/${WES_LIB_KIT}/${GENOME_ASSEMBLY}_${WES_LIB_KIT}
#--------------------------------------------------------------------------------------------------#
    #ln -s ${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/vcf/${SEQ_ID}.mutect2.purecn.filtered.vcf ${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/vcf/${SEQ_ID}.mutect2.keep.germline.filtered.vcf
#--------------------------------------------------------------------------------------------------#
#     ${RUN_SCRIPT}/WES_run.CNV_PureCN.sh \
#         --seqFolder ${SEQ_FOLDER} \
#         --baseDir ${BASE_DIR} \
#         --tumorSeqID ${SEQ_ID} \
#         --runMode 'tonly' \
#         --assembly ${GENOME_ASSEMBLY} \
#         --wesLibKit ${WES_LIB_KIT} \
#         --callMinDepth 20
# #--------------------------------------------------------------------------------------------------#
#     Rscript ${RUN_SCRIPT}/WES_run.CNV_superFreq_rev.R \
#         --baseDir ${BASE_DIR} \
#         --seqFolder ${SEQ_FOLDER} \
#         --tumorSeqID ${SEQ_ID} \
#         --setID ${SEQ_ID}.Tonly \
#         --genomeAssembly ${GENOME_ASSEMBLY} \
#         --wesLibKit ${WES_LIB_KIT} \
#         --mode "exome" \
#         --threads ${THREADS} \
#         --vcfTag "mutect2.keep.germline.filtered"
#--------------------------------------------------------------------------------------------------#
    ${RUN_SCRIPT}/WES_run.CNV_CNVkit.sh \
        --tumorSeqID ${TUMOR_SEQ_ID} \
        --normalSeqID ${NORMAL_SEQ_ID} \
        --baseDir ${BASE_DIR} \
        --seqFolder ${SEQ_FOLDER} \
        --runMode 'nt' \
        --wesLibKit ${WES_LIB_KIT} \
        --assembly ${GENOME_ASSEMBLY} \
        --threads ${THREADS}

