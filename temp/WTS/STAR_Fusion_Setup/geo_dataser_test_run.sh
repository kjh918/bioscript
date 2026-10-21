
#! /bin/bash

    RUN_SINGULARITY="singularity exec -B /storage,/data /storage/images"
    RUN_STAR_FUSION="${RUN_SINGULARITY}/star-fusion-0.13.0.sif STAR-Fusion"
    CTAT_RESOURCE_LIB=/storage/references_and_index/hg19/star_fusion/CTAT_resource_lib
    TREADS=24

    ## GSE230773
    FQ_DIR=/storage/home/kangsm/analysis/star-fusion/geo/GSE230773
    RES_DIR=/storage/home/kangsm/analysis/star-fusion/gse230773_test_run

    for ID in HT1080_CTRL_2D_rep1 HT1080_CTRL_2D_rep2 HT1080_CTRL_2D_rep3 HT1080_CTRL_2D_rep4 HT1080_FUS-DDIT3_2D_rep1 HT1080_FUS-DDIT3_2D_rep2 HT1080_FUS-DDIT3_2D_rep3
    do
        RESULT_DIR=${RES_DIR}/${ID}
        mkdir -p ${RESULT_DIR}

        ${RUN_STAR_FUSION} \
            --genome_lib_dir ${CTAT_RESOURCE_LIB} \
            --left_fq ${FQ_DIR}/${ID}_R1.fastq.gz \
            --right_fq ${FQ_DIR}/${ID}_R2.fastq.gz \
            --output_dir ${RESULT_DIR} \
            --CPU ${TREADS}

    done

    ## GSE243279
    FQ_DIR=/storage/home/kangsm/analysis/star-fusion/geo/GSE243279
    RES_DIR=/storage/home/kangsm/analysis/star-fusion/gse243279_test_run
    mkdir -p ${RES_DIR}

    for ID in TC32_CTRL_rep1 TC32_CTRL_rep2 TC32_CTRL_rep3 TC32_ETS1_High_rep1 TC32_ETS1_High_rep2 TC32_ETS1_High_rep3 TC32_ETS1_Low_rep1 TC32_ETS1_Low_rep2 TC32_ETS1_Low_rep3
    do
        RESULT_DIR=${RES_DIR}/${ID}
        mkdir -p ${RESULT_DIR}

        ${RUN_STAR_FUSION} \
            --genome_lib_dir ${CTAT_RESOURCE_LIB} \
            --left_fq ${FQ_DIR}/${ID}_R1.fastq.gz \
            --right_fq ${FQ_DIR}/${ID}_R2.fastq.gz \
            --output_dir ${RESULT_DIR} \
            --CPU ${TREADS}

    done

    

















  


