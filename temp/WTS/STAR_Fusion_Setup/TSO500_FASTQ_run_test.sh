

    RUN_SINGULARITY="singularity exec -B /storage,/data /storage/images"
    RUN_STAR_FUSION="${RUN_SINGULARITY}/star-fusion-0.13.0.sif STAR-Fusion"
    CTAT_RESOURCE_LIB=/storage/references_and_index/hg19/star_fusion/CTAT_resource_lib
    TREADS=16

    BASE_DIR=/storage/home/kangsm/analysis/star-fusion/tso500_fastq_run

    for ID in "23Cancer01_09" "23Cancer01_10"
    do

        RES_DIR=${BASE_DIR}/${ID}
        if [ ! -d ${RES_DIR} ]; then mkdir -p ${RES_DIR}; fi

        ${RUN_STAR_FUSION} \
            --genome_lib_dir ${CTAT_RESOURCE_LIB} \
            --left_fq ${BASE_DIR}/fastq/${ID}_R1.fastq.gz \
            --right_fq ${BASE_DIR}/fastq/${ID}_R2.fastq.gz \
            --output_dir ${RES_DIR} \
            --CPU ${TREADS}

    done

#---| CTAT-slicing |
    RUN_SINGULARITY="singularity exec -B /storage,/data /storage/images"
    RUN_CTAT_SPLICING="${RUN_SINGULARITY}/ctat-splicing-0.0.2.sif /usr/local/src/CTAT-SPLICING"

    CTAT_RESOURCE_LIB=/storage/references_and_index/hg19/star_fusion/CTAT_resource_lib

    ID="23Cancer01_09"
    DATA_DIR=/storage/home/kangsm/analysis/star-fusion/tso500_fastq_run/23Cancer01_09
    RES_DIR=${DATA_DIR}/ctat_splicing
    mkdir -p ${RES_DIR}

    samtools sort -@20 --output-fmt BAM -o ${DATA_DIR}/Aligned.out.sorted.bam ${DATA_DIR}/Aligned.out.bam
    samtools index  -@20 -b ${DATA_DIR}/Aligned.out.sorted.bam

    ${RUN_CTAT_SPLICING}/STAR_to_cancer_introns.py \
        --ctat_genome_lib ${CTAT_RESOURCE_LIB} \
        --SJ_tab_file  ${DATA_DIR}/SJ.out.tab \
        --chimJ_file ${DATA_DIR}/Chimeric.out.junction \
        --vis \
        --bam_file ${DATA_DIR}/Aligned.out.sorted.bam \
        --output_prefix ${RES_DIR}/${ID} \
        --sample_name ${ID}

#--------------------------------------------------------------------------------------------------#

    RUN_SINGULARITY="singularity exec -B /storage,/data /storage/images"
    RUN_STAR="${RUN_SINGULARITY}/star-2.7.11.sif STAR"
    RUN_SAMTOOLS="${RUN_SINGULARITY}/samtools-1.18.sif samtools"

    STAR_INDEX=/storage/references_and_index/hg19/rsem_index/pe100
    GTF_FILE=/storage/references_and_index/hg19/gtf/hg19.gtf
    BAM_DIR=/storage/home/kangsm/analysis/gcx.wts.pipeline.fusionGenesDetection/TSO500_fastq_test_run/bam
    FASTQ_DIR=/storage/home/kangsm/analysis/gcx.wts.pipeline.fusionGenesDetection/TSO500_fastq_test_run/fastq

    SEQ_ID="23Cancer01_10"

    ${RUN_STAR} \
        --genomeDir ${STAR_INDEX} \
        --readFilesIn ${FASTQ_DIR}/${SEQ_ID}_R1.fastq.gz ${FASTQ_DIR}/${SEQ_ID}_R2.fastq.gz \
        --outFileNamePrefix ${BAM_DIR}/${SEQ_ID}. \
        --sjdbGTFfile ${GTF_FILE} \
        --readFilesCommand zcat \
        --quantMode TranscriptomeSAM \
        --twopassMode Basic \
        --genomeLoad NoSharedMemory \
        --runThreadN 15 \
        --outFilterMatchNminOverLread 0.33 \
        --outFilterScoreMinOverLread 0.33 \
        --readQualityScoreBase 20 \
        --outFilterMismatchNmax 10 \
        --outSAMtype BAM Unsorted \
        --outSAMunmapped Within \
        --chimSegmentMin 10 \
        --chimJunctionOverhangMin 10 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --chimScoreJunctionNonGTAG 0 \
        --chimMultimapNmax 50 \
        --peOverlapNbasesMin 10 \
        --alignSplicedMateMapLminOverLmate 0.5 \
        --outBAMcompression 0 \
        --outFilterMultimapNmax 50 \
        --chimOutType WithinBAM HardClip \
        --chimScoreDropMax 30 \
        --chimScoreSeparation 1 \
        --chimSegmentReadGapMax 3  

    ##
    ${RUN_SAMTOOLS} sort -@ 15 -o ${BAM_DIR}/${SEQ_ID}.Aligned.out.Sorted.bam ${BAM_DIR}/${SEQ_ID}.Aligned.out.bam
    ${RUN_SAMTOOLS} index -@ 15 -b ${BAM_DIR}/${SEQ_ID}.Aligned.out.Sorted.bam 

#--------------------------------------------------------------------------------------------------#

    ARRIBA_RES_DIR=/storage/home/kangsm/analysis/gcx.wts.pipeline.fusionGenesDetection/TSO500_fastq_test_run/arriba
    ARRIBA_RUN="${RUN_SINGULARITY}/arriba-2.4.0.sif /arriba_v2.4.0"

    ARRIBA_REFERENCE=/storage/references_and_index/hg19/arriba
    ARRIBA_BLACKLIST=${ARRIBA_REFERENCE}/blacklist_hg19.tsv.gz
    ARRIBA_KNOWN_FUSIONS=${ARRIBA_REFERENCE}/known_fusions_hg19.tsv.gz
    ARRIBA_PROTEIN_FEATURES=${ARRIBA_REFERENCE}/protein_domains_hg19.gff3
    ARRIBA_CYTOBAND=${ARRIBA_REFERENCE}/cytobands_hg19.tsv
    ARRIBA_GTF=${ARRIBA_REFERENCE}/hg19.gtf
    GENOME_FASTA=/storage/references_and_index/hg19/fasta/human_g1k_v37_decoy.fasta

    ${ARRIBA_RUN}/arriba \
        -x ${BAM_DIR}/${SEQ_ID}.Aligned.out.Sorted.bam \
        -o ${ARRIBA_RES_DIR}/${SEQ_ID}.arriba.fusion.genes.tsv \
        -O ${ARRIBA_RES_DIR}/${SEQ_ID}.arriba.filtered.out.fusion.genes.tsv \
        -a ${GENOME_FASTA} \
        -g ${ARRIBA_GTF} \
        -b ${ARRIBA_BLACKLIST} \
        -k ${ARRIBA_KNOWN_FUSIONS} \
        -p ${ARRIBA_PROTEIN_FEATURES} \
        -i "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y" 

