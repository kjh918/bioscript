
#! /bin/bash

    BASE_DIR="/storage/home/kangsm/analysis"
    SEQ_FOLDER="gcx.rnaseq.maqc"

    RUN_SINGULARITY="singularity exec -B /storage,/data /storage/images"
    THREADS=15
    RAW_FASTQ_DIR=${BASE_DIR}/${SEQ_FOLDER}/fastq

    PE_SAMPLE_LIST=$(cat ${BASE_DIR}/${SEQ_FOLDER}/meta/pe-sample.list)
    for SEQ_ID in ${PE_SAMPLE_LIST}
    do
        FASTQ_DIR=${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/fastq
        QC_DIR=${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/qcfiles
        LOG_DIR=${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/log
        if [ ! -d ${FASTQ_DIR} ]; then mkdir -p ${FASTQ_DIR}; fi
        if [ ! -d ${LOG_DIR} ]; then mkdir -p ${LOG_DIR}; fi
        if [ ! -d ${QC_DIR} ]; then 
            mkdir -p ${QC_DIR}
            mkdir -p ${QC_DIR}/qc_res
        fi
    done

    for SEQ_ID in ${PE_SAMPLE_LIST}
    do
        QC_DIR=${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/qcfiles
        ${RUN_SINGULARITY}/fastqc-0.12.1.sif fastqc --extract \
        --threads ${THREADS} --outdir ${QC_DIR} \
        ${RAW_FASTQ_DIR}/${SEQ_ID}_*.fastq.gz
    done

    for SEQ_ID in ${PE_SAMPLE_LIST}
    do
        FASTQ_DIR=${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/fastq
        QC_DIR=${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/qcfiles
        ${RUN_SINGULARITY}/fastp-0.23.4.sif fastp \
            --in1 ${RAW_FASTQ_DIR}/${SEQ_ID}_R1.fastq.gz --in2 ${RAW_FASTQ_DIR}/${SEQ_ID}_R2.fastq.gz \
            --json ${QC_DIR}/${SEQ_ID}.fastp.json \
            --html ${FASTQ_DIR}/${SEQ_ID}.fastp.html \
            --length_required 70 \
            --qualified_quality_phred 20 \
            --thread ${THREADS} \
            --detect_adapter_for_pe \
            --trim_poly_g \
            --out1 ${FASTQ_DIR}/${SEQ_ID}.trimmed_R1.fastq.gz --out2 ${FASTQ_DIR}/${SEQ_ID}.trimmed_R2.fastq.gz
    done

    GENOME_ASSEMBLY="hg19"
    STAR_INDEX=/storage/references_and_index/${GENOME_ASSEMBLY}/rsem_index/pe100
    GTF_FILE=/storage/references_and_index/${GENOME_ASSEMBLY}/gtf/${GENOME_ASSEMBLY}.gtf
    RUN_STAR="singularity exec -B /storage,/data /storage/images/star-2.7.11.sif STAR"
    RUN_SAMTOOLS="singularity exec -B /storage,/data /storage/images/samtools-1.18.sif samtools"
    MAX_MISMATCH=10

    for SEQ_ID in ${PE_SAMPLE_LIST}
    do
        /storage/home/kangsm/runScripts/WTS_run.FastqToBAM.sh \
            --BASE_DIR "/storage/home/kangsm/analysis" \
            --SEQ_FOLDER "gcx.rnaseq.maqc" \
            --SEQ_ID ${SEQ_ID} \
            --THREADS 20 \
            --READ_LENGTH 100 \
            --STRANDED "false" \
            --GENOME_ASSEMBLY "hg19"
    done

    for SEQ_ID in ${PE_SAMPLE_LIST}
    do
        /storage/home/kangsm/runScripts/WTS_run.Quantification_RSEM.sh \
            --BASE_DIR "/storage/home/kangsm/analysis" \
            --SEQ_FOLDER "gcx.rnaseq.maqc" \
            --SEQ_ID ${SEQ_ID} \
            --THREADS 20 \
            --READ_LENGTH 100 \
            --STRANDED "false" \
            --GENOME_ASSEMBLY "hg19"
    done

    for SEQ_ID in ${PE_SAMPLE_LIST}
    do
        /storage/home/kangsm/runScripts/WTS_run.AlignedBAM_QC.sh \
            --BASE_DIR "/storage/home/kangsm/analysis" \
            --SEQ_FOLDER "gcx.rnaseq.maqc" \
            --SEQ_ID ${SEQ_ID} \
            --THREADS 20 \
            --STRANDED "false" \
            --GENOME_ASSEMBLY "hg19"
    done



    GENOME_FASTA=/storage/references_and_index/hg19/fasta/human_g1k_v37_decoy.fasta
    SEQC_GTF=/storage/references_and_index/${GENOME_ASSEMBLY}/gtf/hg19.genes.gtf
    REF_FLAT=/storage/references_and_index/${GENOME_ASSEMBLY}/extra/${GENOME_ASSEMBLY}_refFlat.txt
    rRNA_INTERVAL=/storage/references_and_index/${GENOME_ASSEMBLY}/extra/${GENOME_ASSEMBLY}_rRNA.interval_list

    GATK4="singularity exec -B /storage,/data /storage/images/gatk-4.4.0.0.sif gatk"
    RNASEQC_RUN="singularity exec -B /storage,/data /storage/images/rna-seqc-2.4.2.sif rnaseqc"
    STRAND_OPTION="NONE"

    for SEQ_ID in ${PE_SAMPLE_LIST}
    do
        QC_DIR=${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/qcfiles
        BAM_DIR=${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/bam
        SORTED_BAM=${BAM_DIR}/${SEQ_ID}.Aligned.out.Sorted.bam 

        ${RNASEQC_RUN} \
            ${SEQC_GTF} \
            ${SORTED_BAM} \
            ${QC_DIR} \
            --sample ${SEQ_ID} \
            --fasta ${GENOME_FASTA} \
            --base-mismatch 15 \
            --gene-length 500 \
            --coverage \
            --detection-threshold=3

    done

    RSEM_INDEX=/storage/references_and_index/hg19/rsem_index/pe100/hg19
    RUN_RSEM="singularity exec -B /storage,/data /storage/images/rsem-1.3.3.sif rsem-calculate-expression"

    for SEQ_ID in ${PE_SAMPLE_LIST}
    do
        QC_DIR=${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/qcfiles
        BAM_DIR=${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/bam
        QUANT_DIR=${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/quant
        if [ ! -d ${QUANT_DIR} ]; then mkdir -p ${QUANT_DIR}; fi

        ${RUN_RSEM} \
            --no-bam-output --estimate-rspd \
            --strandedness ${STRAND_OPTION} \
            --num-threads ${THREADS} \
            --alignments ${BAM_DIR}/${SEQ_ID}.Aligned.toTranscriptome.out.bam \
            ${RSEM_INDEX} \
            ${QUANT_DIR}/${SEQ_ID}  
    done

