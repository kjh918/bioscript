
#-------------------------------------------#
#- WES Variant Calling Analysis Run Script -#
#-                                         -#
#- Date : 2023-09-07                       -#
#- Ver  : v1.0                             -#
#-------------------------------------------#

#!/bin/bash

DX_DIR=/BiO/DATA/WES/GCX230727001-7_OrganoidScience

#---/ DIRECTORIES
FASTQ_RAW_DIR=${DX_DIR}/fastq_raw
FASTQ_MERGED_DIR=${DX_DIR}/fastq_merged
FASTQ_TRIMMED_DIR=${DX_DIR}/fastq_trimmed

if [ ! -d $FASTQ_MERGED_DIR ]; then mkdir -p ${FASTQ_MERGED_DIR}; fi
if [ ! -d $FASTQ_TRIMMED_DIR ]; then mkdir -p ${FASTQ_TRIMMED_DIR}; fi

#---/ RESULT FOLDERS Depend on Genome Version
if [ "$REF_VER" == "hg19" ]; then 
    mkdir -p ${DX_DIR}/hg19
    RES_DIR=${DX_DIR}/hg19
else
    mkdir -p ${DX_DIR}/hg38
    RES_DIR=${DX_DIR}/hg38
fi

mkdir -p ${RES_DIR}/qc_reports && QC_REPORT_DIR=${RES_DIR}/qc_reports
mkdir -p $QC_REPORT_DIR/fastqc_1st
mkdir -p $QC_REPORT_DIR/fastqc_2nd
mkdir -p ${RES_DIR}/tmp && TMP_DIR=${RES_DIR}/tmp
mkdir -p ${RES_DIR}/bam_tmp && BAM_TMP_DIR=${RES_DIR}/bam_tmp
mkdir -p ${RES_DIR}/bam_ready && BAM_FINAL_DIR=${RES_DIR}/bam_ready
mkdir -p ${RES_DIR}/vcf_out && VCF_OUT_DIR=${RES_DIR}/vcf_out

#---/

#---/ SWTOOLS PARAMTERS
SWTOOLS=/BiO/SWTools
BWA_INDEX=/BiO/Reference_and_Index/hg19/bwa_index/gcx.hg19
REF_GENOME=/BiO/Reference_and_Index/hg19/gcx.hg19.fasta
RGID="GCX"
RGLB="TwistExome2"
RGPL="ILLUMINA"
REF_VCF=

SAMPLE_LIST=$(cat /BiO/DATA/WES/GCX230727001-7_OrganoidScience/sample.list)


#---/ MERGE FASTQ FILES INTO 1 FASTQ-FILE
    if[ "FASTQ_MERGE" == "true" ]; then
        for SeqID in '17-100X' '17-200X' '17-400X' '3-100X' '3-200X' '3-400X' 'tissue-200X'
        do
            for LANE in 1 2 3 4
            do 
                R1_Files=$(ls ${FASTQ_RAW_DIR} | grep "${SeqID}_S[0-9]_L0[0-9]${LANE}_R1")
                R2_Files=$(ls ${FASTQ_RAW_DIR} | grep "${SeqID}_S[0-9]_L0[0-9]${LANE}_R2")
                cat ${FASTQ_RAW_DIR}/${R1_Files} >> ${FASTQ_MERGED_DIR}/${SeqID}_R1.fastq.gz
                cat ${FASTQ_RAW_DIR}/${R2_Files} >> ${FASTQ_MERGED_DIR}/${SeqID}_R2.fastq.gz
            done
        done
    else
        FASTQ_MERGED_DIR=${FASTQ_RAW_DIR}
    fi
#---/

SAMPLE_LIST=$(cat /BiO/DATA/WES/GCX230727001-7_OrganoidScience/sample.list)

#---/ FASTQ-QC
    ${SWTOOLS}/fastqc \
        -o ${QC_REPORT_DIR}/fastqc_1st \
        --noextract \
        -t 50 \
        ${FASTQ_MERGED_DIR}/*.fastq.gz
#---/


#---/ FASTQ-TRIMMING & TRIMMING-QC
    parallel -j 6 -k "${SWTOOLS}/fastp \
        --in1 ${FASTQ_MERGED_DIR}/{}_R1.fastq.gz \
        --in2 ${FASTQ_MERGED_DIR}/{}_R2.fastq.gz \
        --html ${QC_REPORT_DIR}/{}.fastq-qc.html \
        --json ${QC_REPORT_DIR}/{}.fastq-qc.json \
        --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        --out1 ${FASTQ_TRIMMED_DIR}/{}.trimmed_R1.fastq.gz \
        --out2 ${FASTQ_TRIMMED_DIR}/{}.trimmed_R2.fastq.gz" ::: ${SAMPLE_LIST}
#---/ 

#---/ FASTQ-QC 2
    ${SWTOOLS}/fastqc \
        -o ${QC_REPORT_DIR}/fastqc_2nd \
        --noextract \
        -t 50 \
        ${FASTQ_TRIMMED_DIR}/*.fastq.gz
#---/


#---/ BWA ALIGNMENT
    parallel -j 6 -k "${SWTOOLS}/bwa mem -aM \
        -R '@RG\tID:${RGID}\tLB:${RGLB}\tSM:{}\tPL:${RGPL}' \
        -t 20 \
        -k 70 \
        ${BWA_INDEX} \
        ${FASTQ_TRIMMED_DIR}/{}.trimmed_R1.fastq.gz \
        ${FASTQ_TRIMMED_DIR}/{}.trimmed_R2.fastq.gz |
        ${SWTOOLS}/samtools view -Shb -o ${BAM_TMP_DIR}/{}.bam - " ::: ${SAMPLE_LIST}
#---/

#---/ ALIGNMENT POST-PROCESSING
  #---/ SAM CONVERT TO BAM
    parallel -j 6 -k "${SWTOOLS}/sambamba view -S \
        -t ${THREADS} \
        -f bam \
        -o ${BAM_TMP_DIR}/{}.bam \
        ${BAM_TMP_DIR}/{}.sam" ::: ${SAMPLE_LIST}

  #---/ BAM SORT
    parallel -j 6 -k "${SWTOOLS}/sambamba sort -t ${THREADS} \
        -o ${BAM_TMP_DIR}/{}.sorted.bam \
        ${BAM_TMP_DIR}/{}.bam" ::: ${SAMPLE_LIST}

  #---/ MARKDUPLICATE (REMOVE DUPLICATES)
    #---/ 전임자 방법 
    parallel -j 6 -k "${SWTOOLS}/sambamba markdup -r -p \
        -t ${THREADS} \
        --overflow-list-size 600000 \
        --tmpdir ${TMP_DIR} \
        ${BAM_TMP_DIR}/{}.sorted.bam \
        ${BAM_FINAL_DIR}/{}.sorted.dedup.bam" ::: ${SAMPLE_LIST}

    #---/ GATK 방법
    parallel -j 7 -k "java -XX:+UseParallelGC -XX:ParallelGCThreads=20 -Xmx60g -jar ${SWTOOLS}/gatk-local.jar MarkDuplicates \
        -I ${BAM_TMP_DIR}/{}.sorted.bam \
        -M ${QC_REPORT_DIR}/{}.marked.dup.metrics.txt \
        --REMOVE_DUPLICATES=true \
        -O ${BAM_TMP_DIR}/{}.sorted.dedup.method2.bam" ::: ${SAMPLE_LIST} 


    #---/ COLLECT INSERT-SIZE METRICS
    parallel -j 6 -k "java -XX:+UseParallelGC -XX:ParallelGCThreads=20 -Xmx20g -jar ${SWTOOLS}/gatk-local.jar CollectInsertSizeMetrics \
        -I ${BAM_FINAL_DIR}/{}.sorted.dedup.bam \
        -O ${BAM_FINAL_DIR}/{}.sorted.dedup.insert.metrics.txt \
        -H ${BAM_FINAL_DIR}/{}.sorted.dedup.insert.size.histogram.pdf" ::: ${SAMPLE_LIST}


  #---/ BASE RE-CALIBRATOR
    parallel -j 7 -k "java -XX:+UseParallelGC -XX:ParallelGCThreads=20 -Xmx60g -jar ${SWTOOLS}/gatk-local.jar BaseRecalibrator \
        -I ${BAM_FINAL_DIR}/{}.sorted.dedup.bam \
        -R /BiO/Reference_and_Index/hg19/gcx.hg19.fasta \
        --known-sites /BiO/Reference_and_Index/hg19/vcf/1000G_omni2.5.hg19.sites.vcf.gz \
        --known-sites /BiO/Reference_and_Index/hg19/vcf/hapmap_3.3.hg19.sites.vcf.gz \
        --known-sites /BiO/Reference_and_Index/hg19/vcf/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz \
        --known-sites /BiO/Reference_and_Index/hg19/vcf/dbsnp_138.hg19.vcf.gz \
        -O ${BAM_FINAL_DIR}/{}.recal.data.table" ::: ${SAMPLE_LIST}

  #---/ APPLY BQSR
    parallel -j 7 -k "java -XX:+UseParallelGC -XX:ParallelGCThreads=20 -Xmx60g -jar ${SWTOOLS}/gatk-local.jar ApplyBQSR \
        -I ${BAM_FINAL_DIR}/{}.sorted.dedup.bam \
        -O ${BAM_FINAL_DIR}/{}.sorted.dedup.recal.bam \
        -bqsr ${BAM_FINAL_DIR}/{}.recal.data.table" ::: ${SAMPLE_LIST}

#---/ 

#---/ SAMTOOLS STAT
    parallel -j 7 -k ${SWTOOLS}/samtools stat \
        ${BAM_FINAL_DIR}/{}.sorted.dedup.recal.bam > ${QC_REPORT_DIR}/{}.final.bam.stats.txt ::: ${SAMPLE_LIST}
     


#---/ VARIANT CALLING
    parallel -j 7 -k ${SWTOOLS}/gatk Mutect2 \
        -I ${BAM_FINAL_DIR}/{}.sorted.dedup.recal.bam \
        -R ${REF_GENOME} \
        --germline-resource /BiO/Reference_and_Index/hg19/vcf/af-only-gnomad.HG19.vcf \
        -O ${VCF_OUT_DIR}/{}.mutect2.single.vcf ::: ${SAMPLE_LIST}

    parallel -j 7 -k "bgzip -c ${VCF_OUT_DIR}/{}.mutect2.single.vcf > ${VCF_OUT_DIR}/{}.mutect2.single.vcf.gz" ::: ${SAMPLE_LIST}
    parallel -j 7 -k "tabix -p vcf ${VCF_OUT_DIR}/{}.mutect2.single.vcf.gz" ::: ${SAMPLE_LIST}


