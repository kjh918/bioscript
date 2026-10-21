
    INPUT_WGS_BAM=
    BED_FILE=
    OUTPUT_BED_WES_BAM=
    OUTPUT_BED_WES_SORTED_BAM=
    BED_WES_FASTQ1=
    BED_WES_FASTQ2=

# 1. create BED intersect BAM 
    intersectBed -abam ${INPUT_WGS_BAM} -b ${BED_FILE} > ${OUTPUT_BED_WES_BAM}

# 2. sorting by name
    samtools sort -n -@ 15 ${OUTPUT_BED_WES_BAM} > ${OUTPUT_BED_WES_SORTED_BAM}

# 3. Bam to Fastq 
    bedtools bamtofastq -i ${OUTPUT_BED_WES_SORTED_BAM} -fq ${BED_WES_FASTQ1} -fq2 ${BED_WES_FASTQ2}

# 4. Compression
    bgzip -@ 10 -c ${BED_WES_FASTQ1} > ${BED_WES_FASTQ1}.gz
    bgzip -@ 10 -c ${BED_WES_FASTQ2} > ${BED_WES_FASTQ2}.gz






    ID="sample_110"
    WGS_BAM=/data/wgs/hg19/KNUCH/sample_110/bam/sample_110.recaled.bam
    BED_FILE=/storage/references_and_index/hg19/bed/twist.exome.2.0/hg19_twist.exome.2.0.probe.bed
    BED_BAM_DIR=/storage/home/kangsm/analysis/gcx.MRD/VARIANT_SELECTION_METHOD_DEVELOPMENT/VALID05_KNUCH/wgs_bam_sampling

    intersectBed -abam ${WGS_BAM} -b ${BED_FILE} > ${BED_BAM_DIR}/sample_110.twist.exome.2.0.bam

    # parallel processing : bedtools
    ID_LIST=$(cat /storage/home/kangsm/analysis/gcx.MRD/VARIANT_SELECTION_METHOD_DEVELOPMENT/VALID05_KNUCH/meta/wgs_bam_sample.all.list)
    parallel -j 18 -k "intersectBed -wa -a /data/wgs/hg19/KNUCH/{}/bam/{}.recaled.bam -b /storage/references_and_index/hg19/bed/twist.exome.2.0/hg19_twist.exome.2.0.probe.bed > /storage/home/kangsm/analysis/gcx.MRD/VARIANT_SELECTION_METHOD_DEVELOPMENT/VALID05_KNUCH/wgs_bam_sampling/{}_bedtools.twist.exome.2.0.bam" ::: $ID_LIST

    # parallel processing : samtools
    ID_LIST=$(cat /storage/home/kangsm/analysis/gcx.MRD/VARIANT_SELECTION_METHOD_DEVELOPMENT/VALID05_KNUCH/meta/wgs_bam_sample.all.list)

    parallel -j 18 -k "samtools view -@ 4 -b -L /storage/references_and_index/hg19/bed/twist.exome.2.0/hg19_twist.exome.2.0.probe.bed /data/wgs/hg19/KNUCH/{}/bam/{}.recaled.bam > /storage/home/kangsm/analysis/gcx.MRD/VARIANT_SELECTION_METHOD_DEVELOPMENT/VALID05_KNUCH/wgs_bam_sampling/{}_twist.exome.2.0.bam" ::: $ID_LIST

    P14_DX
    WGS_DATA_DIR=/storage/home/kangsm/analysis/gcx.MRD/VARIANT_SELECTION_METHOD_DEVELOPMENT/VALID05_KNUCH/wgs_bam_sampling
    

    SAMPLE_ID=$(cat /storage/home/kangsm/analysis/gcx.MRD/VARIANT_SELECTION_METHOD_DEVELOPMENT/VALID05_KNUCH/meta/sample_list)

    parallel -j 18 -k "samtools sort -n -@4 ${WGS_DATA_DIR}/{}.bam > ${WGS_DATA_DIR}/{}.name.sorted.bam" ::: $SAMPLE_ID
    
    parallel -j 18 -k "bedtools bamtofastq -i ${WGS_DATA_DIR}/{}.name.sorted.bam -fq ${WGS_DATA_DIR}/{}_R1.fastq -fq2 ${WGS_DATA_DIR}/{}_R2.fastq" ::: $SAMPLE_ID

    parallel -j 18 -k "bgzip -@ 4 -c ${WGS_DATA_DIR}/{}_R1.fastq > ${WGS_DATA_DIR}/{}_R1.fastq.gz" ::: $SAMPLE_ID
    parallel -j 18 -k "bgzip -@ 4 -c ${WGS_DATA_DIR}/{}_R2.fastq > ${WGS_DATA_DIR}/{}_R2.fastq.gz" ::: $SAMPLE_ID

    
    BASE_DIR=/storage/home/kangsm/analysis/gcx.MRD/VARIANT_SELECTION_METHOD_DEVELOPMENT/VALID05_KNUCH
    
    for ID in P18 P30
    do
        mkdir -p ${BASE_DIR}/${ID}/fastq
        cp ${WGS_DATA_DIR}/${ID}_TUMOR_R1.fastq.gz ${WGS_DATA_DIR}/${ID}_TUMOR_R2.fastq.gz ${WGS_DATA_DIR}/${ID}_NORMAL_R1.fastq.gz ${WGS_DATA_DIR}/${ID}_NORMAL_R2.fastq.gz ${BASE_DIR}/${ID}/fastq/
    done

    for ID in P14 P18 P30
    do
        mkdir -p ${BASE_DIR}/${ID}/${ID}_NORMAL/fastq
        cp ${BASE_DIR}/${ID}/fastq/${ID}_NORMAL_R1.fastq.gz ${BASE_DIR}/${ID}/${ID}_NORMAL/fastq/${ID}_NORMAL.trimmed_R1.fastq.gz
        cp ${BASE_DIR}/${ID}/fastq/${ID}_NORMAL_R2.fastq.gz ${BASE_DIR}/${ID}/${ID}_NORMAL/fastq/${ID}_NORMAL.trimmed_R2.fastq.gz

        mkdir -p ${BASE_DIR}/${ID}/${ID}_DX/fastq
        cp ${BASE_DIR}/wgs_bam_sampling/${ID}_DX_R1.fastq.gz ${BASE_DIR}/${ID}/${ID}_DX/fastq/${ID}_DX.trimmed_R1.fastq.gz
        cp ${BASE_DIR}/wgs_bam_sampling/${ID}_DX_R2.fastq.gz ${BASE_DIR}/${ID}/${ID}_DX/fastq/${ID}_DX.trimmed_R2.fastq.gz
        
        mkdir -p ${BASE_DIR}/${ID}/${ID}_RELAPSE/fastq
        cp ${BASE_DIR}/wgs_bam_sampling/${ID}_RELAPSE_R1.fastq.gz ${BASE_DIR}/${ID}/${ID}_RELAPSE/fastq/${ID}_RELAPSE.trimmed_R1.fastq.gz
        cp ${BASE_DIR}/wgs_bam_sampling/${ID}_RELAPSE_R2.fastq.gz ${BASE_DIR}/${ID}/${ID}_RELAPSE/fastq/${ID}_RELAPSE.trimmed_R2.fastq.gz
        
        mkdir -p ${BASE_DIR}/${ID}/${ID}_CFDNA01/fastq
        cp ${BASE_DIR}/wgs_bam_sampling/${ID}_CFDNA01_R1.fastq.gz ${BASE_DIR}/${ID}/${ID}_CFDNA01/fastq/${ID}_CFDNA01.trimmed_R1.fastq.gz
        cp ${BASE_DIR}/wgs_bam_sampling/${ID}_CFDNA01_R2.fastq.gz ${BASE_DIR}/${ID}/${ID}_CFDNA01/fastq/${ID}_CFDNA01.trimmed_R2.fastq.gz
        
        mkdir -p ${BASE_DIR}/${ID}/${ID}_CFDNA02/fastq
        cp ${BASE_DIR}/wgs_bam_sampling/${ID}_CFDNA02_R1.fastq.gz ${BASE_DIR}/${ID}/${ID}_CFDNA02/fastq/${ID}_CFDNA02.trimmed_R1.fastq.gz
        cp ${BASE_DIR}/wgs_bam_sampling/${ID}_CFDNA02_R2.fastq.gz ${BASE_DIR}/${ID}/${ID}_CFDNA02/fastq/${ID}_CFDNA02.trimmed_R2.fastq.gz

        mkdir -p ${BASE_DIR}/${ID}/${ID}_TUMOR_P100/fastq
        cp ${BASE_DIR}/${ID}/fastq/${ID}_TUMOR_R1.fastq.gz ${BASE_DIR}/${ID}/${ID}_TUMOR_P100/fastq/${ID}_TUMOR_P100.trimmed_R1.fastq.gz
        cp ${BASE_DIR}/${ID}/fastq/${ID}_TUMOR_R2.fastq.gz ${BASE_DIR}/${ID}/${ID}_TUMOR_P100/fastq/${ID}_TUMOR_P100.trimmed_R2.fastq.gz
    
        mkdir -p ${BASE_DIR}/${ID}/${ID}_TUMOR_P80/fastq
        cp ${BASE_DIR}/${ID}/fastq/${ID}_TUMOR_P80_R1.fastq.gz ${BASE_DIR}/${ID}/${ID}_TUMOR_P80/fastq/${ID}_TUMOR_P80.trimmed_R1.fastq.gz
        cp ${BASE_DIR}/${ID}/fastq/${ID}_TUMOR_P80_R2.fastq.gz ${BASE_DIR}/${ID}/${ID}_TUMOR_P80/fastq/${ID}_TUMOR_P80.trimmed_R2.fastq.gz

        mkdir -p ${BASE_DIR}/${ID}/${ID}_TUMOR_P60/fastq
        cp ${BASE_DIR}/${ID}/fastq/${ID}_TUMOR_P60_R1.fastq.gz ${BASE_DIR}/${ID}/${ID}_TUMOR_P60/fastq/${ID}_TUMOR_P60.trimmed_R1.fastq.gz
        cp ${BASE_DIR}/${ID}/fastq/${ID}_TUMOR_P60_R2.fastq.gz ${BASE_DIR}/${ID}/${ID}_TUMOR_P60/fastq/${ID}_TUMOR_P60.trimmed_R2.fastq.gz

        mkdir -p ${BASE_DIR}/${ID}/${ID}_TUMOR_P40/fastq
        cp ${BASE_DIR}/${ID}/fastq/${ID}_TUMOR_P40_R1.fastq.gz ${BASE_DIR}/${ID}/${ID}_TUMOR_P40/fastq/${ID}_TUMOR_P40.trimmed_R1.fastq.gz
        cp ${BASE_DIR}/${ID}/fastq/${ID}_TUMOR_P40_R2.fastq.gz ${BASE_DIR}/${ID}/${ID}_TUMOR_P40/fastq/${ID}_TUMOR_P40.trimmed_R2.fastq.gz
    done














