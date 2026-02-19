WES_run.FastqToBAM.sh
#!/bin/bash

#==============================================================================#
# ScriptID : GMC-NGS-B02-A04
# Author : kangsm
# version : v1.0 
# lastest update date : 2025.12.26
#==============================================================================#

#---| Input Options |------------------------------------------------------------------------------#
    usage() {
        echo "Usage: WES_run.FastqToBAM.sh [ ARGUMENTS ]...
        [ -d | --seqFolder ]      : data processing base directory
        [ -s | --seqID ]          : sample ID (Used as sample's processing product folder name)
        [ --baseDir ]             : Base Work Dir
        [ --targetBed ]           : target BED file 
        [ --threads ]             : N-CPU threads
        [ --assembly ]            : Reference genome assembly
        [ --RGID ]                : ReadGroup ID        (default = FLOWCELL ID)
        [ --RGLB ]                : ReadGroup Library   (default = TWISTEXOME2)
        [ --RGPL ]                : ReadGroup Platform  (default = ILLUMINA)
        [ --RGCN ]                : ReadGroup Platform  (default = GCX)
        [ -h | --help ]           : Print this usage 
        "
    }
#--------------------------------------------------------------------------------------------------#

#---| Set Defualt Values |-------------------------------------------------------------------------#
    BASE_DIR=/data/wes
    THREADS=15
    RGLB="twist.exome.2.0"
    RGPL="ILLUMINA"
    RGCN="GCX"
    GENOME_ASSEMBLY="hg19"
#--------------------------------------------------------------------------------------------------#

#---| Parse Option Arguments |---------------------------------------------------------------------#
    ARGS=$(getopt -a -o d:s:h: --long seqFolder:,seqID:,baseDir:,targetBed:,threads:,assembly:,RGID:,RGLB:,RGPL:,RGCN:,help -- "$@" )
    VALID_ARGS=$?
    if [ "$VALID_ARGS" != "0" ]; then 
        usage >&2 
        exit 2
    fi

    eval set -- "$ARGS"
    while :
    do
    case "$1" in
        -d | --seqFolder )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  SEQ_FOLDER=$2 ; shift 2 ;;
            esac ;;
        -s | --seqID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  SEQ_ID=$2 ; shift 2 ;;
            esac  ;;
        --baseDir )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = /data/wes" ; shift ;;
                *)  BASE_DIR=$2 ; shift 2 ;;
            esac ;;
        --targetBed )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  TARGET_BED=$2 ; shift 2 ;;
            esac ;;
        --threads )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 15" ; shift ;;
                *)  THREADS=$2 ; shift 2 ;;
            esac ;;
        --assembly )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = hg19" ; shift ;;
                *)  GENOME_ASSEMBLY=$2 ; shift 2 ;;
            esac ;;
        --RGID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  RGID=$2 ; shift 2 ;;
            esac ;;
        --RGLB )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = TWISTEXOME2" ; shift ;;
                *)  RGLB=$2 ; shift 2 ;;
            esac ;;
        --RGPL )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = ILLUMINA" ; shift ;;
                *)  RGPL=$2 ; shift 2 ;;
            esac ;;
        --RGCN )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = GCX" ; shift ;;
                *)  RGCN=$2 ; shift 2 ;;
            esac ;;
        -h | --help )
        usage >&2 ; exit 2 ;;   
        --) shift ; break ;;
        *)  usage >&2 ; exit 2 ;;
    esac
    done
#--------------------------------------------------------------------------------------------------#

#---| FOLDERS |------------------------------------------------------------------------------------#
    DATA_DIR=${BASE_DIR}/${SEQ_FOLDER}
    ##  
    FASTQ_DIR=${DATA_DIR}/${SEQ_ID}/fastq
    BAM_DIR=${DATA_DIR}/${SEQ_ID}/bam
    QC_DIR=${DATA_DIR}/${SEQ_ID}/qcfiles
    TMP_DIR=${DATA_DIR}/${SEQ_ID}/tmp
    LOG_DIR=${DATA_DIR}/${SEQ_ID}/log
    ##
    if [ ! -d ${BAM_DIR} ]; then mkdir -p ${BAM_DIR}; fi
    if [ ! -d ${TMP_DIR} ]; then mkdir -p ${TMP_DIR}; fi
    if [ ! -d ${LOG_DIR} ]; then mkdir -p ${LOG_DIR}; fi
    if [ ! -d ${QC_DIR} ]; then 
        mkdir -p ${QC_DIR}
        mkdir -p ${QC_DIR}/qc_res
    fi
#--------------------------------------------------------------------------------------------------#

#---| LOGFILE |------------------------------------------------------------------------------------#
    LOGFILE=${LOG_DIR}/${SEQ_ID}.fastq.to.analysisReady.bam.$(date '+%Y%m%d').log
    if [ ! -f ${LOGFILE} ]; then touch ${LOGFILE}; fi
#--------------------------------------------------------------------------------------------------#

#---| SWTOOLS RUN |--------------------------------------------------------------------------------#
    RUN_SINGULARITY="singularity exec -B /storage,/data /storage/images"
    RUN_PICARD="java -XX:ParallelGCThreads=14 -Xmx16384m -jar /storage/apps/bin/picard.jar"
    GATK3="${RUN_SINGULARITY}/gatk-3.8-1.sif java -Xmx16384m -jar /usr/GenomeAnalysisTK.jar"
    GATK4="${RUN_SINGULARITY}/gatk-4.4.0.0.sif gatk"
    GATK4_JAR_16G="${RUN_SINGULARITY}/gatk-4.4.0.0.sif java -XX:ParallelGCThreads=14 -Xmx16384m -jar /gatk/gatk-package-4.4.0.0-local.jar"
#--------------------------------------------------------------------------------------------------#

#---| REFERENCES |---------------------------------------------------------------------------------#
    if [ "${GENOME_ASSEMBLY}" = "hg38" ]; then
        # BWA INDEX #
        BWA_INDEX=/storage/references_and_index/hg38/fasta/Homo_sapiens_assembly38.fasta
        # GENOME FASTA #
        GENOME_FASTA=/storage/references_and_index/hg38/fasta/Homo_sapiens_assembly38.fasta
        # REFERENCE VCF #
        # SITES
        VCF_DBSNP=/storage/references_and_index/hg38/vcf/Homo_sapiens_assembly38.dbsnp138.vcf.gz
        VCF_1KG_SNP=/storage/references_and_index/hg38/vcf/1000G_phase1.snps.high_confidence.hg38.vcf.gz
        # INDEL
        VCF_INDEL=/storage/references_and_index/hg38/vcf/Homo_sapiens_assembly38.known_indels.vcf.gz
        VCF_1KG_MILLS=/storage/references_and_index/hg38/vcf/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
    else
        # BWA INDEX #
        BWA_INDEX=/storage/references_and_index/hg19/fasta/human_g1k_v37_decoy.fasta
        # GENOME FASTA #
        GENOME_FASTA=/storage/references_and_index/hg19/fasta/human_g1k_v37_decoy.fasta
        # REFERENCE VCF #
        # SITES
        VCF_DBSNP=/storage/references_and_index/hg19/vcf/dbsnp_138.b37.vcf.gz
        VCF_1KG_SNP=/storage/references_and_index/hg19/vcf/1000G_phase1.snps.high_confidence.b37.vcf.gz
        # INDEL
        VCF_INDEL=/storage/references_and_index/hg19/vcf/Homo_sapiens_assembly19.known_indels.vcf.gz
        VCF_1KG_INDEL=/storage/references_and_index/hg19/vcf/1000G_phase1.indels.b37.vcf.gz
        VCF_1KG_MILLS=/storage/references_and_index/hg19/vcf/Mills_and_1000G_gold_standard.indels.b37.vcf.gz
    fi
    ## TARGET_BED <- From Arguments
#--------------------------------------------------------------------------------------------------#    

#==================================================================================================# 
    echo " " >> ${LOGFILE}
    echo "WES FastqToBAM log" >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " RUN DATE             : $(date '+%Y-%m-%d %H:%M:%S')" >> ${LOGFILE}
    echo " SEQ FOLDER           : ${DATA_DIR}" >> ${LOGFILE}
    echo " SEQ ID               : ${SEQ_ID}" >> ${LOGFILE}
    echo " READ GROUP ID        : ${RGID}" >> ${LOGFILE}
    echo " LIBRARY              : ${RGLB}" >> ${LOGFILE}
    echo " REFERENCE            : ${GENOME_ASSEMBLY}" >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " " >> ${LOGFILE}   
#==================================================================================================# 

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | FastqToBAM START." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | uBAM Generation START." >> ${LOGFILE}
#==================================================================================================# 

#---| uBMA (Unligned BAM) GENERATION |-------------------------------------------------------------#
    ${RUN_PICARD} FastqToSam \
        --FASTQ ${FASTQ_DIR}/${SEQ_ID}.trimmed_R1.fastq.gz \
        --FASTQ2 ${FASTQ_DIR}/${SEQ_ID}.trimmed_R2.fastq.gz \
        --OUTPUT ${BAM_DIR}/${SEQ_ID}.fastqtosam.bam \
        --READ_GROUP_NAME ${RGID} \
        --PLATFORM ${RGPL} \
        --LIBRARY_NAME ${RGLB} \
        --SAMPLE_NAME ${SEQ_ID} \
        --SEQUENCING_CENTER ${RGCN} \
        --TMP_DIR ${TMP_DIR}
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | uBAM Generation FINISHED." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | BWA aligned BAM Generation START." >> ${LOGFILE}
#==================================================================================================# 

#---| aBAM (Aligned BAM) GENERATION |--------------------------------------------------------------#
    ${RUN_SINGULARITY}/bwa-0.7.17.sif bwa mem \
        -M -t ${THREADS} \
        -R "@RG\tID:${RGID}\tPL:${RGPL}\tLB:${RGLB}\tSM:${SEQ_ID}\tCN:${RGCN}" \
        ${BWA_INDEX} \
        ${FASTQ_DIR}/${SEQ_ID}.trimmed_R1.fastq.gz \
        ${FASTQ_DIR}/${SEQ_ID}.trimmed_R2.fastq.gz > ${BAM_DIR}/${SEQ_ID}.bwa.mem.sam
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | BWA aligned BAM Generation FINISHED." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | BAM Merge Alignment START." >> ${LOGFILE}
#==================================================================================================# 

#---| BAM MERGE ALIGNMENT |------------------------------------------------------------------------#
    ${RUN_PICARD} MergeBamAlignment \
        --REFERENCE_SEQUENCE ${GENOME_FASTA} \
        --UNMAPPED_BAM ${BAM_DIR}/${SEQ_ID}.fastqtosam.bam \
        --ALIGNED_BAM ${BAM_DIR}/${SEQ_ID}.bwa.mem.sam \
        --OUTPUT ${BAM_DIR}/${SEQ_ID}.merged.bam \
        --CREATE_INDEX true \
        --ADD_MATE_CIGAR true \
        --CLIP_ADAPTERS false \
        --CLIP_OVERLAPPING_READS true \
        --INCLUDE_SECONDARY_ALIGNMENTS true \
        --MAX_INSERTIONS_OR_DELETIONS -1 \
        --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
        --ATTRIBUTES_TO_RETAIN XS \
        --TMP_DIR ${TMP_DIR}
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | BAM Merge Alignment FINISHED." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | MARK-DUPLICATES START." >> ${LOGFILE}
#==================================================================================================# 

#---| MARK-DUPLICATES |----------------------------------------------------------------------------#
    ${GATK4} MarkDuplicates \
        --java-options "-XX:ParallelGCThreads=14 -Xmx16384m" \
        --INPUT ${BAM_DIR}/${SEQ_ID}.merged.bam \
        --METRICS_FILE ${QC_DIR}/${SEQ_ID}.mark.duplicates.metrics.txt \
        --OUTPUT ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.bam \
        --CREATE_INDEX true
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | MARK-DUPLICATES FINISHED." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | LOCAL RE-ALIGNMENT START." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') |  - Re-align Target Creation Start." >> ${LOGFILE}
#==================================================================================================# 

#---| LOCAL RE-ALIGNMENT : TargetCreation |---------------------------------------------------------#
    if [ "${GENOME_ASSEMBLY}" = "hg38" ]; then
        ${GATK3} -T RealignerTargetCreator \
            -R ${GENOME_FASTA} \
            -L ${TARGET_BED} \
            -known ${VCF_INDEL} \
            -known ${VCF_1KG_MILLS} \
            -I ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.bam \
            -o ${QC_DIR}/${SEQ_ID}.realignertargetcreator.intervals
    else
        ${GATK3} -T RealignerTargetCreator \
            -R ${GENOME_FASTA} \
            -L ${TARGET_BED} \
            -known ${VCF_INDEL} \
            -known ${VCF_1KG_INDEL} \
            -known ${VCF_1KG_MILLS} \
            -I ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.bam \
            -o ${QC_DIR}/${SEQ_ID}.realignertargetcreator.intervals
    fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') |  - Re-align Target Creation Finished." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') |  - Indel Re-align Start." >> ${LOGFILE}
#==================================================================================================# 

#---| LOCAL RE-ALIGNMENT : INDEL Re-align |--------------------------------------------------------#
    if [ "${GENOME_ASSEMBLY}" = "hg38" ]; then
        ${GATK3} -T IndelRealigner \
            -R ${GENOME_FASTA} \
            -targetIntervals ${QC_DIR}/${SEQ_ID}.realignertargetcreator.intervals \
            -known ${VCF_INDEL} \
            -known ${VCF_1KG_MILLS} \
            -I ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.bam \
            -o ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.realign.bam  
    else
        ${GATK3} -T IndelRealigner \
            -R ${GENOME_FASTA} \
            -targetIntervals ${QC_DIR}/${SEQ_ID}.realignertargetcreator.intervals \
            -known ${VCF_INDEL} \
            -known ${VCF_1KG_INDEL} \
            -known ${VCF_1KG_MILLS} \
            -I ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.bam \
            -o ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.realign.bam  
    fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') |  - Indel Re-align Finished." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | BASE RECALIBRATION START." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') |  - Base Recalibration Start." >> ${LOGFILE}
#==================================================================================================# 

#---| BASE RECALIBRATION : Base Recalibration Table |----------------------------------------------#
    if [ "${GENOME_ASSEMBLY}" = "hg38" ]; then
        ${GATK4} BaseRecalibrator \
            --java-options "-XX:ParallelGCThreads=14 -Xmx16384m" \
            --input ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.realign.bam \
            --reference ${GENOME_FASTA} \
            --known-sites ${VCF_DBSNP} \
            --known-sites ${VCF_1KG_SNP} \
            --known-sites ${VCF_INDEL} \
            --known-sites ${VCF_1KG_MILLS} \
            --output ${QC_DIR}/${SEQ_ID}.recal.table.txt
    else
        ${GATK4} BaseRecalibrator \
            --java-options "-XX:ParallelGCThreads=14 -Xmx16384m" \
            --input ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.realign.bam \
            --reference ${GENOME_FASTA} \
            --known-sites ${VCF_DBSNP} \
            --known-sites ${VCF_1KG_SNP} \
            --known-sites ${VCF_INDEL} \
            --known-sites ${VCF_1KG_INDEL} \
            --known-sites ${VCF_1KG_MILLS} \
            --output ${QC_DIR}/${SEQ_ID}.recal.table.txt
    fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') |  - Base Recalibration Finished." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') |  - Apply BQRS Start." >> ${LOGFILE}
#==================================================================================================# 

#---| BASE RECALIBRATION : Apply BQRS |------------------------------------------------------------#
    ${GATK4} ApplyBQSR \
        --java-options "-XX:ParallelGCThreads=14 -Xmx16384m" \
        --input ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.realign.bam \
        --output ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.realign.recal.bam \
        --bqsr-recal-file ${QC_DIR}/${SEQ_ID}.recal.table.txt
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') |  - Apply BQRS Finished." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | BASE RECALIBRATION FINISHED." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | Rename Final-Bam to analysisReady BAM" >> ${LOGFILE}
#==================================================================================================# 

#---| CREATE ANALYSIS-READY-BAM |------------------------------------------------------------------#
    cp ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.realign.recal.bam ${BAM_DIR}/${SEQ_ID}.analysisReady.bam
    cp ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.realign.recal.bai ${BAM_DIR}/${SEQ_ID}.analysisReady.bam.bai
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo " " >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | FastqToBAM DONE." >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " " >> ${LOGFILE}
#==================================================================================================# 

