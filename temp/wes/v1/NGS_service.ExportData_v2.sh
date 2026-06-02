
#! /bin/bash

# SCRIPT VERSION : v2.0
# DATE           : 2025-03-06
# AUTHOR         : kangsm@gencurix.com
# CHANGES        : 'SEQ-ID' 표기를 'Sample Name'으로 변경하는 부분 추가

#---| Input Options |------------------------------------------------------------------------------#
    usage() {
        echo "Usage: GCX-WES-TumorOnlyVariantCalling.sh [ options ]
        [ --seqFolder ]      : data processing base directory
        [ --baseDir ]        : work base folder. default = /data/wes
        [ --checksum ]       : integrity check. 'md5', 'sha256' or 'none'. default = 'none'
        [ -h | --help ]      : Print this usage 
        "
    }
#--------------------------------------------------------------------------------------------------#

#---| Set Defualt Values |-------------------------------------------------------------------------#
    BASE_DIR=/data/wes
    CHECKSUM="none"
#--------------------------------------------------------------------------------------------------#

#---| Parse Option Arguments |---------------------------------------------------------------------#
    ARGS=$(getopt -a -o h: --long seqFolder:,baseDir:,checksum:,help -- "$@" )
    VALID_ARGS=$?
    if [ "$VALID_ARGS" != "0" ]; then 
        usage >&2 
        exit 2
    fi

    eval set -- "$ARGS"
    while :
    do
    case "$1" in
        --seqFolder )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  SEQ_FOLDER=$2 ; shift 2 ;;
            esac ;;
        --baseDir )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = /data/wes" ; shift ;;
                *)  BASE_DIR=$2 ; shift 2 ;;
            esac ;;
        --checksum )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = sha256" ; shift ;;
                *)  CHECKSUM=$2 ; shift 2 ;;
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
    RAW_FASTQ=${BASE_DIR}/fastq
    LOG_DIR=${DATA_DIR}/meta
    EXPORT_DIR=${DATA_DIR}/export_data
    #--------------------------------------------------------------------------#
    if [ ! -d ${EXPORT_DIR} ]; then
        mkdir -p ${EXPORT_DIR}
        mkdir -p ${EXPORT_DIR}/fastq
        mkdir -p ${EXPORT_DIR}/bam
        mkdir -p ${EXPORT_DIR}/vcf
        mkdir -p ${EXPORT_DIR}/maf
        mkdir -p ${EXPORT_DIR}/report
    fi
#--------------------------------------------------------------------------------------------------#

#---| LOGFILE |------------------------------------------------------------------------------------#
    LOGFILE=${LOG_DIR}/${SEQ_FOLDER}.export.data.$(date '+%Y%m%d').log
    if [ ! -f ${LOGFILE} ]; then touch ${LOGFILE}; fi
#--------------------------------------------------------------------------------------------------#

#---| RUN CHECKSUM CMD |---------------------------------------------------------------------------#
    if [ "${CHECKSUM}" = "md5" ]; then
        RUN_CHECKSUM="md5sum"
        FILE_EXT="md5sum"
    elif [ "${CHECKSUM}" = "sha256" ]; then
        RUN_CHECKSUM="shasum -a 256"
        FILE_EXT="sha256"
    fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo " " >> ${LOGFILE}
    echo "DATA EXPORT log" >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " RUN DATE             : $(date '+%Y-%m-%d %H:%M:%S')" >> ${LOGFILE}
    echo " SEQ FOLDER           : ${DATA_DIR}" >> ${LOGFILE}
    echo " FILE INTEGRITY CHECK : ${CHECKSUM}" >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " " >> ${LOGFILE}   
#==================================================================================================# 

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | EXPORT ARCHIVE GENERATION START." >> ${LOGFILE}
#==================================================================================================# 

#---| FILE COPY |----------------------------------------------------------------------------------#
    ID_LIST_FILE=${DATA_DIR}/meta/${SEQ_FOLDER}_export_data.list
    TASK_IDS=$(seq 1 $(cat $ID_LIST_FILE | wc -l))    
    #--------------------------------------------------------------------------#
    echo " "
    for IDS in $TASK_IDS
    do
        TASK_ID=$(sed -n -e "${IDS} p" ${ID_LIST_FILE})
        SEQ_ID=$(echo ${TASK_ID} | awk '{print $2}')
        SAMPLE_ID=$(echo ${TASK_ID} | awk '{print $3}')
        EXPORT_TONLY=$(echo ${TASK_ID} | awk '{print $4}')
        EXPORT_TS_N=$(echo ${TASK_ID} | awk '{print $5'})
        EXPORT_ORG_N=$(echo ${TASK_ID} | awk '{print $6}')
        EXPORT_GERMLINE=$(echo ${TASK_ID} | awk '{print $7}')
        EXPORT_TONLY_GF=$(echo ${TASK_ID} | awk '{print $8}')
        
        # start file copy ------------------------------------------------------
        echo -e "---[] ${SEQ_ID} : ${SAMPLE_ID} export ]-----------------------"

        #---/ fastq /-----------------------------------------------------------
        # echo -e " \n"
        if [ -f ${RAW_FASTQ}/${SEQ_ID}_R1.fastq.gz ]; then
            echo -e "|--->>> ${SEQ_ID}_R1.fastq.gz copying..."
            cp ${RAW_FASTQ}/${SEQ_ID}_R1.fastq.gz ${EXPORT_DIR}/fastq/${SAMPLE_ID}_R1.fastq.gz
            echo -e "|--->>> ${SEQ_ID}_R2.fastq.gz copying..."
            cp ${RAW_FASTQ}/${SEQ_ID}_R2.fastq.gz ${EXPORT_DIR}/fastq/${SAMPLE_ID}_R2.fastq.gz
        else
            echo -e "|--->>> No FASTQ of ${SEQ_ID} found."
            echo -e "${SEQ_ID} - ${SAMPLE_ID} : NO FASTQ FILES. NOT EXPORTED." >> ${LOGFILE}   
        fi

        #---/ bam /-------------------------------------------------------------
        # echo -e " \n"
        if [ -f ${DATA_DIR}/${SEQ_ID}/bam/${SEQ_ID}.analysisReady.bam ]; then
            echo -e "|--->>> ${SEQ_ID}.analysisReady.bam copying..."
            cp ${DATA_DIR}/${SEQ_ID}/bam/${SEQ_ID}.analysisReady.bam ${EXPORT_DIR}/bam/${SAMPLE_ID}.bam
            echo -e "|--->>> ${SEQ_ID}.analysisReady.bam.bai copying..."
            cp ${DATA_DIR}/${SEQ_ID}/bam/${SEQ_ID}.analysisReady.bam.bai ${EXPORT_DIR}/bam/${SAMPLE_ID}.bam.bai
        else
            echo -e "|--->>> No BAM of ${SEQ_ID} found."
            echo -e "${SEQ_ID} - ${SAMPLE_ID} : NO BAM FILE. NOT EXPORTED." >> ${LOGFILE}   
        fi

        #---/ vcf and maf /-----------------------------------------------------
        # echo -e " \n"
        if [ "${EXPORT_TONLY}" = "YES" ]; then
            TAG_TONLY="mutect2.bias.filtered"
            
            echo -e "|--->>> ${SEQ_ID} tumor only call VCF copying..."
            VCF_FILE=${DATA_DIR}/${SEQ_ID}/vcf/${SEQ_ID}.${TAG_TONLY}.vcf
            if [ -f ${VCF_FILE} ]; then 
                cp ${VCF_FILE} ${EXPORT_DIR}/vcf/${SAMPLE_ID}.tumor.only.vcf
                sed -i -e "s/${SEQ_ID}/${SAMPLE_ID}/g" ${EXPORT_DIR}/vcf/${SAMPLE_ID}.tumor.only.vcf
            fi

            echo -e "|--->>> ${SEQ_ID} tumor only call annotated MAF copying..."
            MAF_FILE=${DATA_DIR}/${SEQ_ID}/vcf/${SEQ_ID}.${TAG_TONLY}.vep.annotated.maf
            if [ -f ${MAF_FILE} ]; then 
                cp ${MAF_FILE} ${EXPORT_DIR}/maf/${SAMPLE_ID}.tumor.only.maf
                sed -i -e "s/${SEQ_ID}/${SAMPLE_ID}/g" ${EXPORT_DIR}/maf/${SAMPLE_ID}.tumor.only.maf
            fi

            echo -e "|--->>> ${SEQ_ID} tumor only call somatic variants MAF copying..."
            SOMATIC_MAF_FILE=${DATA_DIR}/${SEQ_ID}/vcf/${SEQ_ID}.${TAG_TONLY}.somatic.variants.only.maf
            if [ -f ${SOMATIC_MAF_FILE} ]; then 
                cp ${SOMATIC_MAF_FILE} ${EXPORT_DIR}/maf/${SAMPLE_ID}.tumor.only.somatic.variants.only.maf
                sed -i -e "s/${SEQ_ID}/${SAMPLE_ID}/g" ${EXPORT_DIR}/maf/${SAMPLE_ID}.tumor.only.somatic.variants.only.maf
            fi
        fi

        # echo -e " \n"
        if [ "${EXPORT_TS_N}" = "YES" ]; then
            TAG_TS_N="mutect2.NT.bias.filtered"

            echo -e "|--->>> ${SEQ_ID} matched normal tissue call VCF copying..."
            VCF_FILE=${DATA_DIR}/${SEQ_ID}/vcf/${SEQ_ID}.${TAG_TS_N}.vcf
            if [ -f ${VCF_FILE} ]; then 
                cp ${VCF_FILE} ${EXPORT_DIR}/vcf/${SAMPLE_ID}.matched.normal.TISSUE.vcf
                sed -i -e "s/${SEQ_ID}/${SAMPLE_ID}/g" ${EXPORT_DIR}/vcf/${SAMPLE_ID}.matched.normal.TISSUE.vcf
            fi

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
        fi

        # echo -e " \n"
        if [ "${EXPORT_ORG_N}" = "YES" ]; then
            TAG_ORG_N="mutect2.ORG.NT.bias.filtered"

            echo -e "|--->>> ${SEQ_ID} matched normal organoid call VCF copying..."
            VCF_FILE=${DATA_DIR}/${SEQ_ID}/vcf/${SEQ_ID}.${TAG_ORG_N}.vcf
            if [ -f ${VCF_FILE} ]; then cp 
                ${VCF_FILE} ${EXPORT_DIR}/vcf/${SAMPLE_ID}.matched.normal.ORGANOID.vcf
                sed -i -e "s/${SEQ_ID}/${SAMPLE_ID}/g" ${EXPORT_DIR}/vcf/${SAMPLE_ID}.matched.normal.ORGANOID.vcf
            fi

            echo -e "|--->>> ${SEQ_ID} matched normal organoid call annotated MAF copying..."
            MAF_FILE=${DATA_DIR}/${SEQ_ID}/vcf/${SEQ_ID}.${TAG_ORG_N}.vep.annotated.maf
            if [ -f ${MAF_FILE} ]; then 
                cp ${MAF_FILE} ${EXPORT_DIR}/maf/${SAMPLE_ID}.matched.normal.ORGANOID.maf
                sed -i -e "s/${SEQ_ID}/${SAMPLE_ID}/g" ${EXPORT_DIR}/maf/${SAMPLE_ID}.matched.normal.ORGANOID.maf
            fi

            echo -e "|--->>> ${SEQ_ID} matched normal orgnoaid call somatic variants MAF copying..."
            SOMATIC_MAF_FILE=${DATA_DIR}/${SEQ_ID}/vcf/${SEQ_ID}.${TAG_ORG_N}.somatic.variants.only.maf
            if [ -f ${SOMATIC_MAF_FILE} ]; then 
                cp ${SOMATIC_MAF_FILE} ${EXPORT_DIR}/maf/${SAMPLE_ID}.matched.normal.ORGANOID.somatic.variants.only.maf
                sed -i -e "s/${SEQ_ID}/${SAMPLE_ID}/g" ${EXPORT_DIR}/maf/${SAMPLE_ID}.matched.normal.ORGANOID.somatic.variants.only.maf
            fi
        fi

        if [ "${EXPORT_GERMLINE}" = "YES" ]; then
            TAG_GERMLINE="deepvariant.germline.variant.filtered"

            echo -e "|--->>> ${SEQ_ID} deepvariant call germline variants VCF copying..."
            VCF_FILE=${DATA_DIR}/${SEQ_ID}/vcf/${SEQ_ID}.${TAG_GERMLINE}.vcf
            if [ -f ${VCF_FILE} ]; then 
                cp ${VCF_FILE} ${EXPORT_DIR}/vcf/${SAMPLE_ID}.germline.vcf
                sed -i -e "s/${SEQ_ID}/${SAMPLE_ID}/g" ${EXPORT_DIR}/vcf/${SAMPLE_ID}.germline.vcf
            fi

            echo -e "|--->>> ${SEQ_ID} deepvariant call germline variants annotated MAF copying..."
            MAF_FILE=${DATA_DIR}/${SEQ_ID}/vcf/${SEQ_ID}.${TAG_GERMLINE}.vep.annotated.maf
            if [ -f ${MAF_FILE} ]; then 
                cp ${MAF_FILE} ${EXPORT_DIR}/maf/${SAMPLE_ID}.germline.maf
                sed -i -e "s/${SEQ_ID}/${SAMPLE_ID}/g" ${EXPORT_DIR}/maf/${SAMPLE_ID}.germline.maf
            fi

            echo -e "|--->>> ${SEQ_ID} deepvariant call germline SYNONYMOUS variants MAF copying..."
            SYN_MAF_FILE=${DATA_DIR}/${SEQ_ID}/vcf/${SEQ_ID}.${TAG_GERMLINE}.synonymous.maf
            if [ -f ${SYN_MAF_FILE} ]; then 
                cp ${SYN_MAF_FILE} ${EXPORT_DIR}/maf/${SAMPLE_ID}.germline.synonymous.variants.only.maf
                sed -i -e "s/${SEQ_ID}/${SAMPLE_ID}/g" ${EXPORT_DIR}/maf/${SAMPLE_ID}.germline.synonymous.variants.only.maf
            fi

            echo -e "|--->>> ${SEQ_ID} deepvariant call germline NON-SYNONYMOUS variants MAF copying..."
            NONSYN_MAF_FILE=${DATA_DIR}/${SEQ_ID}/vcf/${SEQ_ID}.${TAG_GERMLINE}.non_synonymous.maf
            if [ -f ${NONSYN_MAF_FILE} ]; then 
                cp ${NONSYN_MAF_FILE} ${EXPORT_DIR}/maf/${SAMPLE_ID}.germline.nonsynonymous.variants.only.maf
                sed -i -e "s/${SEQ_ID}/${SAMPLE_ID}/g" ${EXPORT_DIR}/maf/${SAMPLE_ID}.germline.nonsynonymous.variants.only.maf
            fi
        fi
       
        #======================================================================================#
            echo "$(date '+%Y-%m-%d %H:%M:%S') | ${SEQ_ID} : ${SAMPLE_ID} : ${SAMPLE_TYPE} : Exported = FASTQ, BAM, VCF, MAF" >> ${LOGFILE}
        #======================================================================================#
        
    done

    cp ${BASE_DIR}/${SEQ_FOLDER}/analysis/export/*.png ${EXPORT_DIR}/report/
    cp ${BASE_DIR}/${SEQ_FOLDER}/analysis/export/*.xlsx ${EXPORT_DIR}/report/

#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---" >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | File Copy Finished." >> ${LOGFILE}
#==================================================================================================# 

#---| INTEGRITY CHECK SAMPLES |--------------------------------------------------------------------#
    if [ "${CHECKSUM}" = "md5" ] || [ "${CHECKSUM}" = "sha256" ]; then
        cd ${EXPORT_DIR}/fastq
        ${RUN_CHECKSUM} *.fastq.gz > FASTQ.${FILE_EXT}
        cd ${EXPORT_DIR}/bam
        ${RUN_CHECKSUM} *.bam > BAM.${FILE_EXT}
        ${RUN_CHECKSUM} *.bai >> BAM.${FILE_EXT}
        cd ${EXPORT_DIR}/vcf
        ${RUN_CHECKSUM} *.vcf > VCF.${FILE_EXT}
        cd ${EXPORT_DIR}/maf
        ${RUN_CHECKSUM} *.maf > MAF.${FILE_EXT}
    fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | File Integrity Checking Finished." >> ${LOGFILE}
#==================================================================================================# 

#==================================================================================================# 
    echo " " >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | EXPORT ARCHIVE GENERATION FINISHED." >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " " >> ${LOGFILE}
#==================================================================================================#
