
#! /bin/bash

# script version v0.1

#---| Input Options |------------------------------------------------------------------------------#
    usage() {
        echo "Usage: WES_run.CNV_CNVkit.sh [ ARGUMENTS ]
          --tumorSeqID      : input tumor bam file (with absolute path).
          --normalSeqID     : input matched normal bam file (with absolute path, required if runMode = 'paired').
          --baseDir         : data processing base directory
          --seqFolder       : sample ID (used as result folder name).
          --runMode         : analysis mode. 'tonly' or 'nt'. default = 'tonly'
          --wesLibKit       : target region bed. 'sureselect','sureselect_utr','twist'. default = 'twist'
          --assembly        : reference genome version. 'hg19' or 'hg38'.
          --threads         : n-cpu
          --nType           : type of normal sample
          --superfreqSetID  : type of normal sample
          -h | --help       : Print this usage 
        "
    }
#--------------------------------------------------------------------------------------------------#

#---| Set Defualt Values |-------------------------------------------------------------------------#
    BASE_DIR=/data/wes
    WES_LIB_KIT="twist.exome.2.0"
    GENOME_ASSEMBLY="hg19"
    RUN_MODE="tonly"
    THREADS=15
    NORMAL_TYPE="TS"
#--------------------------------------------------------------------------------------------------#

#---| Parse Option Arguments |---------------------------------------------------------------------#

    ARGS=$(getopt -a -o h: --long tumorSeqID:,normalSeqID:,baseDir:,seqFolder:,runMode:,wesLibKit:,assembly:,threads:,nType:,superfreqSetID:,help: -- "$@" )
    VALID_ARGS=$?
    if [ "$VALID_ARGS" != "0" ]; then 
        usage >&2 
        exit 2
    fi

    eval set -- "$ARGS"
    while :
    do
    case "$1" in
        --tumorSeqID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  TUMOR_SEQ_ID=$2 ; shift 2 ;;
            esac ;;
        --normalSeqID )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = NULL" ; shift ;;
                *)  NORMAL_SEQ_ID=$2 ; shift 2 ;;
            esac  ;;
        --baseDir )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = /data/wes" ; shift ;;
                *)  BASE_DIR=$2 ; shift 2 ;;
            esac ;;
        --seqFolder )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  SEQ_FOLDER=$2 ; shift 2 ;;
            esac ;;
        --runMode )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = tonly" ; shift ;;
                *)  RUN_MODE=$2 ; shift 2 ;;
            esac ;;
        --wesLibKit )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = twist.exome.2.0" ; shift ;;
                *)  WES_LIB_KIT=$2 ; shift 2 ;;
            esac ;;
        --assembly )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = hg19" ; shift ;;
                *)  GENOME_ASSEMBLY=$2 ; shift 2 ;;
            esac ;;
        --threads )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 15" ; shift ;;
                *)  THREADS=$2 ; shift 2 ;;
            esac ;;
        --nType )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = TS" ; shift ;;
                *)  NORMAL_TYPE=$2 ; shift 2 ;;
            esac ;;
        --superfreqSetID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  SUPERFREQ_SET_ID=$2 ; shift 2 ;;
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
    TUMOR_BAM_DIR=${DATA_DIR}/${TUMOR_SEQ_ID}/bam
    NORMAL_BAM_DIR=${DATA_DIR}/${NORMAL_SEQ_ID}/bam
    LOG_DIR=${DATA_DIR}/${TUMOR_SEQ_ID}/log
    #--------------------------------------------------------------------------#  
    CNV_RES_DIR=${DATA_DIR}/${TUMOR_SEQ_ID}/cnv
    if [ ! -d ${CNV_RES_DIR} ]; then mkdir -p ${CNV_RES_DIR}; fi
    #--------------------------------------------------------------------------#
    CNVKIT_RES_DIR=${CNV_RES_DIR}/cnvkit
    if [ ! -d ${CNVKIT_RES_DIR} ]; then mkdir -p ${CNVKIT_RES_DIR}; fi
    #--------------------------------------------------------------------------#
    if [ ! -d ${LOG_DIR} ]; then mkdir -p ${LOG_DIR}; fi
#--------------------------------------------------------------------------------------------------#

#---| LOGFILE |------------------------------------------------------------------------------------#
    LOGFILE=${LOG_DIR}/${TUMOR_SEQ_ID}.CNV.cnvkit.$(date '+%Y%m%d').log
    if [ ! -f ${LOGFILE} ]; then touch ${LOGFILE}; fi
#--------------------------------------------------------------------------------------------------#

#---| SWTOOLS RUN |--------------------------------------------------------------------------------#   
    RUN_SINGULARITY="singularity exec -B /storage,/data /storage/images"
    RUN_CNVKIT="${RUN_SINGULARITY}/cnvkit-0.9.10.sif cnvkit.py"
    RUN_SCRIPT=/storage/home/kangsm/runScripts
#--------------------------------------------------------------------------------------------------#

#---| REFERENCES |---------------------------------------------------------------------------------#
    if [ "${GENOME_ASSEMBLY}" = "hg38" ]; then
        GENOME_FASTA=/storage/references_and_index/hg38/fasta/Homo_sapiens_assembly38.fasta
    else
        GENOME_FASTA=/storage/references_and_index/hg19/fasta/human_g1k_v37_decoy.fasta
    fi
    #--------------------------------------------------------------------------#
    CNVKIT_RESOURCES=/storage/references_and_index/cnv/cnvkit
    CNVKIT_ACCESS_BED=${CNVKIT_RESOURCES}/${GENOME_ASSEMBLY}_cnvkit.access.bed
    REF_FLAT=${CNVKIT_RESOURCES}/${GENOME_ASSEMBLY}.refFlat.txt
    WES_BAIT_BED=/storage/references_and_index/${GENOME_ASSEMBLY}/bed/${WES_LIB_KIT}/${GENOME_ASSEMBLY}_${WES_LIB_KIT}.probe.bed
    CNVKIT_PON_DIR=/data/pon/pon_cnv.cnvkit/${GENOME_ASSEMBLY}_${WES_LIB_KIT}
#--------------------------------------------------------------------------------------------------#

#---| RUN MODE TAG |-------------------------------------------------------------------------------#
    if [ "${RUN_MODE}" = "nt" ]; then
        NORMAL_SAMPPLE_TAG=${NORMAL_SEQ_ID}
    else
        NORMAL_SAMPPLE_TAG="No Normal Sample (Tumor-Only Mode)"
    fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo " " >> ${LOGFILE}
    echo "WES CNV-CNVkit log" >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " RUN DATE             : $(date '+%Y-%m-%d %H:%M:%S')" >> ${LOGFILE}
    echo " SEQ FOLDER           : ${DATA_DIR}" >> ${LOGFILE}
    echo " TUMOR_SEQ_ID         : ${TUMOR_SEQ_ID}" >> ${LOGFILE}
    echo " NORMAL_SEQ_ID        : ${NORMAL_SAMPPLE_TAG}" >> ${LOGFILE}
    echo " GENOME ASSEMBLY      : ${GENOME_ASSEMBLY}" >> ${LOGFILE}
    echo " RUN MODE             : ${RUN_MODE}" >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " " >> ${LOGFILE}
#==================================================================================================# 

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | CNV-CNVKIT START." >> ${LOGFILE} 
#==================================================================================================# 

#---| NORMAL TYPE |--------------------------------------------------------------------------------#
    if [ "${RUN_MODE}" = "nt" ]; then
        if [ "${NORMAL_TYPE}" = "ORG" ]; then
            MATCH_TAG="ORG.NT"
        else
            MATCH_TAG="NT"
        fi
    fi
#--------------------------------------------------------------------------------------------------#

#---| CNVKIT VCF |---------------------------------------------------------------------------------# 
    if [ "${RUN_MODE}" = "nt" ]; then
        #==========================================================================================# 
             echo "$(date '+%Y-%m-%d %H:%M:%S') | [ CNVkit VCF PREPARE ] Start." >> ${LOGFILE}
        #==========================================================================================#
        ##
        KEEP_GERMLINE_VCF=${DATA_DIR}/${TUMOR_SEQ_ID}/vcf/${TUMOR_SEQ_ID}.mutect2.${MATCH_TAG}.keep.germline.filtered.vcf
        ##
        if [ -f ${KEEP_GERMLINE_VCF} ]; then
            CNVKIT_VCF_TAG=${DATA_DIR}/${TUMOR_SEQ_ID}/vcf/${TUMOR_SEQ_ID}.${MATCH_TAG}.cnvkit
            #------------------------------------------------------------------#
            cp ${KEEP_GERMLINE_VCF} ${CNVKIT_VCF_TAG}.vcf
            #------------------------------------------------------------------#
            INSERT_LINE=$(grep -n "##INFO" ${CNVKIT_VCF_TAG}.vcf | head -1 | awk -F ":" '{print $1}')
            sed -i -e "${INSERT_LINE} i\##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic event\">" ${CNVKIT_VCF_TAG}.vcf
            sed -i -e "${INSERT_LINE} i\##PEDIGREE=<Derived=${TUMOR_SEQ_ID},Original=${NORMAL_SEQ_ID}>" ${CNVKIT_VCF_TAG}.vcf
            #------------------------------------------------------------------#
            sed -i -e "s/PASS\tAS_FilterStatus/PASS\tSOMATIC;AS_FilterStatus/" ${CNVKIT_VCF_TAG}.vcf
            #------------------------------------------------------------------#
            grep ^# ${CNVKIT_VCF_TAG}.vcf > ${CNVKIT_VCF_TAG}.pf.vcf
            bcftools query -i "FMT/DP>25" -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\t%FORMAT" ${CNVKIT_VCF_TAG}.vcf >> ${CNVKIT_VCF_TAG}.pf.vcf
            bcftools view -v snps ${CNVKIT_VCF_TAG}.pf.vcf > ${CNVKIT_VCF_TAG}.ready.vcf

            #======================================================================================#
                echo "$(date '+%Y-%m-%d %H:%M:%S') | [ CNVkit VCF PREPARE ] Finished." >> ${LOGFILE}
            #======================================================================================#
        else
            #======================================================================================#
                echo "$(date '+%Y-%m-%d %H:%M:%S') | No PureCN VCF. CNVkit run Without VCF." >> ${LOGFILE}
                echo "$(date '+%Y-%m-%d %H:%M:%S') | [ CNVkit VCF PREPARE ] Finished." >> ${LOGFILE}
            #======================================================================================#
        fi
    fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | [ LOAD Pre-Analysis DATA ] Start." >> ${LOGFILE}
#==================================================================================================# 

#---| PLOIDY,PURITY,SEX DATA |---------------------------------------------------------------------# 
    if [ "${RUN_MODE}" = "nt" ]; then
        ASCAT_RES_FILE=${DATA_DIR}/${TUMOR_SEQ_ID}/cnv/ascat/${TUMOR_SEQ_ID}.${MATCH_TAG}.ASCAT.summary.txt
        if [ ! -f ${ASCAT_RES_FILE} ]; then
            #======================================================================================#
                echo "$(date '+%Y-%m-%d %H:%M:%S') | NO ASCAT DATA !" >> ${LOGFILE}
                echo "$(date '+%Y-%m-%d %H:%M:%S') | CNVkit RUN STOP." >> ${LOGFILE}
                echo ">>>>> NO ASCAT DATA ! STOPPED !"
            #======================================================================================#
            exit 0 
        else
            ASCAT_DATA=$(cat ${ASCAT_RES_FILE})
            PURITY=$(echo $ASCAT_DATA | awk '{print $2}')
            PLOIDY=$(echo $ASCAT_DATA | awk '{print $4}')
            PLOIDY_INT=$(echo $ASCAT_DATA | awk '{print $5}')
        fi
        
        SUPERFREQ_RES_FILE=${DATA_DIR}/${TUMOR_SEQ_ID}/cnv/superfreq/${TUMOR_SEQ_ID}.${SUPERFREQ_SET_ID}.superFreq.summary.txt
        if [ ! -f ${SUPERFREQ_RES_FILE} ]; then
            #==========================================================================================#
                echo "$(date '+%Y-%m-%d %H:%M:%S') | NO superFreq DATA !" >> ${LOGFILE}
                echo "$(date '+%Y-%m-%d %H:%M:%S') | CNVkit RUN STOP." >> ${LOGFILE}
                echo ">>>>> NO SUPERFREQ DATA!  STOPPED !"
            #==========================================================================================#
            exit 2 
        else
            SUPERFREQ_DATA=$(cat ${SUPERFREQ_RES_FILE})
            SEX=$(echo $SUPERFREQ_DATA | awk '{print $4}')
        fi
    else
        PURECN_RES_FILE=${DATA_DIR}/${TUMOR_SEQ_ID}/cnv/purecn/${TUMOR_SEQ_ID}.Tonly.PureCN.summary.txt
        if [ ! -f ${PURECN_RES_FILE} ]; then
            SUPERFREQ_RES_FILE=${DATA_DIR}/${TUMOR_SEQ_ID}/cnv/superfreq/${TUMOR_SEQ_ID}.${SUPERFREQ_SET_ID}.superFreq.summary.txt
            if [ ! -f ${SUPERFREQ_RES_FILE} ]; then
                #======================================================================================#
                    echo "$(date '+%Y-%m-%d %H:%M:%S') | NO PureCN or superFreq DATA !" >> ${LOGFILE}
                    echo "$(date '+%Y-%m-%d %H:%M:%S') | CNVkit RUN STOP." >> ${LOGFILE}
                    echo ">>>>> NO PureCN or superFreq DATA ! STOPPED !"
                #======================================================================================#
                exit 0 
            else
                SUPERFREQ_DATA=$(cat ${SUPERFREQ_RES_FILE})
                PURITY=$(echo $SUPERFREQ_DATA | awk '{print $5}')
                PLOIDY=$(echo $SUPERFREQ_DATA | awk '{print $2}')
                PLOIDY_INT=$(echo $SUPERFREQ_DATA | awk '{print $3}')
            fi
        else
            PURECN_DATA=$(cat ${PURECN_RES_FILE})
            PURITY=$(echo $PURECN_DATA | awk '{print $2}')
            PLOIDY=$(echo $PURECN_DATA | awk '{print $3}')
            PLOIDY_INT=$(echo $PURECN_DATA | awk '{print $4}')
        fi
        SUPERFREQ_RES_FILE=${DATA_DIR}/${TUMOR_SEQ_ID}/cnv/superfreq/${TUMOR_SEQ_ID}.${SUPERFREQ_SET_ID}.superFreq.summary.txt
        if [ ! -f ${SUPERFREQ_RES_FILE} ]; then
            #==========================================================================================#
                echo "$(date '+%Y-%m-%d %H:%M:%S') | NO superFreq DATA !" >> ${LOGFILE}
                echo "$(date '+%Y-%m-%d %H:%M:%S') | CNVkit RUN STOP." >> ${LOGFILE}
                echo ">>>>> NO superFreq DATA!  STOPPED !"
            #==========================================================================================#
            exit 2 
        else
            SUPERFREQ_DATA=$(cat ${SUPERFREQ_RES_FILE})
            SEX=$(echo $SUPERFREQ_DATA | awk '{print $4}')
        fi
    fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================#
    echo "$(date '+%Y-%m-%d %H:%M:%S') | [ LOAD Pre-Analysis DATA ] Finished." >> ${LOGFILE}
    echo ">>>>> ALL NECESSARY PRE-ANALYSIS DATA FOUND."
#==================================================================================================#
    

#---| AUTOBIN |------------------------------------------------------------------------------------#
    if [ "${RUN_MODE}" = "nt" ]; then

        #==========================================================================================# 
            echo "$(date '+%Y-%m-%d %H:%M:%S') | RUN MODE = NT-Paired Mode." >> ${LOGFILE}
            echo "$(date '+%Y-%m-%d %H:%M:%S') | [ AUTOBIN ] Start." >> ${LOGFILE}
        #==========================================================================================# 

        ${RUN_CNVKIT} autobin \
            --method hybrid \
            --fasta ${GENOME_FASTA} \
            --access ${CNVKIT_ACCESS_BED} \
            --targets ${WES_BAIT_BED} \
            --annotate ${REF_FLAT} \
            --target-output-bed ${CNVKIT_RES_DIR}/${GENOME_ASSEMBLY}_${WES_LIB_KIT}.${MATCH_TAG}.target.bed \
            --antitarget-output-bed ${CNVKIT_RES_DIR}/${GENOME_ASSEMBLY}_${WES_LIB_KIT}.${MATCH_TAG}.antitarget.bed \
            ${NORMAL_BAM_DIR}/${NORMAL_SEQ_ID}.analysisReady.bam 
        #----------------------------------------------------------------------#
        TARGET_BED=${CNVKIT_RES_DIR}/${GENOME_ASSEMBLY}_${WES_LIB_KIT}.${MATCH_TAG}.target.bed
        ANTITARGET_BED=${CNVKIT_RES_DIR}/${GENOME_ASSEMBLY}_${WES_LIB_KIT}.${MATCH_TAG}.antitarget.bed

        #==========================================================================================# 
            echo "$(date '+%Y-%m-%d %H:%M:%S') | [ AUTOBIN ] Finished." >> ${LOGFILE}
        #==========================================================================================# 
    else
        #==========================================================================================# 
            echo "$(date '+%Y-%m-%d %H:%M:%S') | RUN MODE = T-only Mode." >> ${LOGFILE}
            echo "$(date '+%Y-%m-%d %H:%M:%S') | [ AUTOBIN ] Omitted." >> ${LOGFILE}
            echo "$(date '+%Y-%m-%d %H:%M:%S') | Pooled Normal will be used as Reference." >> ${LOGFILE}
        #==========================================================================================#

        PON_TARGET_BED=${CNVKIT_PON_DIR}/targets/${GENOME_ASSEMBLY}_${WES_LIB_KIT}.pon.target.bed
        PON_ANTITARGET_BED=${CNVKIT_PON_DIR}/targets/${GENOME_ASSEMBLY}_${WES_LIB_KIT}.pon.antitarget.bed

        #==========================================================================================# 
            echo "$(date '+%Y-%m-%d %H:%M:%S') | PoN Target & AntiTarget Assigned for Default T-only Mode Analysis." >> ${LOGFILE}
        #==========================================================================================# 
    fi
#--------------------------------------------------------------------------------------------------#

#---| REFERENCCE |---------------------------------------------------------------------------------#
    if [ "${RUN_MODE}" = "nt" ]; then
        
        #==========================================================================================# 
            echo "$(date '+%Y-%m-%d %H:%M:%S') | [ REFERENCE ] Matched Normal Sample Coverage Analysis Start." >> ${LOGFILE}
        #==========================================================================================#

        ${RUN_CNVKIT} coverage \
            --fasta ${GENOME_FASTA} \
            --min-mapq 20 \
            --output ${CNVKIT_RES_DIR}/${NORMAL_SEQ_ID}.${MATCH_TAG}.target.coverage.cnn \
            --processes ${THREADS} \
            ${NORMAL_BAM_DIR}/${NORMAL_SEQ_ID}.analysisReady.bam \
            ${TARGET_BED}
        #----------------------------------------------------------------------#
        ${RUN_CNVKIT} coverage \
            --fasta ${GENOME_FASTA} \
            --min-mapq 20 \
            --output ${CNVKIT_RES_DIR}/${NORMAL_SEQ_ID}.${MATCH_TAG}.antitarget.coverage.cnn \
            --processes ${THREADS} \
            ${NORMAL_BAM_DIR}/${NORMAL_SEQ_ID}.analysisReady.bam \
            ${ANTITARGET_BED}
        #----------------------------------------------------------------------#

        #==========================================================================================# 
            echo "$(date '+%Y-%m-%d %H:%M:%S') | [ REFERENCE ] Matched Normal Sample Coverage Analysis Finished." >> ${LOGFILE}
            echo "$(date '+%Y-%m-%d %H:%M:%S') | [ REFERENCE ] Matched Normal Reference Generate Start." >> ${LOGFILE}
        #==========================================================================================#

        ${RUN_CNVKIT} reference \
            --fasta ${GENOME_FASTA} \
            --output ${CNVKIT_RES_DIR}/matched.normal.${MATCH_TAG}.reference.cnn \
            ${CNVKIT_RES_DIR}/${NORMAL_SEQ_ID}.${MATCH_TAG}.*.coverage.cnn
        #----------------------------------------------------------------------#
        REF_CNN_NT=${CNVKIT_RES_DIR}/matched.normal.${MATCH_TAG}.reference.cnn

        #==========================================================================================# 
            echo "$(date '+%Y-%m-%d %H:%M:%S') | [ REFERENCE ] Matched Normal Reference Generate Finished." >> ${LOGFILE}
        #==========================================================================================#
    else
        #==========================================================================================# 
            echo "$(date '+%Y-%m-%d %H:%M:%S') | [ REFERENCE ] RUN MODE = T-only" >> ${LOGFILE}
            echo "$(date '+%Y-%m-%d %H:%M:%S') | [ REFERENCE ] PoN Reference will be used." >> ${LOGFILE}
        #==========================================================================================#   

        REF_CNN_PON=${CNVKIT_PON_DIR}/${GENOME_ASSEMBLY}_${WES_LIB_KIT}.pon.reference.cnn
    
        #==========================================================================================# 
            echo "$(date '+%Y-%m-%d %H:%M:%S') | PoN Reference Assigned for Default T-only Mode Analysis." >> ${LOGFILE}
        #==========================================================================================# 
    fi
    
#--------------------------------------------------------------------------------------------------#

#---| TUMOR BAM COVERAGE |-------------------------------------------------------------------------#
    if [ "${RUN_MODE}" = "nt" ]; then

        #==========================================================================================# 
            echo "$(date '+%Y-%m-%d %H:%M:%S') | [ COVERAGE ] Tumor Sample Coverage Analysis with NT-mode Start." >> ${LOGFILE}
        #==========================================================================================#
        
        ${RUN_CNVKIT} coverage \
            --fasta ${GENOME_FASTA} \
            --min-mapq 20 \
            --output ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.${MATCH_TAG}.target.coverage.cnn \
            --processes ${THREADS} \
            ${TUMOR_BAM_DIR}/${TUMOR_SEQ_ID}.analysisReady.bam \
            ${TARGET_BED}
        #--------------------------------------------------------------------------#
        ${RUN_CNVKIT} coverage \
            --fasta ${GENOME_FASTA} \
            --min-mapq 20 \
            --output ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.${MATCH_TAG}.antitarget.coverage.cnn \
            --processes ${THREADS} \
            ${TUMOR_BAM_DIR}/${TUMOR_SEQ_ID}.analysisReady.bam \
            ${ANTITARGET_BED}

        #==========================================================================================# 
            echo "$(date '+%Y-%m-%d %H:%M:%S') | [ COVERAGE ] Tumor Sample Coverage Analysis with NT-mode Finished." >> ${LOGFILE}
        #==========================================================================================#
    else
        #==========================================================================================# 
            echo "$(date '+%Y-%m-%d %H:%M:%S') | [ COVERAGE ] Tumor Sample Coverage Analysis with T-only mode Start." >> ${LOGFILE}
        #==========================================================================================#

        #--------------------------------------------------------------------------#
        ${RUN_CNVKIT} coverage \
            --fasta ${GENOME_FASTA} \
            --min-mapq 20 \
            --output ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.Tonly.target.coverage.cnn \
            --processes ${THREADS} \
            ${TUMOR_BAM_DIR}/${TUMOR_SEQ_ID}.analysisReady.bam \
            ${PON_TARGET_BED}
        #--------------------------------------------------------------------------#
        ${RUN_CNVKIT} coverage \
            --fasta ${GENOME_FASTA} \
            --min-mapq 20 \
            --output ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.Tonly.antitarget.coverage.cnn \
            --processes ${THREADS} \
            ${TUMOR_BAM_DIR}/${TUMOR_SEQ_ID}.analysisReady.bam \
            ${PON_ANTITARGET_BED}    

        #==========================================================================================# 
            echo "$(date '+%Y-%m-%d %H:%M:%S') | [ COVERAGE ] Tumor Sample Coverage Analysis with T-only mode Finished." >> ${LOGFILE}
        #==========================================================================================#
    fi
#--------------------------------------------------------------------------------------------------#


#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | [ FIX TUMOR COVERAGE ] Start." >> ${LOGFILE}
#==================================================================================================# 

#---| ADJUST TUMOR COVERAGE |----------------------------------------------------------------------#
    if [ "${RUN_MODE}" = "nt" ]; then    
        ${RUN_CNVKIT} fix \
            --output ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.${MATCH_TAG}.cnr \
            --sample-id ${TUMOR_SEQ_ID} \
            ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.${MATCH_TAG}.target.coverage.cnn \
            ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.${MATCH_TAG}.antitarget.coverage.cnn \
            ${REF_CNN_NT}
    else
        ${RUN_CNVKIT} fix \
            --output ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.Tonly.cnr \
            --sample-id ${TUMOR_SEQ_ID} \
            ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.Tonly.target.coverage.cnn \
            ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.Tonly.antitarget.coverage.cnn \
            ${REF_CNN_PON}
    fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | [ FIX TUMOR COVERAGE ] Finished." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | [ SEGMENTATION ] Start." >> ${LOGFILE}
#==================================================================================================# 

#---| SEGMENTATION |-------------------------------------------------------------------------------#
    if [ "${RUN_MODE}" = "nt" ]; then    
        if [ -f ${CNVKIT_VCF} ]; then
            ${RUN_CNVKIT} segment \
                --output ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.${MATCH_TAG}.cns \
                --dataframe ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.${MATCH_TAG}.cbs.txt \
                --method cbs \
                --smooth-cbs \
                --drop-low-coverage \
                --drop-outliers 4 \
                --processes ${THREADS} \
                --vcf ${CNVKIT_VCF_TAG}.ready.vcf \
                --sample-id ${TUMOR_SEQ_ID} \
                --normal-id ${NORMAL_SEQ_ID} \
                ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.${MATCH_TAG}.cnr
        else
            ${RUN_CNVKIT} segment \
                --output ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.${MATCH_TAG}.cns \
                --dataframe ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.${MATCH_TAG}.cbs.txt \
                --method cbs \
                --smooth-cbs \
                --drop-low-coverage \
                --drop-outliers 4 \
                --processes ${THREADS} \
                --sample-id ${TUMOR_SEQ_ID} \
                --normal-id ${NORMAL_SEQ_ID} \
                ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.${MATCH_TAG}.cnr
        fi
    else
        ${RUN_CNVKIT} segment \
            --output ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.Tonly.cns \
            --dataframe ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.Tonly.cbs.txt \
            --method cbs \
            --drop-outliers 5 \
            --processes ${THREADS} \
            --sample-id ${TUMOR_SEQ_ID} \
            ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.Tonly.cnr
    fi
#--------------------------------------------------------------------------------------------------#

#---| SEGMETRICS |---------------------------------------------------------------------------------#
    if [ "${RUN_MODE}" = "nt" ]; then    
        ${RUN_CNVKIT} segmetrics \
            --segments ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.${MATCH_TAG}.cns \
            --output ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.${MATCH_TAG}.segmetrics.txt \
            --drop-low-coverage \
            --mean --median --mode --mse --ci --smooth-bootstrap \
            ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.${MATCH_TAG}.cnr
    else
        ${RUN_CNVKIT} segmetrics \
            --segments ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.Tonly.cns \
            --output ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.Tonly.segmetrics.txt \
            --drop-low-coverage \
            --mean --median --mode --mse --ci --smooth-bootstrap \
            ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.Tonly.cnr
    fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | [ SEGMENTATION ] Finished." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | [ CNV CALL ] Start." >> ${LOGFILE}
#==================================================================================================# 

#---| CNV CALL |-----------------------------------------------------------------------------------#
    if [ "${RUN_MODE}" = "nt" ]; then       
        if [ -f ${CNVKIT_VCF} ]; then
            ${RUN_CNVKIT} call \
                --output ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.${MATCH_TAG}.call.cns \
                --sample-id ${TUMOR_SEQ_ID} \
                --normal-id ${NORMAL_SEQ_ID} \
                --vcf ${CNVKIT_VCF_TAG}.ready.vcf \
                --purity ${PURITY} \
                --sample-sex ${SEX} \
                ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.${MATCH_TAG}.cns
        else
            ${RUN_CNVKIT} call \
                --output ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.${MATCH_TAG}.call.cns \
                --sample-id ${TUMOR_SEQ_ID} \
                --normal-id ${NORMAL_SEQ_ID} \
                --purity ${PURITY} \
                --sample-sex ${SEX} \
                ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.${MATCH_TAG}.cns
        fi
    else
        ${RUN_CNVKIT} call \
            --output ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.Tonly.call.cns \
            --sample-id ${TUMOR_SEQ_ID} \
            --purity ${PURITY} \
            --sample-sex ${SEX} \
            ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.Tonly.cns
    fi
#--------------------------------------------------------------------------------------------------#

#---| USED PARAMTERS FOR SUMMARY AND VISUALIZATION |-----------------------------------------------#
    if [ "${RUN_MODE}" = "nt" ]; then       
        echo -e "ID\tPLOIDY\tPLOIDY_INT\tPURITY\tSEX" > ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.${MATCH_TAG}.call.params
        echo -e "${TUMOR_SEQ_ID}\t${PLOIDY}\t${PLOIDY_INT}\t${PURITY}\t${SEX}" >>${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.${MATCH_TAG}.call.params
    else
        echo -e "ID\tPLOIDY\tPLOIDY_INT\tPURITY\tSEX" > ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.Tonly.call.params
        echo -e "${TUMOR_SEQ_ID}\t${PLOIDY}\t${PLOIDY_INT}\t${PURITY}\t${SEX}" >>${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.Tonly.call.params
    fi
#--------------------------------------------------------------------------------------------------#

#---| LogRatio + BAF output File |-----------------------------------------------------------------#
    if [ "${RUN_MODE}" = "nt" ]; then
        ${RUN_CNVKIT} export nexus-ogt \
            -i ${TUMOR_SEQ_ID} \
            -n ${NORMAL_SEQ_ID} \
            -o ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.${MATCH_TAG}.logR.baf.txt \
            ${CNVKIT_RES_DIR}/${TUMOR_SEQ_ID}.${MATCH_TAG}.cnr \
            ${CNVKIT_VCF_TAG}.ready.vcf
    fi
#--------------------------------------------------------------------------------------------------#

#---| RE-CALL CN and CREATE SUMMARY |--------------------------------------------------------------#
    Rscript ${RUN_SCRIPT}/WES_run.CNV_CNVkit_Recall_and_Summary.R \
        --BASE_DIR ${BASE_DIR} \
        --SEQ_FOLDER ${SEQ_FOLDER} \
        --SEQ_ID ${TUMOR_SEQ_ID} \
        --ASSEMBLY ${GENOME_ASSEMBLY} \
        --THREADS ${THREADS} \
        --CNVKIT_RUN_MODE ${RUN_MODE} \
        --NORMAL_TYPE ${NORMAL_TYPE}
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | [ CNV CALL ] Finished." >> ${LOGFILE}
#==================================================================================================# 

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | CNV-CNVKIT FINISHED." >> ${LOGFILE}
    echo " " >> ${LOGFILE}
    echo "WES CNV-CNVkit DONE." >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " " >> ${LOGFILE}
#==================================================================================================# 

