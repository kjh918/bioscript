

#!/bin/bash

# script version v0.5

#---| Input Options |------------------------------------------------------------------------------#
    usage() {
        echo "Usage: WES_run.CNV_PureCN.sh [ ARGUMENTS ]...
        [ --seqFolder ]      : data processing base directory
        [ --tumorSeqID ]     : Tumor Sample ID 
        [ --normalSeqID ]    : Matched Normal Sample ID
        [ --runMode ]        : Run mode, 'nt' or 'tonly'.
        [ --baseDir ]        : Base Work Dir
        [ --assembly ]       : Genome assembly ( hg19 = GRCh37 )
        [ --wesLibKit]       : WES Library Kit ( default = twist exome 2.0)
        [ --callMinDepth ]   : Callable BED minimum depth cut-off (default = 25)
        [ --nType ]          : VCF TAG
        [ --forceCoverage ]  : force re-run BAM coverage. default = false
        [ -h | --help ]      : Print this usage         
        "
    }
#--------------------------------------------------------------------------------------------------#

#---| Set Defualt Values |-------------------------------------------------------------------------#
    BASE_DIR=/data/wes
    GENOME_ASSEMBLY="hg19"
    WES_LIB_KIT="twist.exome.2.0"
    CALL_MIN_DEPTH=20
    RUN_MODE="tonly"
    NORMAL_TYPE="TS"
    FORCE_COVERAGE="false"
#--------------------------------------------------------------------------------------------------#

#---| Parse Option Arguments |---------------------------------------------------------------------#
    ARGS=$(getopt -a -o h: --long seqFolder:,tumorSeqID:,normalSeqID:,wesLibKit:,baseDir:,assembly:,runMode:,callMinDepth:,forceCoverage:,nType:,help -- "$@" )
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
        -t | --tumorSeqID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  TUMOR_SEQ_ID=$2 ; shift 2 ;;
            esac  ;;
        -n | --normalSeqID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  NORMAL_SEQ_ID=$2 ; shift 2 ;;
            esac  ;;
        --baseDir )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = /data/wes" ; shift ;;
                *)  BASE_DIR=$2 ; shift 2 ;;
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
        --callMinDepth )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 20" ; shift ;;
                *)  CALL_MIN_DEPTH=$2 ; shift 2 ;;
            esac ;;
        --runMode )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = tonly" ; shift ;;
                *)  RUN_MODE=$2 ; shift 2 ;;
            esac ;;
        --nType )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = TS" ; shift ;;
                *)  NORMAL_TYPE=$2 ; shift 2 ;;
            esac ;;
        --forceCoverage )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = false" ; shift ;;
                *)  FORCE_COVERAGE=$2 ; shift 2 ;;
            esac ;;
        -h | --help )
        usage >&2 ; exit 2 ;;   
        --) shift ; break ;;
        *)  usage >&2 ; exit 2 ;;
    esac
    done
#---/

#---| FOLDERS |------------------------------------------------------------------------------------#
    DATA_DIR=${BASE_DIR}/${SEQ_FOLDER}
    TUMOR_BAM_DIR=${DATA_DIR}/${TUMOR_SEQ_ID}/bam
    NORMAL_BAM_DIR=${DATA_DIR}/${NORMAL_SEQ_ID}/bam
    VCF_DIR=${DATA_DIR}/${TUMOR_SEQ_ID}/vcf
    LOG_DIR=${DATA_DIR}/${TUMOR_SEQ_ID}/log
    #--------------------------------------------------------------------------#  
    CNV_RES_DIR=${DATA_DIR}/${TUMOR_SEQ_ID}/cnv
    if [ ! -d ${CNV_RES_DIR} ]; then mkdir -p ${CNV_RES_DIR}; fi
    #--------------------------------------------------------------------------#  
    PURECN_RES_DIR=${CNV_RES_DIR}/purecn
    if [ ! -d ${PURECN_RES_DIR} ]; then mkdir -p ${PURECN_RES_DIR}; fi
    #--------------------------------------------------------------------------#  
    PURECN_RESOURCES=/storage/references_and_index/cnv/purecn
#--------------------------------------------------------------------------------------------------#

#---| SWTOOLS RUN |--------------------------------------------------------------------------------#   
    RUN_PURECN=/storage/apps/R/4.3.1/lib/R/library/PureCN/extdata
    RUN_SCRIPT=/storage/home/kangsm/runScripts
    RUN_SINGULARITY="singularity exec -B /storage,/data /storage/images"
    GATK3="${RUN_SINGULARITY}/gatk-3.8-1.sif java -Xmx16384m -jar /usr/GenomeAnalysisTK.jar"
#--------------------------------------------------------------------------------------------------#

#---| LOGFILE |------------------------------------------------------------------------------------#
    LOGFILE=${LOG_DIR}/${TUMOR_SEQ_ID}.CNV.purecn.$(date '+%Y%m%d').log
    if [ ! -f ${LOGFILE} ]; then touch ${LOGFILE}; fi
#--------------------------------------------------------------------------------------------------#

#---| REFERENCES |---------------------------------------------------------------------------------#
    if [ "${GENOME_ASSEMBLY}" = "hg38" ]; then
        GENOME_FASTA=/storage/references_and_index/hg38/fasta/Homo_sapiens_assembly38.fasta
    else
        GENOME_FASTA=/storage/references_and_index/hg19/fasta/human_g1k_v37_decoy.fasta
    fi
    #--------------------------------------------------------------------------#
    SIMPLE_REPEAT=${PURECN_RESOURCES}/snp_blacklist/${GENOME_ASSEMBLY}_simple.repeat.bed
    INTERVALS=${PURECN_RESOURCES}/interval_files/${GENOME_ASSEMBLY}_${WES_LIB_KIT}.baits_intervals.txt
    PON_RDS=${PURECN_RESOURCES}/normalDB/${GENOME_ASSEMBLY}_${WES_LIB_KIT}/normalDB_${GENOME_ASSEMBLY}_${WES_LIB_KIT}_${GENOME_ASSEMBLY}.rds
#--------------------------------------------------------------------------------------------------#

#---| NORMAL PANEL TAG |---------------------------------------------------------------------------#
    if [ "${RUN_MODE}" = "nt" ]; then 
        NORMAL_PANEL=${NORMAL_SEQ_ID}
        RUN_MODE_TAG="NT-Mode"
    else
        NORMAL_PANEL="Panel of Normal"
        RUN_MODE_TAG="Tumor-Only Mode"
    fi
#--------------------------------------------------------------------------------------------------#

#---| VCF TAG |------------------------------------------------------------------------------------#
    if [ "${RUN_MODE}" = "tonly" ]; then
        VCF_TAG="mutect2.keep.germline.filtered"
        RES_FOLDER_TAG="Tonly"
    else
        if [ "${NORMAL_TYPE}" = "ORG" ]; then
            VCF_TAG="mutect2.ORG.NT.keep.germline.filtered"
            RES_FOLDER_TAG="ORG.NT"
        else
            VCF_TAG="mutect2.NT.keep.germline.filtered"
            RES_FOLDER_TAG="NT"
        fi
    fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo " " >> ${LOGFILE}
    echo "WES PureCN Purity and CNV Analysis log" >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " RUN DATE             : $(date '+%Y-%m-%d %H:%M:%S')" >> ${LOGFILE}
    echo " SEQ FOLDER           : ${DATA_DIR}" >> ${LOGFILE}
    echo " TUMOR  SEQ ID        : ${TUMOR_SEQ_ID}" >> ${LOGFILE}
    echo " NORMAL SEQ ID        : ${NORMAL_PANEL}" >> ${LOGFILE}
    echo " REFERENCE            : ${GENOME_ASSEMBLY}" >> ${LOGFILE}
    echo " RUN MODE             : ${RUN_MODE_TAG}" >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " " >> ${LOGFILE}
#==================================================================================================# 

#---| CHECK PURECN VCF |---------------------------------------------------------------------------#
    ##
    TUMOR_VCF=${VCF_DIR}/${TUMOR_SEQ_ID}.${VCF_TAG}.vcf
    ##
    #--------------------------------------------------------------------------#    
    if [ ! -f ${TUMOR_VCF} ]; then
        #==========================================================================================# 
            echo "$(date '+%Y-%m-%d %H:%M:%S') | NO VCF for PureCN Analysis. STOPPED." >> ${LOGFILE}
        #==========================================================================================# 

        echo "NO VCF for PureCN. Stopped."
        exit 0
    fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | PURECN ANALYSIS START." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | [ Internal Segmentation : Tumor ] Start." >> ${LOGFILE}
#==================================================================================================# 

#---| INTERNAL SEGMENTATION : TUMOR |--------------------------------------------------------------#
    TUMOR_COVERAGE=${PURECN_RES_DIR}/${TUMOR_SEQ_ID}.analysisReady_coverage_loess.txt.gz
    #--------------------------------------------------------------------------#    
    if [ ! -f ${TUMOR_COVERAGE} ]; then
        Rscript ${RUN_PURECN}/Coverage_rev.R \
            --out-dir ${PURECN_RES_DIR} \
            --bam ${TUMOR_BAM_DIR}/${TUMOR_SEQ_ID}.analysisReady.bam \
            --intervals ${INTERVALS}
    else
        if [ "${FORCE_COVERAGE}" = "true" ]; then

            #==================================================================================# 
                echo "$(date '+%Y-%m-%d %H:%M:%S') | Re-Run Tumor Sample Coverage." >> ${LOGFILE}
            #==================================================================================#

            Rscript ${RUN_PURECN}/Coverage_rev.R \
                --out-dir ${PURECN_RES_DIR} \
                --bam ${TUMOR_BAM_DIR}/${TUMOR_SEQ_ID}.analysisReady.bam \
                --intervals ${INTERVALS}
        fi
    fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | [ Internal Segmentation : Tumor ] Finished." >> ${LOGFILE}
#==================================================================================================# 

#---| INTERNAL SEGMENTATION : NORMAL |-------------------------------------------------------------#
    if [ "${RUN_MODE}" = "nt" ]; then 
        #==========================================================================================# 
            echo "$(date '+%Y-%m-%d %H:%M:%S') | Run Mode = NT-Mode. Check for Normal Sample Coverage." >> ${LOGFILE}
        #==========================================================================================# 

        NORMAL_COVERAGE=${DATA_DIR}/${NORMAL_SEQ_ID}/cnv/purecn/${NORMAL_SEQ_ID}.analysisReady_coverage_loess.txt.gz
        #----------------------------------------------------------------------#
        if [ ! -f ${NORMAL_COVERAGE} ]; then
            #======================================================================================# 
                echo "$(date '+%Y-%m-%d %H:%M:%S') | NO Normal Sample Coverage." >> ${LOGFILE}
            #======================================================================================#

            NORMAL_COV_DIR=${DATA_DIR}/${NORMAL_SEQ_ID}/cnv/purecn
            if [ ! -d ${NORMAL_COV_DIR} ]; then mkdir -p ${NORMAL_COV_DIR}; fi

            #======================================================================================# 
                echo "$(date '+%Y-%m-%d %H:%M:%S') | [ Internal Segmentation : Normal ] Start." >> ${LOGFILE}
            #======================================================================================#

            Rscript ${RUN_PURECN}/Coverage_rev.R \
                --out-dir ${DATA_DIR}/${NORMAL_SEQ_ID}/cnv/purecn \
                --bam ${DATA_DIR}/${NORMAL_SEQ_ID}/bam/${NORMAL_SEQ_ID}.analysisReady.bam \
                --intervals ${INTERVALS}

            #======================================================================================# 
                echo "$(date '+%Y-%m-%d %H:%M:%S') | [ Internal Segmentation : Normal ] Finished." >> ${LOGFILE}
            #======================================================================================#
        else
            #======================================================================================# 
                echo "$(date '+%Y-%m-%d %H:%M:%S') | FOUND Normal Sample Coverage." >> ${LOGFILE}
            #======================================================================================#

            if [ "${FORCE_COVERAGE}" = "true" ]; then
                #==================================================================================# 
                    echo "$(date '+%Y-%m-%d %H:%M:%S') | Re-Run Normal Sample Coverage." >> ${LOGFILE}
                #==================================================================================#

                Rscript ${RUN_PURECN}/Coverage_rev.R \
                    --out-dir ${DATA_DIR}/${NORMAL_SEQ_ID}/cnv/purecn \
                    --bam ${DATA_DIR}/${NORMAL_SEQ_ID}/bam/${NORMAL_SEQ_ID}.analysisReady.bam \
                    --intervals ${INTERVALS}
            fi           
        fi
    else
        #==========================================================================================# 
            echo "$(date '+%Y-%m-%d %H:%M:%S') | Run Mode = Tumor-Ony Mode. Check Normal Sample Coverage Skipped." >> ${LOGFILE}
        #==========================================================================================# 
    fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | PureCN Analysis Start." >> ${LOGFILE}
#==================================================================================================# 

#---| PureCN Analysis |----------------------------------------------------------------------------#
    if [ "${RUN_MODE}" = "nt" ]; then
        #==========================================================================================# 
            echo "$(date '+%Y-%m-%d %H:%M:%S') | [ PureCN : NT-Mode ] Start." >> ${LOGFILE}
        #==========================================================================================#

        mkdir -p ${PURECN_RES_DIR}/${RES_FOLDER_TAG}
        #----------------------------------------------------------------------#
        Rscript ${RUN_PURECN}/PureCN_rev.R \
            --out ${PURECN_RES_DIR}/${RES_FOLDER_TAG} \
            --tumor ${TUMOR_COVERAGE} \
            --normal ${NORMAL_COVERAGE} \
            --normaldb ${PON_RDS} \
            --fun-segmentation PSCBS \
            --sampleid ${TUMOR_SEQ_ID} \
            --vcf ${TUMOR_VCF} \
            --intervals ${INTERVALS} \
            --genome ${GENOME_ASSEMBLY} \
            --snp-blacklist ${SIMPLE_REPEAT}

        #==========================================================================================# 
            echo "$(date '+%Y-%m-%d %H:%M:%S') | [ PureCN : NT-Mode ] Finished." >> ${LOGFILE}
        #==========================================================================================#
    else
        #==========================================================================================# 
            echo "$(date '+%Y-%m-%d %H:%M:%S') | [ PureCN : Tunor-Only Mode ] Start." >> ${LOGFILE}
        #==========================================================================================#

        mkdir -p ${PURECN_RES_DIR}/${RES_FOLDER_TAG}
        #----------------------------------------------------------------------#
        Rscript ${RUN_PURECN}/PureCN_rev.R \
            --out ${PURECN_RES_DIR}/${RES_FOLDER_TAG} \
            --tumor ${TUMOR_COVERAGE} \
            --normaldb ${PON_RDS} \
            --fun-segmentation PSCBS \
            --sampleid ${TUMOR_SEQ_ID} \
            --vcf ${TUMOR_VCF} \
            --intervals ${INTERVALS} \
            --genome ${GENOME_ASSEMBLY} \
            --snp-blacklist ${SIMPLE_REPEAT}

        #==========================================================================================# 
            echo "$(date '+%Y-%m-%d %H:%M:%S') | [ PureCN : Tunor-Only Mode ] Finished." >> ${LOGFILE}
        #==========================================================================================#
    fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | PureCN Analysis Finished." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | [ Callable Region Processing ] Start." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') |  - Callable Region Minimum Depth = ${CALL_MIN_DEPTH}" >> ${LOGFILE}
#==================================================================================================# 

#---| CALLABLE REGION PROCESSING |-----------------------------------------------------------------#
    CALLABLE_CDS_BED=${PURECN_RES_DIR}/${TUMOR_SEQ_ID}.minDP.${CALL_MIN_DEPTH}.status.CDS.bed
    #--------------------------------------------------------------------------#
    if [ ! -f ${CALLABLE_CDS_BED} ]; then
        ${GATK3} -T CallableLoci \
            -R ${GENOME_FASTA} \
            -I:tumor ${TUMOR_BAM_DIR}/${TUMOR_SEQ_ID}.analysisReady.bam  \
            --summary ${PURECN_RES_DIR}/${TUMOR_SEQ_ID}.minDP.${CALL_MIN_DEPTH}.callable.table.txt \
            -o ${PURECN_RES_DIR}/${TUMOR_SEQ_ID}.minDP.${CALL_MIN_DEPTH}.callable.status.bed \
            --minDepth ${CALL_MIN_DEPTH}
        #----------------------------------------------------------------------#
        grep CALLABLE ${PURECN_RES_DIR}/${TUMOR_SEQ_ID}.minDP.${CALL_MIN_DEPTH}.callable.status.bed > \
        ${PURECN_RES_DIR}/${TUMOR_SEQ_ID}.minDP.${CALL_MIN_DEPTH}.status.filtered.bed
        #----------------------------------------------------------------------#
        Rscript ${RUN_PURECN}/FilterCallableLoci.R \
            --genome ${GENOME_ASSEMBLY} \
            --in-file ${PURECN_RES_DIR}/${TUMOR_SEQ_ID}.minDP.${CALL_MIN_DEPTH}.status.filtered.bed \
            --out-file ${CALLABLE_CDS_BED}
    fi   
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | [ Callable Region Processing ] Finished." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | [ Puriry and CNV Calling ] Start." >> ${LOGFILE}
#==================================================================================================# 

#---| PURITY & CNV CALL |--------------------------------------------------------------------------#
        Rscript ${RUN_PURECN}/Dx.R \
            --out ${PURECN_RES_DIR}/${RES_FOLDER_TAG}/${TUMOR_SEQ_ID} \
            --rds ${PURECN_RES_DIR}/${RES_FOLDER_TAG}/${TUMOR_SEQ_ID}.rds \
            --callable ${CALLABLE_CDS_BED} \
            --exclude ${SIMPLE_REPEAT} \
            --signatures \
            --force    
#--------------------------------------------------------------------------------------------------#

#---| PURITY & PLOIDY RESULT TEXT FILE |-----------------------------------------------------------#
    Rscript ${RUN_SCRIPT}/WES_run.CNV_PureCN_Summary.R \
        --baseDir ${BASE_DIR} \
        --seqFolder ${SEQ_FOLDER} \
        --seqID ${TUMOR_SEQ_ID} \
        --resultFolder ${RES_FOLDER_TAG}
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | [ Puriry and CNV Calling ] Finished." >> ${LOGFILE}
    echo " " >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | WES PureCN Purity and CNV Analysis DONE." >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " " >> ${LOGFILE}
#==================================================================================================# 



