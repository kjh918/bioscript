
#!/bin/bash

#==============================================================================#
# ScriptID : GMC-NGS-B02-A05
# Author : kangsm
# version : v1.0 
# lastest update date : 2025.12.26
#==============================================================================#

#---| Input Options |------------------------------------------------------------------------------#
    usage() {
        echo "Usage: WES_run.AnalysisReadyBAM_QC.sh [ ARGUMENTS ]...
        [ -d | --seqFolder ]      : data processing base directory
        [ -s | --seqID ]          : sample ID (Used as sample's processing product folder name)
        [ --baseDir ]             : Base Work Dir
        [ --threads ]             : N-CPU threads
        [ --assembly ]            : Reference genome assembly
        [ --intervalBait ]        : Probe interval_list
        [ --intervalTarget ]      : Target interval_list
        [ --targetBed ]           : Target BED file 
        [ --removeTempBAM ]       : Delete Temp BAM files  (default = true)
        [ --runMultiqc ]          : run MultiQC in this step or not ( default = true )
        [ -h | --help ]           : Print this usage 
        "
    }
#--------------------------------------------------------------------------------------------------#

#---| Set Defualt Values |-------------------------------------------------------------------------#
    BASE_DIR=/storage/home/kangsm/analysis
    THREADS=15
    GENOME_ASSEMBLY="hg19"
    REMOVE_TMP_BAM="true"
    RUN_MULTIQC="true"
#--------------------------------------------------------------------------------------------------#

#---| Parse Option Arguments |---------------------------------------------------------------------#
    ARGS=$(getopt -a -o d:s:h: --long seqFolder:,seqID:,baseDir:,threads:,assembly:,intervalBait:,intervalTarget:,removeTempBAM:,runMultiqc:,targetBed:,help -- "$@" )
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
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  BASE_DIR=$2 ; shift 2 ;;
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
        --intervalBait )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  BAIT_INTERVALS=$2 ; shift 2 ;;
            esac ;;
        --intervalTarget )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  TARGET_INTERVALS=$2 ; shift 2 ;;
            esac ;;
        --targetBed )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  TARGET_BED=$2 ; shift 2 ;;
            esac ;;
        --removeTempBAM )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = true" ; shift ;;
                *)  REMOVE_TMP_BAM=$2 ; shift 2 ;;
            esac ;;
        --runMultiqc )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = true" ; shift ;;
                *)  RUN_MULTIQC=$2 ; shift 2 ;;
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
    LOGFILE=${LOG_DIR}/${SEQ_ID}.analysisReady.bam.qc.$(date '+%Y%m%d').log
    if [ ! -f ${LOGFILE} ]; then touch ${LOGFILE}; fi
#--------------------------------------------------------------------------------------------------#

#---| SWTOOLS RUN |--------------------------------------------------------------------------------#
    RUN_SINGULARITY="singularity exec -B /storage,/data /storage/images"
    GATK4="${RUN_SINGULARITY}/gatk-4.4.0.0.sif gatk"
    GATK4_JAR_16G="${RUN_SINGULARITY}/gatk-4.4.0.0.sif java -XX:ParallelGCThreads=14 -Xmx16384m -jar /gatk/gatk-package-4.4.0.0-local.jar"
#--------------------------------------------------------------------------------------------------#

#---| REFERENCES |---------------------------------------------------------------------------------#
    if [ "${GENOME_ASSEMBLY}" = "hg38" ]; then
        # GENOME FASTA #
        GENOME_FASTA=/storage/references_and_index/hg38/fasta/Homo_sapiens_assembly38.fasta
    else
        # GENOME FASTA #
        GENOME_FASTA=/storage/references_and_index/hg19/fasta/human_g1k_v37_decoy.fasta
    fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo " " >> ${LOGFILE}
    echo "WES AnalysisReady BAM QC log" >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " RUN DATE             : $(date '+%Y-%m-%d %H:%M:%S')" >> ${LOGFILE}
    echo " SEQ FOLDER           : ${DATA_DIR}" >> ${LOGFILE}
    echo " SEQ ID               : ${SEQ_ID}" >> ${LOGFILE}
    echo " REFERENCE            : ${GENOME_ASSEMBLY}" >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " " >> ${LOGFILE}   
#==================================================================================================# 

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | ANALYSIS-READY BAM QC START." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') |   - [ QC1 : Sequecning Artifact ] Start." >> ${LOGFILE}
#==================================================================================================# 

#---| CHECK ANALYSIS-READY BAM |-------------------------------------------------------------------#
    if [ ! -f ${BAM_DIR}/${SEQ_ID}.analysisReady.bam ]; then
        cp ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.realign.recal.bam ${BAM_DIR}/${SEQ_ID}.analysisReady.bam
        cp ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.realign.recal.bai ${BAM_DIR}/${SEQ_ID}.analysisReady.bam.bai
    fi
#--------------------------------------------------------------------------------------------------#    

#---| COLLECT SEQUENCING ARTIFACTS |---------------------------------------------------------------#
    ${GATK4} CollectSequencingArtifactMetrics \
        --java-options "-XX:ParallelGCThreads=14 -Xmx16384m" \
        --INPUT ${BAM_DIR}/${SEQ_ID}.analysisReady.bam  \
        --OUTPUT ${QC_DIR}/${SEQ_ID}.artifact \
        --FILE_EXTENSION .txt \
        --REFERENCE_SEQUENCE ${GENOME_FASTA}
#--------------------------------------------------------------------------------------------------#    

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') |   - [ QC1 : Sequecning Artifact ] Finished." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') |   - [ QC2 : Alignment Summary ] Start." >> ${LOGFILE}
#==================================================================================================# 

#---| BAM FILES SUFFIX |---------------------------------------------------------------------------#
    BAM_SUFFIX="merged.bam merged.dup.marked.bam merged.dup.marked.realign.bam merged.dup.marked.realign.recal.bam"
#--------------------------------------------------------------------------------------------------#    
 
#---| ALIGNMENT SUMMARY |--------------------------------------------------------------------------#
    parallel -j 4 -k "${GATK4_JAR_16G} CollectAlignmentSummaryMetrics \
        --INPUT ${BAM_DIR}/${SEQ_ID}.{} \
        --OUTPUT ${QC_DIR}/${SEQ_ID}.{}.AlignSummary.metrics.txt" ::: ${BAM_SUFFIX}
#--------------------------------------------------------------------------------------------------#    

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') |   - [ QC2 : Alignment Summary ] Finished." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') |   - [ QC3 : HS-metrics ] Start." >> ${LOGFILE}
#==================================================================================================# 

#---| COLLECT HS-METRIC |--------------------------------------------------------------------------#
    parallel -j 4 -k "${GATK4_JAR_16G} CollectHsMetrics \
        --INPUT ${BAM_DIR}/${SEQ_ID}.{} \
        --OUTPUT ${QC_DIR}/${SEQ_ID}.{}.HS.metrics.txt \
        --BAIT_INTERVALS ${BAIT_INTERVALS} \
        --TARGET_INTERVALS ${TARGET_INTERVALS}" ::: ${BAM_SUFFIX}  
#--------------------------------------------------------------------------------------------------#    

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') |   - [ QC3 : HS-metrics ] Finished." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') |   - [ QC4 : InsertSize-metrics ] Start." >> ${LOGFILE}
#==================================================================================================# 

#---| COLLECT INSERT-SIZE-METRIC |-----------------------------------------------------------------#
    parallel -j 4 -k "${GATK4_JAR_16G} CollectInsertSizeMetrics \
        --INPUT ${BAM_DIR}/${SEQ_ID}.{} \
        --OUTPUT ${QC_DIR}/${SEQ_ID}.{}.InsertSize.metrics.txt \
        --Histogram_FILE ${QC_DIR}/${SEQ_ID}.{}.InsertSize.histogram.pdf \
        --REFERENCE_SEQUENCE ${GENOME_FASTA}" ::: ${BAM_SUFFIX}  
#--------------------------------------------------------------------------------------------------#    

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') |   - [ QC4 : InsertSize-metrics ] Finished." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') |   - [ QC5 : mosdepth ] Start." >> ${LOGFILE}
#==================================================================================================# 

#---| RUN MOSDEPTH |-------------------------------------------------------------------------------#
    ${RUN_SINGULARITY}/mosdepth-0.3.6.sif /opt/mosdepth \
        --threads ${THREADS} \
        --by ${TARGET_BED} \
        --mapq 20 \
        ${QC_DIR}/${SEQ_ID} \
        ${BAM_DIR}/${SEQ_ID}.analysisReady.bam 
#--------------------------------------------------------------------------------------------------#    

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') |   - [ QC5 : mosdepth ] Finished." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') |   - [ QC6 : alfred ] Start." >> ${LOGFILE}
#==================================================================================================# 

#---| RUN ALFRED |---------------------------------------------------------------------------------#
    ${RUN_SINGULARITY}/alfred-0.2.6.sif /opt/alfred/bin/alfred qc \
        --reference ${GENOME_FASTA} \
        --bed ${TARGET_BED} \
        --outfile ${QC_DIR}/${SEQ_ID}.alfred.qc.tsv.gz \
        ${BAM_DIR}/${SEQ_ID}.analysisReady.bam 
    ##
    zgrep ^CM ${QC_DIR}/${SEQ_ID}.alfred.qc.tsv.gz | cut -f 2- > ${QC_DIR}/${SEQ_ID}.alfred.chr.map.stats.txt
    zgrep ^TC ${QC_DIR}/${SEQ_ID}.alfred.qc.tsv.gz | cut -f 2- > ${QC_DIR}/${SEQ_ID}.alfred.target.coverage.txt
#--------------------------------------------------------------------------------------------------#    

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') |   - [ QC6 : alfred ] Finished." >> ${LOGFILE}
#==================================================================================================# 

#---| QC RESULT SUMMARY |--------------------------------------------------------------------------#
    if [ "${RUN_MULTIQC}" = "true" ]; then
        #==========================================================================================#
            echo "$(date '+%Y-%m-%d %H:%M:%S') | MULTIQC RESULT SUMMARY START." >> ${LOGFILE}
        #==========================================================================================#
        ##
        MULTIQC_CONFIG=/storage/home/kangsm/runScripts/NGS_config.MultiQC_Custom.yaml
        ##
        ${RUN_SINGULARITY}/multiqc-1.16.sif multiqc --force \
            --filename ${SEQ_ID}_qc \
            --outdir ${QC_DIR}/qc_res \
            --config ${MULTIQC_CONFIG} \
            --data-dir ${QC_DIR}
        ##
        #==========================================================================================# 
            echo "$(date '+%Y-%m-%d %H:%M:%S') | MULTIQC RESULT SUMMARY FINISHED." >> ${LOGFILE}
        #==========================================================================================#
    else
        #==========================================================================================# 
            echo "$(date '+%Y-%m-%d %H:%M:%S') | MULTIQC RESULT SUMMARY SKIPPED." >> ${LOGFILE}
        #==========================================================================================#
    fi
#--------------------------------------------------------------------------------------------------#    

#---| DELETE TEMP-BAM FILES |----------------------------------------------------------------------#
    if [ "${REMOVE_TMP_BAM}" = "true" ]; then
        rm ${BAM_DIR}/${SEQ_ID}.merged.bam   
        rm ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.bam
        rm ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.realign.bam
        rm ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.realign.recal.bam
        rm ${BAM_DIR}/${SEQ_ID}.merged.bai   
        rm ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.bai
        rm ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.realign.bai
        rm ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.realign.recal.bai
        rm ${BAM_DIR}/${SEQ_ID}.bwa.mem.sam
        rm ${BAM_DIR}/${SEQ_ID}.fastqtosam.bam
        #==========================================================================================# 
            echo "$(date '+%Y-%m-%d %H:%M:%S') | REMOVE TMP-BAM FILES DONE." >> ${LOGFILE}
        #==========================================================================================#
    else
        #==========================================================================================# 
            echo "$(date '+%Y-%m-%d %H:%M:%S') | REMOVE TMP-BAM FILES SKIPPED." >> ${LOGFILE}
        #==========================================================================================#
    fi
#--------------------------------------------------------------------------------------------------#    

#==================================================================================================# 
    echo " " >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | ANALYSIS-READY BAM QC DONE." >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " " >> ${LOGFILE}
#==================================================================================================# 


