
#!/bin/bash

#==============================================================================#
# ScriptID : GMC-NGS-B02-A05
# Author : kangsm
# version : v1.0 
# lastest update date : 2025.12.26
#==============================================================================#

#---| Input Options |------------------------------------------------------------------------------#
    usage() {
        echo "
        Usage: WGS_P04-1_BAM-QC.sh [ ARGUMENTS ]...
        
        --BASE_DIR          [ REQUIRED ] analysis base folder by NGS data type
        
        --BATCH_ID          [ REQUIRED ] dataset batch-id
        
        --SEQ_ID            [ REQURIED ] sample-id 
                
        --THREADS           parallel processing cpu cores 
                            default = 15
        
        --GENOME_ASSEMBLY   Reference genome assembly version
                            default = hg38
        
        --BAM_TYPE          bam file type. 
                            hmftools or gatk.
                            default = hmftools
       
        --RUN_MULTIQC       run Multi-QC in this step or not
                            default = true
        
        -h | --help    : Print this usage 
        "
    }
#--------------------------------------------------------------------------------------------------#

#---| Parse Option Arguments |---------------------------------------------------------------------#
    ARGS=$(getopt -a -o d:s:h: --long BATCH_ID:,SEQ_ID:,BASE_DIR:,THREADS:,GENOME_ASSEMBLY:,BAM_TYPE:,RUN_MULTIQC:,help -- "$@" )
    VALID_ARGS=$?
    if [ "$VALID_ARGS" != "0" ]; then 
        usage >&2 
        exit 2
    fi

    eval set -- "$ARGS"
    while :
    do
    case "$1" in
        -d | --BATCH_ID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  BatchID=$2 ; shift 2 ;;
            esac ;;
        -s | --SEQ_ID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  SeqID=$2 ; shift 2 ;;
            esac  ;;
        --BASE_DIR )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  BaseDir=$2 ; shift 2 ;;
            esac ;;
        --THREADS )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 15" ; shift ;;
                *)  Threads=$2 ; shift 2 ;;
            esac ;;
        --GENOME_ASSEMBLY )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = hg19" ; shift ;;
                *)  GenomeAssembly=$2 ; shift 2 ;;
            esac ;;
        --BAM_TYPE )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = gatk" ; shift ;;
                *)  BamType=$2 ; shift 2 ;;
            esac ;;
        --RUN_MULTIQC )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = true" ; shift ;;
                *)  RunMultiqc=$2 ; shift 2 ;;
            esac ;;
        -h | --help )
        usage >&2 ; exit 2 ;;   
        --) shift ; break ;;
        *)  usage >&2 ; exit 2 ;;
    esac
    done
#--------------------------------------------------------------------------------------------------#

#---| Param Check and Set Default |-----------------------------------------------------------------
    if [ -z "${BaseDir}" ]; then BaseDir="/data/wgs"; fi
    if [ -z "${Threads}" ]; then Threads=15; fi
    if [ -z "${GenomeAssembly}" ]; then GenomeAssembly="hg38"; fi
    if [ -z "${BamType}" ]; then BamType="hmftools"; fi
    if [ -z "${ReadLength}" ]; then ReadLength=150; fi
    if [ -z "${RunMultiqc}" ]; then RunMultiqc=true; fi
#---------------------------------------------------------------------------------------------------

#---| Manual Params |------------------------------------------------------------------------------#
    # BaseDir="/data/wgs/RUO_Service"
    # BatchID="WGS-26-02"
    # SeqID="WGS-26-02_01"
    # Threads=12
    # GenomeAssembly="hg19"
    # BamType="hmftools"
    # RunMultiqc="true"
#--------------------------------------------------------------------------------------------------#


#---| FOLDERS |------------------------------------------------------------------------------------#
    DataDir=${BaseDir}/${BatchID}
    ##  
    BamDir=${DataDir}/${SeqID}/bam
    QcDir=${DataDir}/${SeqID}/qcfiles
    TmpDir=${DataDir}/${SeqID}/tmp
    LogDir=${DataDir}/${SeqID}/log
    ##
    if [ ! -d ${BamDir} ]; then mkdir -p ${BamDir}; fi
    if [ ! -d ${TmpDir} ]; then mkdir -p ${TmpDir}; fi
    if [ ! -d ${LogDir} ]; then mkdir -p ${LogDir}; fi
    if [ ! -d ${QcDir} ]; then mkdir -p ${QcDir}; fi
#--------------------------------------------------------------------------------------------------#

#---| LogFile |------------------------------------------------------------------------------------#
    LogFile=${LogDir}/${SeqID}.analysisReady.bam.qc.$(date '+%Y%m%d').log
    if [ ! -f ${LogFile} ]; then touch ${LogFile}; fi
#--------------------------------------------------------------------------------------------------#

#---| SWTOOLS RUN |--------------------------------------------------------------------------------#
    export LC_ALL=C
    export LANG=C
    GATK_VERSION="4.6.2.0"
    RunSingularity="singularity exec -B /storage,/data /storage/images"
    RunGatk4="${RunSingularity}/gatk-${GATK_VERSION}.sif gatk"
    RunMosdepth="${RunSingularity}/mosdepth-0.3.6.sif /opt/mosdepth"
    RunAlfredQc="${RunSingularity}/alfred-0.2.6.sif /opt/alfred/bin/alfred qc"
#--------------------------------------------------------------------------------------------------#

#---| RESOURCES |----------------------------------------------------------------------------------#
    WGS_ChromIntervals=/storage/references_and_index/${GenomeAssembly}/bed/${GenomeAssembly}.wgs.call.target.region.interval_list
    WGS_ChromBED=/storage/references_and_index/${GenomeAssembly}/bed/${GenomeAssembly}.wgs.call.target.region.bed 
#--------------------------------------------------------------------------------------------------#

# #==================================================================================================# 
#     echo " " >> ${LogFile}
#     echo "WGS AnalysisReady BAM QC log" >> ${LogFile}
#     echo "----------------------------------------------------------------------" >> ${LogFile}
#     echo " RUN DATE             : $(date '+%Y-%m-%d %H:%M:%S')" >> ${LogFile}
#     echo " BATCH ID             : ${DataDir}" >> ${LogFile}
#     echo " SEQ ID               : ${SeqID}" >> ${LogFile}
#     echo " REFERENCE            : ${GenomeAssembly}" >> ${LogFile}
#     echo "----------------------------------------------------------------------" >> ${LogFile}
#     echo " " >> ${LogFile}   
# #==================================================================================================# 

# #==================================================================================================# 
#     echo "$(date '+%Y-%m-%d %H:%M:%S') | ANALYSIS-READY BAM QC START." >> ${LogFile}
#     echo "$(date '+%Y-%m-%d %H:%M:%S') |   - [ QC1 : Sequecning Artifact ] Start." >> ${LogFile}
# #==================================================================================================# 

### INPUT DATA
    if [[ "${BamType}" = "hmftools" ]]; then
        InputBam="${BaseDir}/${BatchID}/${SeqID}/hmftools/redux/${SeqID}.redux.bam"
        if [[ "${GenomeAssembly}" = "hg19" ]]; then
            GenomeFasta="/storage/references_and_index/hg19/hmftools/Homo_sapiens.GRCh37.GATK.illumina.fasta"
        else
            GenomeFasta="/storage/references_and_index/hg38/hmftools/GRCh38_masked_exclusions_alts_hlas.fasta"
        fi
        OutputPrefix="${SeqID}.redux"
        MultiqcOutputDir="${QcDir}/${OutputPrefix}.qc_res"
        mkdir -p "${MultiqcOutputDir}"
    else
        GenomeFasta=/storage/references_and_index/${GenomeAssembly}/fasta/${GenomeAssembly}.primary_decoy_HLA_no.ALT.fa
        ${RunGatk4} --java-options "-Xmx16g" PrintReads \
            --reference "${GenomeFasta}" \
            --input "${BamDir}/${SeqID}.analysisReady.bam" \
            --intervals "${WGS_ChromIntervals}" \
            --output "${BamDir}/${SeqID}.analysisReady.primary.only.bam" \
            --create-output-bam-index true 
        InputBam="${BamDir}/${SeqID}.analysisReady.primary.only.bam"        
        OutputPrefix="${SeqID}.gatk"
        MultiqcOutputDir="${QcDir}/${OutputPrefix}.qc_res"
        mkdir -p "${MultiqcOutputDir}"
    fi
###
### GATK WGS-METRICS
    if [[ "${BamType}" = "hmftools" ]]; then
        ln -s ${DataDir}/${SeqID}/hmftools/redux/${OutputPrefix}.wgs.metrics.txt ${QcDir}/${OutputPrefix}.wgs.metrics.txt
    else
        ${RunGatk4} --java-options "-Xmx16g" CollectWgsMetrics \
            --INPUT "${InputBam}" \
            --OUTPUT "${QcDir}/${OutputPrefix}.wgs.metrics.txt" \
            --INTERVALS "${WGS_ChromIntervals}" \
            --REFERENCE_SEQUENCE "${GenomeFasta}" \
            --VALIDATION_STRINGENCY LENIENT
    fi
###
### ALIGNMENT SUMMARY
    ${RunGatk4} --java-options "-Xmx16g" CollectAlignmentSummaryMetrics \
        --INPUT "${InputBam}" \
        --REFERENCE_SEQUENCE "${GenomeFasta}" \
        --OUTPUT "${QcDir}/${OutputPrefix}.align.summary.metrics.txt" \
        --VALIDATION_STRINGENCY LENIENT
###
### INSERT-SIZE METRICS
    ${RunGatk4} --java-options "-Xmx16g" CollectInsertSizeMetrics \
        --INPUT "${InputBam}" \
        --OUTPUT "${QcDir}/${OutputPrefix}.insert.size.metrics.txt" \
        --Histogram_FILE "${QcDir}/${OutputPrefix}.insert.size.histogram.pdf" \
        --REFERENCE_SEQUENCE "${GenomeFasta}"
###
### MOSDEPTH
    ${RunMosdepth} \
        --threads "${Threads}" \
        --by "${WGS_ChromBED}" \
        --mapq 20 --flag 1796 --no-per-base \
        "${QcDir}/${OutputPrefix}" \
        "${InputBam}" 
###
### ALFRED
    ${RunAlfredQc} \
        --reference "${GenomeFasta}" \
        --outfile "${QcDir}/${OutputPrefix}.alfred.qc.tsv.gz" \
        --jsonout "${QcDir}/${OutputPrefix}.alfred.qc.json.gz" \
        "${InputBam}"
    ##
    zgrep ^CM "${QcDir}/${OutputPrefix}.alfred.qc.tsv.gz" | cut -f 2- > "${QcDir}/${OutputPrefix}.alfred.chr.map.stats.txt"
    zgrep ^TC "${QcDir}/${OutputPrefix}.alfred.qc.tsv.gz" | cut -f 2- > "${QcDir}/${OutputPrefix}.alfred.target.coverage.txt"
###
### QC RESULT SUMMARY 
    if [[ "${RunMultiqc}" = "true" || "${RunMultiqc}" = "True" ]]; then
        MultiqcConfig="/storage/home/kangsm/runScripts/Pipeline_WGS/WGS_P05-2_QC-POST_PROCESSING-MULTIQC-CONFIG.yaml"
        ${RunSingularity}/multiqc-1.16.sif multiqc --force \
            --filename "${SeqID}_qc" \
            --outdir "${MultiqcOutputDir}" \
            --config "${MultiqcConfig}" \
            --data-dir "${QcDir}"
    fi
###


