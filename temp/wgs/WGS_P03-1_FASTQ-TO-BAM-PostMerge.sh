
#!/bin/bash

#==============================================================================#
# ScriptID : GMC-NGS-B02-A04
# Author : kangsm
# version : v1.0 
# lastest update date : 2025.12.26
#==============================================================================#

#---| Input Options |-------------------------------------------------------------------------------
    usage() {
        echo "
        Usage: WGS_P03-1_FASTQ-TO-BAM.sh [ ARGUMENTS ]...
        
        --BASE_DIR          [ REQUIRED ] analysis base folder by NGS data type
        
        --BATCH_ID          [ REQUIRED ] dataset batch-id
        
        --SEQ_ID            [ REQURIED ] sample-id
        
        --THREADS           parallel processing cpu cores 
                            default = 15
        
        --GENOME_ASSEMBLY   reference genome GENOME_ASSEMBLY version
                            default = hg38
               
        -h | --help         Print this usage 
        "
    }
#---------------------------------------------------------------------------------------------------

#---| Parse Option Arguments |----------------------------------------------------------------------
    ARGS=$(getopt -a -o h: --long BATCH_ID:,SEQ_ID:,BASE_DIR:,THREADS:,GENOME_ASSEMBLY:,help -- "$@" )
    VALID_ARGS=$?
    if [ "$VALID_ARGS" != "0" ]; then 
        usage >&2 
        exit 2
    fi

    eval set -- "$ARGS"
    while :
    do
    case "$1" in
        --BASE_DIR )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  BaseDir=$2 ; shift 2 ;;
            esac ;;
        --BATCH_ID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  BatchID=$2 ; shift 2 ;;
            esac ;;
        --SEQ_ID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  SeqID=$2 ; shift 2 ;;
            esac  ;;
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
        -h | --help )
        usage >&2 ; exit 2 ;;   
        --) shift ; break ;;
        *)  usage >&2 ; exit 2 ;;
    esac
    done
#---------------------------------------------------------------------------------------------------

#---| Check Params and Set Default -----------------------------------------------------------------
    if [ -z "${BaseDir}" ]; then exit 2; fi
    if [ -z "${BatchID}" ]; then exit 2; fi
    if [ -z "${SeqID}" ]; then exit 2; fi
    if [ -z "${Threads}" ]; then Threads=15; fi
    if [ -z "${GenomeAssembly}" ]; then GenomeAssembly="hg38"; fi
#---------------------------------------------------------------------------------------------------

#---| Manual Params |-------------------------------------------------------------------------------
    # BaseDir=/data/wgs
    # BatchID=cancer_cell_lines
    # SeqID=
    # Threads=15
    # GenomeAssembly="hg38"
#---------------------------------------------------------------------------------------------------

#---| FOLDERS |-------------------------------------------------------------------------------------
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
    if [ ! -d ${QcDir} ]; then mkdir -p ${QcDir}/qc_res; fi
#---------------------------------------------------------------------------------------------------

#---| SWTOOLS RUN |---------------------------------------------------------------------------------
    GATK_VERSION="4.6.2.0"
    RunSingularity="singularity exec -B /storage,/data /storage/images"
    RunGatk4="singularity exec -B /storage,/data /storage/images/gatk-${GATK_VERSION}.sif gatk"
    RunBwa="singularity exec -B /storage,/data /storage/images/bwa-0.7.17.sif bwa mem"
#---------------------------------------------------------------------------------------------------

#---| RESOURCES |-----------------------------------------------------------------------------------
    GenomeFasta=/storage/references_and_index/${GenomeAssembly}/fasta/${GenomeAssembly}.primary_decoy_HLA_no.ALT.fa
    BwaIndex=/storage/references_and_index/${GenomeAssembly}/fasta/${GenomeAssembly}.primary_decoy_HLA_no.ALT.fa
    Vcf_dbsnp=/storage/references_and_index/${GenomeAssembly}/vcf/${GenomeAssembly}.dbsnp138.vcf.gz
    Vcf_hcsnp=/storage/references_and_index/${GenomeAssembly}/vcf/${GenomeAssembly}.1000g.phase1.snps.high.confidence.vcf.gz
    Vcf_knwonindels=/storage/references_and_index/${GenomeAssembly}/vcf/${GenomeAssembly}.known.indels.vcf.gz
    Vcf_gsindels=/storage/references_and_index/${GenomeAssembly}/vcf/${GenomeAssembly}.gold.standard.indels.vcf.gz
    Vcf_gnomad=/storage/references_and_index/${GenomeAssembly}/vcf/${GenomeAssembly}.gnomad.af.only.vcf.gz
    Vcf_omni=/storage/references_and_index/${GenomeAssembly}/vcf/${GenomeAssembly}.omni.vcf.gz
    Vcf_hapmap=/storage/references_and_index/${GenomeAssembly}/vcf/${GenomeAssembly}.hapmap.vcf.gz
    WgsBqsrIntervals=/storage/references_and_index/${GenomeAssembly}/bed/${GenomeAssembly}.wgs.parallel.bqsr.interval_list
#---------------------------------------------------------------------------------------------------    

#---| RunLog |--------------------------------------------------------------------------------------
    RunLog=${LogDir}/${SeqID}.wgs.fastq-to-bam.post-merge.$(date '+%Y%m%d').log
    echo " " > ${RunLog}
    echo " WGS FASTQ-TO-BAM RUN LOG" >> ${RunLog}
    echo "----------------------------------------------------------------------" >> ${RunLog}
    echo " RUN DATE             $(date '+%Y-%m-%d %H:%M:%S')" >> ${RunLog}
    echo " BATCH ID             ${BatchID}" >> ${RunLog}
    echo " SEQ ID               ${SeqID}" >> ${RunLog}
    echo " REFERENCE VERSION    ${GenomeAssembly}" >> ${RunLog}
    echo "----------------------------------------------------------------------" >> ${RunLog}
    echo " " >> ${RunLog}   
#---------------------------------------------------------------------------------------------------

#---| MARK-DUPLICATES |-----------------------------------------------------------------------------

    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Mark-Duplicates -----------------------------------------###" >> ${RunLog}

    ${RunGatk4} --java-options "-Xmx32g" MarkDuplicates \
        --INPUT ${BamDir}/${SeqID}.merged.bam \
        --METRICS_FILE ${QcDir}/${SeqID}.mark.duplicates.metrics.txt \
        --OUTPUT ${BamDir}/${SeqID}.merged.dup.marked.bam \
        --CREATE_INDEX true \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
        --TAGGING_POLICY All \
        --MAX_RECORDS_IN_RAM 2000000 2>&1 | tee -a $RunLog

    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] MarkDuplicates Fnished." >> ${RunLog}
    echo " " >> ${RunLog}  
#---------------------------------------------------------------------------------------------------

#---| BASE RECALIBRATION : Base Recalibration Table |-----------------------------------------------

    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Base Recalibration --------------------------------------###" >> ${RunLog}

    BqsrTmp_IntervalsDir=${QcDir}/bqsr_intervals
    BqsrTmp_RecalTableDir=${QcDir}/recal_table_tmp
    BqsrTmp_BqsrBamDir=${BamDir}/bqsr_tmp
    mkdir -p ${BqsrTmp_BqsrBamDir}
    mkdir -p ${BqsrTmp_RecalTableDir} 
    mkdir -p ${BqsrTmp_IntervalsDir}

    echo "--- INTERVAL SPLIT ---" >> ${RunLog}
    ${RunGatk4} SplitIntervals \
        --reference ${GenomeFasta} \
        --intervals ${WgsBqsrIntervals} \
        --scatter-count ${Threads} \
        --output ${BqsrTmp_IntervalsDir} 2>&1 | tee -a $RunLog

    ls ${QcDir}/bqsr_intervals | grep "scattered.interval_list$" | sort > ${QcDir}/bqsr.shards.list   
    sed -i 's/.interval_list//g' ${QcDir}/bqsr.shards.list

    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Intervals Processing Fnished." >> ${RunLog}

    echo "--- RECALIBRATION TABLES CREATION ---" >> ${RunLog}
    parallel -j ${Threads} -k "${RunGatk4} BaseRecalibrator \
            --reference ${GenomeFasta} \
            --input ${BamDir}/${SeqID}.merged.dup.marked.bam \
            --known-sites ${Vcf_dbsnp} \
            --known-sites ${Vcf_hcsnp} \
            --known-sites ${Vcf_knwonindels} \
            --known-sites ${Vcf_gsindels} \
            --known-sites ${Vcf_gnomad} \
            --known-sites ${Vcf_omni} \
            --known-sites ${Vcf_hapmap} \
            --intervals ${BqsrTmp_IntervalsDir}/{}.interval_list \
            --output ${BqsrTmp_RecalTableDir}/{}.recal.table 
    " :::: ${QcDir}/bqsr.shards.list 
    
    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Recal-Tables Creation Parallel Processing Fnished." >> ${RunLog}
    echo "--- RECALIBRATION TABLES INTEGRATION ---" >> ${RunLog}

    if [ -f "${QcDir}/recal.table.args" ]; then rm ${QcDir}/recal.table.args; fi
    
    for LINE in $(seq 1 $(cat ${QcDir}/bqsr.shards.list | wc -l))
    do
        Shards=$(sed -n -e "${LINE} p" ${QcDir}/bqsr.shards.list | awk '{print $1}')
        echo -e "-I ${BqsrTmp_RecalTableDir}/${Shards}.recal.table" >> ${QcDir}/recal.table.args
    done

    echo -e "-O ${QcDir}/${SeqID}.recal.table" >> ${QcDir}/recal.table.args
    
    ${RunGatk4} GatherBQSRReports --arguments_file ${QcDir}/recal.table.args 2>&1 | tee -a $RunLog

    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Recal-Tables Merged." >> ${RunLog}
    echo "--- APPLY BQSR (PARALLEL PROCESS) ---" >> ${RunLog}
    
    parallel -j ${Threads} -k "${RunGatk4} ApplyBQSR \
        --reference ${GenomeFasta} \
        --input ${BamDir}/${SeqID}.merged.dup.marked.bam \
        --bqsr-recal-file ${QcDir}/${SeqID}.recal.table \
        --intervals ${BqsrTmp_IntervalsDir}/{}.interval_list \
        -O ${BqsrTmp_BqsrBamDir}/{}.bqsr.bam
    " :::: ${QcDir}/bqsr.shards.list 

    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] ApplyBQSR Parallel Processing Finished." >> ${RunLog}
    echo "--- MERGE BQSR BAM FILES ---" >> ${RunLog}
    
    find ${BqsrTmp_BqsrBamDir} -name "*.bqsr.bam" | sort > ${QcDir}/bqsr_bams.list
    samtools merge -@ ${Threads} -o ${BamDir}/${SeqID}.merged.dup.marked.recal.bam -b ${QcDir}/bqsr_bams.list

    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] BQSR-BAMs merged." >> ${RunLog}
    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Base Recalibration Finished." >> ${RunLog}
    echo " " >> ${RunLog}  
#---------------------------------------------------------------------------------------------------

#---| CONVERT Final BAM to ANALYSIS-READY-BAM |-----------------------------------------------------

    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] ANALYSIS-READY BAM CREATION -----------------------------###" >> ${RunLog}

    cp ${BamDir}/${SeqID}.merged.dup.marked.recal.bam ${BamDir}/${SeqID}.analysisReady.bam
    samtools index -@ ${Threads} ${BamDir}/${SeqID}.analysisReady.bam

    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] analysisReady-BAM Created." >> ${RunLog}
    echo " " >> ${RunLog}  
#---------------------------------------------------------------------------------------------------

#---| REMOVE TEMP-FILES |---------------------------------------------------------------------------

    # echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] REMOVE TEMP FILES ---------------------------------------###" >> ${RunLog}

    # if [ -f "${BamDir}/${SeqID}.analysisReady.bam" ]; then
    #     rm ${BamDir}/${SeqID}.merged.bam   
    #     rm ${BamDir}/${SeqID}.merged.bai 
    # # fi
    # # if [ -f "${BamDir}/${SeqID}.merged.dup.marked.bam" ]; then
    #     rm ${BamDir}/${SeqID}.merged.dup.marked.bam
    #     rm ${BamDir}/${SeqID}.merged.dup.marked.bai
    # # fi
    # # if [ -f "${BamDir}/${SeqID}.merged.dup.marked.bam" ]; then
    #     rm ${BamDir}/${SeqID}.merged.dup.marked.recal.bam
    # # fi
    # # if [ -f "${BamDir}/${SeqID}.bwa.mem.sam" ]; then 
    #     rm ${BamDir}/${SeqID}.bwa.mem.sam
    # # fi
    # # if [ -f "${BamDir}/${SeqID}.fastqtosam.bam" ]; then 
    #     rm ${BamDir}/${SeqID}.fastqtosam.bam
    # # fi
    # # if [ -d "${BamDir}/bqsr_tmp" ]; then 
    #     rm -r ${BamDir}/bqsr_tmp
    # fi

    # echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Temp Files Removed." >> ${RunLog}
    # echo " " >> ${RunLog}  
#---------------------------------------------------------------------------------------------------

#==================================================================================================# 
    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] WGS FASTQ-TO-BAM PROCESS FINISHED." >> ${RunLog}
    echo " " >> ${RunLog}
#==================================================================================================# 

