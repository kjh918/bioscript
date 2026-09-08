
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
        
        --RGID              [ REQURIED ] read-group-id

        --RGSM              read-group-sample-id
                            default = 'RGID'
                        
        --RGLB              read-group-library
                            default = WGS
        
        --RGPL              read-group-platform 
                            default = ILLUMINA
        
        --RGCN              read-group-center
                            default = GCX
        
        -h | --help         Print this usage 
        "
    }
#---------------------------------------------------------------------------------------------------

#---| Parse Option Arguments |----------------------------------------------------------------------
    ARGS=$(getopt -a -o h: --long BATCH_ID:,SEQ_ID:,BASE_DIR:,THREADS:,GENOME_ASSEMBLY:,RGID:,RGSM:,RGLB:,RGPL:,RGCN:,help -- "$@" )
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
        --RGID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  RGID=$2 ; shift 2 ;;
            esac ;;
        --RGSM )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 'RGID'" ; shift ;;
                *)  RGSM=$2 ; shift 2 ;;
            esac ;;
        --RGLB )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = WGS" ; shift ;;
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
#---------------------------------------------------------------------------------------------------

#---| Check Params and Set Default -----------------------------------------------------------------
    if [ -z "${BaseDir}" ]; then exit 2; fi
    if [ -z "${BatchID}" ]; then exit 2; fi
    if [ -z "${SeqID}" ]; then exit 2; fi
    if [ -z "${RGID}" ]; then exit 2; fi
    if [ -z "${Threads}" ]; then Threads=15; fi
    if [ -z "${GenomeAssembly}" ]; then GenomeAssembly="hg38"; fi
    if [ -z "${RGSM}" ]; then RGSM=${RGID}; fi
    if [ -z "${RGLB}" ]; then RGLB="WGS"; fi
    if [ -z "${RGPL}" ]; then RGPL="ILLUMINA"; fi
    if [ -z "${RGCN}" ]; then RGCN="GCX"; fi
#---------------------------------------------------------------------------------------------------

#---| Manual Params |-------------------------------------------------------------------------------
    # BaseDir=/data/wgs
    # BatchID=cancer_cell_lines
    # SeqID=
    # Threads=15
    # GenomeAssembly="hg38"
    # RGID=SRR22509543
    # RGSM=
    # RGPL="NGS"
    # RGLB="WGS"
    # RGCN="GCX"
#---------------------------------------------------------------------------------------------------

#---| FOLDERS |-------------------------------------------------------------------------------------
    DataDir=${BaseDir}/${BatchID}
    ##  
    FastqDir=${DataDir}/${SeqID}/fastq
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
    WgsCallTargetIntervals=/storage/references_and_index/${GenomeAssembly}/bed/${GenomeAssembly}.wgs.call.target.region.interval_list
#---------------------------------------------------------------------------------------------------    

#---| RunLog |--------------------------------------------------------------------------------------
    RunLog=${LogDir}/${SeqID}.wgs.fastq-to-bam.$(date '+%Y%m%d').log
    echo " " > ${RunLog}
    echo " WGS FASTQ-TO-BAM RUN LOG" >> ${RunLog}
    echo "----------------------------------------------------------------------" >> ${RunLog}
    echo " RUN DATE             $(date '+%Y-%m-%d %H:%M:%S')" >> ${RunLog}
    echo " BATCH ID             ${BatchID}" >> ${RunLog}
    echo " SEQ ID               ${SeqID}" >> ${RunLog}
    echo " READ GROUP ID        ${RGID}" >> ${RunLog}
    echo " SAMPLE NAME          ${RGSM}" >> ${RunLog}
    echo " REFERENCE VERSION    ${GenomeAssembly}" >> ${RunLog}
    echo "----------------------------------------------------------------------" >> ${RunLog}
    echo " " >> ${RunLog}   
#---------------------------------------------------------------------------------------------------

#---| uBMA (Unligned BAM) GENERATION |--------------------------------------------------------------
    
    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] uBAM Generation -----------------------------------------###" >> ${RunLog}

    ${RunGatk4} --java-options "-Xmx8g" FastqToSam \
        --FASTQ ${FastqDir}/${SeqID}.trimmed_R1.fastq.gz \
        --FASTQ2 ${FastqDir}/${SeqID}.trimmed_R2.fastq.gz \
        --OUTPUT ${BamDir}/${SeqID}.fastqtosam.bam \
        --QUALITY_FORMAT Standard \
        --READ_GROUP_NAME ${RGID} \
        --PLATFORM ${RGPL} \
        --LIBRARY_NAME ${RGLB} \
        --SAMPLE_NAME ${RGSM} \
        --SEQUENCING_CENTER ${RGCN} \
        --TMP_DIR ${TmpDir} 2>&1 | tee -a $RunLog

    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] FastqToSam Fnished." >> ${RunLog}
    echo " " >> ${RunLog}  
#---------------------------------------------------------------------------------------------------

#---| aBAM (Aligned BAM) GENERATION |---------------------------------------------------------------
    
    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] BWA-BAM Generation --------------------------------------###" >> ${RunLog}

    ${RunBwa} -M -Y -q -K 100000000 -t ${Threads} \
        ${BwaIndex} \
        ${FastqDir}/${SeqID}.trimmed_R1.fastq.gz \
        ${FastqDir}/${SeqID}.trimmed_R2.fastq.gz > ${BamDir}/${SeqID}.bwa.mem.sam 

    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] BWA-MEM Fnished." >> ${RunLog}
    echo " " >> ${RunLog}  
#---------------------------------------------------------------------------------------------------

#---| BAM MERGE ALIGNMENT |-------------------------------------------------------------------------

    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Merge BAM Alignment -------------------------------------###" >> ${RunLog}

    ${RunGatk4} --java-options "-Xmx24g" MergeBamAlignment \
        --REFERENCE_SEQUENCE ${GenomeFasta} \
        --UNMAPPED_BAM ${BamDir}/${SeqID}.fastqtosam.bam \
        --ALIGNED_BAM ${BamDir}/${SeqID}.bwa.mem.sam \
        --OUTPUT ${BamDir}/${SeqID}.merged.bam \
        --EXPECTED_ORIENTATIONS FR \
        --PRIMARY_ALIGNMENT_STRATEGY BestMapq \
        --SORT_ORDER coordinate \
        --CREATE_INDEX true \
        --INCLUDE_SECONDARY_ALIGNMENTS true \
        --ALIGNER_PROPER_PAIR_FLAGS false \
        --CLIP_ADAPTERS true \
        --ADD_MATE_CIGAR true \
        --MAX_INSERTIONS_OR_DELETIONS 1 \
        --VALIDATION_STRINGENCY LENIENT \
        --TMP_DIR ${TmpDir} 2>&1 | tee -a $RunLog

    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] MergeBamAlignment Fnished." >> ${RunLog}
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

    # direct way
    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Create recal-table for WGS-Call-Target-Intervals Start." >> ${RunLog}

    ${RunGatk4} BaseRecalibrator \
        --reference ${GenomeFasta} \
        --input ${BamDir}/${SeqID}.merged.dup.marked.bam \
        --known-sites ${Vcf_dbsnp} \
        --known-sites ${Vcf_hcsnp} \
        --known-sites ${Vcf_knwonindels} \
        --known-sites ${Vcf_gsindels} \
        --known-sites ${Vcf_gnomad} \
        --known-sites ${Vcf_omni} \
        --known-sites ${Vcf_hapmap} \
        --intervals ${WgsCallTargetIntervals} \
        --output ${QcDir}/${SeqID}.recal.table 2>&1 | tee -a $RunLog

    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Create recal-table for WGS-Call-Target-Intervals Finished." >> ${RunLog}
    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Apply recal-table to Bam Start." >> ${RunLog}

    ${RunGatk4} ApplyBQSR \
        --reference ${GenomeFasta} \
        --input ${BamDir}/${SeqID}.merged.dup.marked.bam \
        --bqsr-recal-file ${QcDir}/${SeqID}.recal.table \
        --intervals ${WgsCallTargetIntervals} \
        -O ${BamDir}/${SeqID}.merged.dup.marked.recal.bam 2>&1 | tee -a $RunLog

    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Apply recal-table to Bam Finished." >> ${RunLog}

    # parallel way
    # BqsrTmp_IntervalsDir=${QcDir}/bqsr_intervals
    # BqsrTmp_RecalTableDir=${QcDir}/recal_table_tmp
    # BqsrTmp_BqsrBamDir=${BamDir}/bqsr_tmp
    # mkdir -p ${BqsrTmp_BqsrBamDir} ${BqsrTmp_RecalTableDir} ${BqsrTmp_IntervalsDir}

    # echo "--- INTERVAL SPLIT ---" >> ${RunLog}
    # ${RunGatk4} SplitIntervals \
    #     --reference ${GenomeFasta} \
    #     --intervals ${WgsBqsrIntervals} \
    #     --scatter-count ${Threads} \
    #     --output ${BqsrTmp_IntervalsDir} 2>&1 | tee -a $RunLog

    # ls ${QcDir}/bqsr_intervals | grep "scattered.interval_list$" | sort > ${QcDir}/bqsr.shards.list   
    # sed -i 's/.interval_list//g' ${QcDir}/bqsr.shards.list

    # echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Intervals Processing Fnished." >> ${RunLog}

    # echo "--- RECALIBRATION TABLES CREATION ---" >> ${RunLog}
    # parallel -j ${Threads} -k "${RunGatk4} BaseRecalibrator \
    #         --reference ${GenomeFasta} \
    #         --input ${BamDir}/${SeqID}.merged.dup.marked.bam \
    #         --known-sites ${Vcf_dbsnp} \
    #         --known-sites ${Vcf_hcsnp} \
    #         --known-sites ${Vcf_knwonindels} \
    #         --known-sites ${Vcf_gsindels} \
    #         --known-sites ${Vcf_gnomad} \
    #         --known-sites ${Vcf_omni} \
    #         --known-sites ${Vcf_hapmap} \
    #         --intervals ${BqsrTmp_IntervalsDir}/{}.interval_list \
    #         --output ${BqsrTmp_RecalTableDir}/{}.recal.table 
    # " :::: ${QcDir}/bqsr.shards.list 
    
    # echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Recal-Tables Creation Parallel Processing Fnished." >> ${RunLog}
    # echo "--- RECALIBRATION TABLES INTEGRATION ---" >> ${RunLog}

    # if [ -f "${QcDir}/recal.table.args" ]; then rm ${QcDir}/recal.table.args; fi
    
    # for LINE in $(seq 1 $(cat ${QcDir}/bqsr.shards.list | wc -l))
    # do
    #     Shards=$(sed -n -e "${LINE} p" ${QcDir}/bqsr.shards.list | awk '{print $1}')
    #     echo -e "-I ${BqsrTmp_RecalTableDir}/${Shards}.recal.table" >> ${QcDir}/recal.table.args
    # done

    # echo -e "-O ${QcDir}/${SeqID}.recal.table" >> ${QcDir}/recal.table.args
    
    # ${RunGatk4} GatherBQSRReports --arguments_file ${QcDir}/recal.table.args 2>&1 | tee -a $RunLog

    # echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Recal-Tables Merged." >> ${RunLog}
    # echo "--- APPLY BQSR (PARALLEL PROCESS) ---" >> ${RunLog}
    
    # parallel -j ${Threads} -k "${RunGatk4} ApplyBQSR \
    #     --reference ${GenomeFasta} \
    #     --input ${BamDir}/${SeqID}.merged.dup.marked.bam \
    #     --bqsr-recal-file ${QcDir}/${SeqID}.recal.table \
    #     --intervals ${BqsrTmp_IntervalsDir}/{}.interval_list \
    #     -O ${BqsrTmp_BqsrBamDir}/{}.bqsr.bam
    # " :::: ${QcDir}/bqsr.shards.list 

    # echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] ApplyBQSR Parallel Processing Finished." >> ${RunLog}
    # echo "--- MERGE BQSR BAM FILES ---" >> ${RunLog}
    
    # find ${BqsrTmp_BqsrBamDir} -name "*.bqsr.bam" | sort > ${QcDir}/bqsr_bams.list
    # samtools merge -@ ${Threads} -o ${BamDir}/${SeqID}.merged.dup.marked.recal.bam -b ${QcDir}/bqsr_bams.list

    # echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] BQSR-BAMs merged." >> ${RunLog}
    # echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Base Recalibration Finished." >> ${RunLog}
    # echo " " >> ${RunLog}  
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

    # rm ${BamDir}/${SeqID}.merged.bam   
    # rm ${BamDir}/${SeqID}.merged.bai 
    # rm ${BamDir}/${SeqID}.merged.dup.marked.bam
    # rm ${BamDir}/${SeqID}.merged.dup.marked.bai
    # rm ${BamDir}/${SeqID}.merged.dup.marked.recal.bam
    # rm ${BamDir}/${SeqID}.bwa.mem.sam
    # rm ${BamDir}/${SeqID}.fastqtosam.bam
    # rm -r ${BamDir}/bqsr_tmp

    # echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Temp Files Removed." >> ${RunLog}
    # echo " " >> ${RunLog}  
#---------------------------------------------------------------------------------------------------

#==================================================================================================# 
    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] WGS FASTQ-TO-BAM PROCESS FINISHED." >> ${RunLog}
    echo " " >> ${RunLog}
#==================================================================================================# 

