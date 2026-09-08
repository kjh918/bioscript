
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

        --RGSM              [ REQUIRED ] read-group-sample-id
                            default = 'RGID'
                        
        --RGLB              read-group-library
                            default = WGS
        
        --RGPL              read-group-platform 
                            default = ILLUMINA
        
        --RGCN              read-group-center
                            default = GCX
        
        --RUN_HMFTOOLS_BAM  run bwa-mem2 alignment for hmftools analysis
                            default = false

        -h | --help         Print this usage 
        "
    }
#---------------------------------------------------------------------------------------------------

#---| Parse Option Arguments |----------------------------------------------------------------------
    ARGS=$(getopt -a -o h: --long BATCH_ID:,SEQ_ID:,BASE_DIR:,THREADS:,GENOME_ASSEMBLY:,RGID:,RGSM:,RGLB:,RGPL:,RGCN:,RUN_HMFTOOLS_BAM:,help -- "$@" )
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
        --RUN_HMFTOOLS_BAM )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = false" ; shift ;;
                *)  RunHmftoolsBam=$2 ; shift 2 ;;
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
    if [ -z "${RGSM}" ]; then exit 2; fi
    if [ -z "${Threads}" ]; then Threads=15; fi
    if [ -z "${GenomeAssembly}" ]; then GenomeAssembly="hg38"; fi
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
    WgsBqsrIntervals=/storage/references_and_index/${GenomeAssembly}/bed/${GenomeAssembly}.wgs.bqsr.interval_list
#---------------------------------------------------------------------------------------------------    

#---| RunLog |--------------------------------------------------------------------------------------
    RunLog=${LogDir}/${SeqID}.wgs.fastq-to-bam.pre-merge.$(date '+%Y%m%d').log
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
        --READ_GROUP_NAME ${SeqID} \
        --PLATFORM ${RGPL} \
        --LIBRARY_NAME ${RGLB} \
        --SAMPLE_NAME ${RGSM} \
        --SEQUENCING_CENTER ${RGCN} \
        --TMP_DIR ${TmpDir} 2>&1 | tee -a $RunLog

    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] FastqToSam Finished." >> ${RunLog}
    echo " " >> ${RunLog}  
#---------------------------------------------------------------------------------------------------

#---| aBAM (Aligned BAM) GENERATION |---------------------------------------------------------------
    
    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] BWA-BAM Generation --------------------------------------###" >> ${RunLog}

    ${RunBwa} -M -Y -q -K 100000000 -t ${Threads} \
        ${BwaIndex} \
        ${FastqDir}/${SeqID}.trimmed_R1.fastq.gz \
        ${FastqDir}/${SeqID}.trimmed_R2.fastq.gz > ${BamDir}/${SeqID}.bwa.mem.sam 
    
    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] BWA-MEM Finished." >> ${RunLog}
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

    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] MergeBamAlignment Finished." >> ${RunLog}

    rm ${BamDir}/${SeqID}.fastqtosam.bam
    rm ${BamDir}/${SeqID}.bwa.mem.sam

    echo " " >> ${RunLog}  
#---------------------------------------------------------------------------------------------------

#---| HMFTOOLS BWA-MEM2 ALIGNMENT |-----------------------------------------------------------------

    # if [ "${RunHmftoolsBam}" = "true" ]; then

    #     echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] RUN BWA-MEM2 FOR HMFTOOLS -------------------------------###" >> ${RunLog}

    #     RunBwaMem2="${RunSingularity}/bwa-mem2-v2.2.1.sif bwa-mem2"

    #     if [ "${GenomeAssembly}" = "hg19" ]; then
    #         HMF_GenomeFasta="/storage/references_and_index/hg19/hmftools/Homo_sapiens.GRCh37.GATK.illumina.fasta"
    #         HMF_BwaMem2Index="/storage/references_and_index/hg19/hmftools/bwa-mem2_index-2.2.1/Homo_sapiens.GRCh37.GATK.illumina.fasta"
    #     else
    #         HMF_GenomeFasta="/storage/references_and_index/hg38/hmftools/GRCh38_masked_exclusions_alts_hlas.fasta"
    #         HMF_BwaMem2Index="/storage/references_and_index/hg38/hmftools/bwa-mem2_index-2.2.1/GRCh38_masked_exclusions_alts_hlas.fasta"
    #     fi
    
    #     case "${Threads}" in
    #         16) SortThreads=4 ;;
    #         24) SortThreads=6 ;;
    #         32) SortThreads=8 ;;
    #     esac
    #     BwaThreads=$(( Threads - SortThreads ))

    #     ${RunBwaMem2} mem -M -Y -q -K 100000000 -t ${BwaThreads} \
    #         -R "@RG\tID:${RGID}\tSM:${RGSM}\tLB:WGS\tPL:ILLUMINA\tCN:GCX" \
    #         ${HMF_BwaMem2Index} \
    #         ${FastqDir}/${SeqID}.trimmed_R1.fastq.gz \
    #         ${FastqDir}/${SeqID}.trimmed_R2.fastq.gz \
    #         | samtools sort -@ ${SortThreads} -o ${BamDir}/${SeqID}.hmftools.bwa.mem2.bam - 

    #     echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] BWA-MEM2 Finished." >> ${RunLog}

    #     samtools index -@ ${Threads} ${BamDir}/${SeqID}.hmftools.bwa.mem2.bam 

    #     echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Samtools Index Finished." >> ${RunLog}
    #     echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] BWA-MEM2 FOR HMFTOOLS Finished." >> ${RunLog}
    #     echo " " >> ${RunLog}  
    # fi
#---------------------------------------------------------------------------------------------------

#==================================================================================================# 
    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] WGS FASTQ-TO-BAM PROCESS FINISHED." >> ${RunLog}
    echo " " >> ${RunLog}
#==================================================================================================# 
