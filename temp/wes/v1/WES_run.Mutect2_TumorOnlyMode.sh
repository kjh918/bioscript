
#!/bin/bash

# SCRIPT VERSION : v0.1
# GUIDELINE ID   : GMC-NGS-B02-A09
# DATE           : 2024-03-11
# AUTHOR         : kangsm@gencurix.com

#---| Input Options |------------------------------------------------------------------------------#
    usage() {
        echo "Usage: Mutect2_TumorOnlyMode.sh [ ARGUMENTS ]...
        [ -d | --seqFolder ]      : data processing base directory
        [ -s | --seqID ]          : sample ID (Used as sample's processing product folder name)
        [ -h | --help ]           : Print this usage 
        [ --interval ]            : Mutect2 interval file (default = twist exome 2.0)
        [ --baseDir ]             : Base Work Dir
        [ --threads ]             : N-CPU threads
        [ --seqDepth ]            : Sequqncing Depth ( max depth if samples > 1 ) 
        [ --assembly ]            : Genome assembly version 
        [ --normalTonlyCall ]     : Run general Tumor-Only Calling
        [ --keepGermlineCall ]    : Run Keeping Germline Variants Calling (for purity & ploidy analysis)
        "
    }
#--------------------------------------------------------------------------------------------------#

#---| Set Defualt Values |-------------------------------------------------------------------------#
    THREADS=15
    BASE_DIR=/data/wes
    TARGET_INTERVAL=/storage/references_and_index/hg19/bed/Twist_WES_2.0/hg19_twist.exome.2.0.target.interval_list
    SEQ_DEPTH=200
    GENOME_ASSEMBLY=hg19
    NORMAL_CALL="true"
    KEEP_GERMLINE="false"
#--------------------------------------------------------------------------------------------------#


#---| Parse Option Arguments |---------------------------------------------------------------------#
    ARGS=$(getopt -a -o d:s:l:h: --long seqFolder:,seqID:,interval:,baseDir:,threads:,seqDepth:,assembly:,keepGermlineCall:,normalTonlyCall:,help -- "$@" )
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
        --interval )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = TwistExome2.0" ; shift ;;
                *)  TARGET_INTERVAL=$2 ; shift 2 ;;
            esac ;;
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
        --seqDepth )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 200" ; shift ;;
                *)  SEQ_DEPTH=$2 ; shift 2 ;;
            esac ;;
        --assembly )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = hg19" ; shift ;;
                *)  GENOME_ASSEMBLY=$2 ; shift 2 ;;
            esac ;;
        --normalTonlyCall )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = true" ; shift ;;
                *)  NORMAL_CALL=$2 ; shift 2 ;;
            esac ;;
        --keepGermlineCall )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = false" ; shift ;;
                *)  KEEP_GERMLINE=$2 ; shift 2 ;;
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
    BAM_DIR=${DATA_DIR}/${SEQ_ID}/bam
    VCF_DIR=${DATA_DIR}/${SEQ_ID}/vcf
    QC_DIR=${DATA_DIR}/${SEQ_ID}/qcfiles
    TMP_DIR=${DATA_DIR}/${SEQ_ID}/tmp
    LOG_DIR=${DATA_DIR}/${SEQ_ID}/log
    ##
    if [ ! -d ${VCF_DIR} ]; then mkdir -p ${VCF_DIR}; fi
    if [ ! -d ${TMP_DIR} ]; then mkdir -p ${TMP_DIR}; fi
    if [ ! -d ${LOG_DIR} ]; then mkdir -p ${LOG_DIR}; fi
#--------------------------------------------------------------------------------------------------#

#---| SWTOOLS RUN |--------------------------------------------------------------------------------#   
    RUN_SINGULARITY="singularity exec -B /storage,/data /storage/images"
    GATK4="${RUN_SINGULARITY}/gatk-4.4.0.0.sif gatk"
    SOBDetector="java -jar /storage/apps/SOBDetector-1.0.4/SOBDetector_v1.0.4.jar"
#--------------------------------------------------------------------------------------------------#

#---| LOGFILE |------------------------------------------------------------------------------------#
    LOGFILE=${LOG_DIR}/${SEQ_ID}.mutect2.tumor.only.mode.$(date '+%Y%m%d').log
    if [ ! -f ${LOGFILE} ]; then touch ${LOGFILE}; fi
#--------------------------------------------------------------------------------------------------#

#---| REFERENCES |---------------------------------------------------------------------------------#
    if [ "${GENOME_ASSEMBLY}" = "hg38" ]; then
        # GENOME FASTA #
        GENOME_FASTA=/storage/references_and_index/hg38/fasta/Homo_sapiens_assembly38.fasta
        # REFERENCE VCF #
        # SITES
        VCF_DBSNP=/storage/references_and_index/hg38/vcf/Homo_sapiens_assembly38.dbsnp138.vcf.gz
        VCF_1KG_SNP=/storage/references_and_index/hg38/vcf/1000G_phase1.snps.high_confidence.hg38.vcf.gz
        # INDEL
        VCF_INDEL=/storage/references_and_index/hg38/vcf/Homo_sapiens_assembly38.known_indels.vcf.gz
        VCF_1KG_MILLS=/storage/references_and_index/hg38/vcf/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
        # POPULATION
        VCF_GNOMAD=/storage/references_and_index/hg38/vcf/af-only-gnomad.hg38.vcf.gz
        VCF_PON=/storage/references_and_index/hg38/vcf/gatk_1kg_panel_of_normal_hg38.vcf.gz
    else
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
        # POPULATION
        VCF_GNOMAD=/storage/references_and_index/hg19/vcf/af-only-gnomad.raw.sites.grch37.vcf.gz
        VCF_PON=/storage/references_and_index/hg19/vcf/gatk_mutect2_wes_panel_of_normal_GRCh37.vcf.gz
    fi        
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo " " >> ${LOGFILE}
    echo "WES SomaticVariantCalling Mutect2 TumorOnlyMode log" >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " RUN DATE             : $(date '+%Y-%m-%d %H:%M:%S')" >> ${LOGFILE}
    echo " SEQ FOLDER           : ${DATA_DIR}" >> ${LOGFILE}
    echo " SEQ ID               : ${SEQ_ID}" >> ${LOGFILE}
    echo " GENOME ASSEMBLY      : ${GENOME_ASSEMBLY}" >> ${LOGFILE}
    echo " VARIANT CALLING MODE : Mutect2 Tumor-Only Mode" >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " " >> ${LOGFILE}
#==================================================================================================# 

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | VARIANT-CALLING START." >> ${LOGFILE}
#==================================================================================================# 

#---| TUMOR-ONLY VARIANT CALLING WITH MUTECT2 : GENERAL |------------------------------------------#
    F1R2_MAX_DEPTH=$(echo ${SEQ_DEPTH} | awk '{print $1*2}')  
    ## 
    if [ "${NORMAL_CALL}" = "true" ]; then

    #==============================================================================================# 
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MUTECT2 - General T-only Mode ]  START." >> ${LOGFILE}
    #==============================================================================================#

        ${GATK4} Mutect2 \
            --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads=14 -Xmx32768m" \
            -pairHMM AVX_LOGLESS_CACHING_OMP \
            --reference ${GENOME_FASTA} \
            --intervals ${TARGET_INTERVAL} \
            --input ${BAM_DIR}/${SEQ_ID}.analysisReady.bam \
            --germline-resource ${VCF_GNOMAD} \
            --panel-of-normals ${VCF_PON} \
            --f1r2-tar-gz ${QC_DIR}/${SEQ_ID}.f1r2.tar.gz \
            --f1r2-max-depth ${F1R2_MAX_DEPTH} \
            --min-base-quality-score 20 \
            --base-quality-score-threshold 20 \
            --output ${VCF_DIR}/${SEQ_ID}.mutect2.vcf

    #==============================================================================================# 
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MUTECT2 - General T-only Mode ]  FINISHED." >> ${LOGFILE}
    #==============================================================================================# 

    fi
#--------------------------------------------------------------------------------------------------#

#---| TUMOR-ONLY VARIANT CALLING WITH MUTECT2 : KEEP GERMLINE |------------------------------------#
    if [ "${KEEP_GERMLINE}" = "true" ]; then

    #==============================================================================================# 
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MUTECT2 - Keep Germline Variants Mode ]  START." >> ${LOGFILE}
    #==============================================================================================# 

        ${GATK4} Mutect2 \
            --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads=14 -Xmx32768m" \
            -pairHMM AVX_LOGLESS_CACHING_OMP \
            --reference ${GENOME_FASTA} \
            --intervals ${TARGET_INTERVAL} \
            --input ${BAM_DIR}/${SEQ_ID}.analysisReady.bam \
            --germline-resource ${VCF_GNOMAD} \
            --panel-of-normals ${VCF_PON} \
            --f1r2-tar-gz ${QC_DIR}/${SEQ_ID}.keep.germline.f1r2.tar.gz \
            --f1r2-max-depth ${F1R2_MAX_DEPTH} \
            --min-base-quality-score 20 \
            --base-quality-score-threshold 20 \
            --interval-padding 50 \
            --genotype-germline-sites true \
            --genotype-pon-sites true \
            --max-mnp-distance 0 \
            --output ${VCF_DIR}/${SEQ_ID}.mutect2.keep.germline.vcf

    #==============================================================================================# 
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MUTECT2 - Keep Germline Variants Mode ]  FINISHED." >> ${LOGFILE}
    #==============================================================================================# 

    fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | [ Read-Orientation-Model ] START." >> ${LOGFILE}
#==================================================================================================# 

#---| Orientation-Bias Filter |--------------------------------------------------------------------#
    if [ "${NORMAL_CALL}" = "true" ]; then

        ${GATK4} LearnReadOrientationModel \
            --java-options "-XX:ParallelGCThreads=14 -Xmx8192m" \
            -I ${QC_DIR}/${SEQ_ID}.f1r2.tar.gz \
            -O ${QC_DIR}/${SEQ_ID}.read-orientation-model.tar.gz

    fi
    #--------------------------------------------------------------------------#
    if [ "${KEEP_GERMLINE}" = "true" ]; then
        
        ${GATK4} LearnReadOrientationModel \
            --java-options "-XX:ParallelGCThreads=14 -Xmx8192m" \
            -I ${QC_DIR}/${SEQ_ID}.keep.germline.f1r2.tar.gz \
            -O ${QC_DIR}/${SEQ_ID}.keep.germline.read-orientation-model.tar.gz
    
    fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | [ Read-Orientation-Model ]  FINISHED." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | [ Pileup Summary ] START." >> ${LOGFILE}
#==================================================================================================# 

#---| Generate Pileup Summaries |------------------------------------------------------------------#
    PILEUP_TABLE=${QC_DIR}/${SEQ_ID}.targeted_sequencing.table
    #--------------------------------------------------------------------------#
    if [ ! -f ${PILEUP_TABLE} ]; then
        ${GATK4} GetPileupSummaries \
            --java-options "-XX:ParallelGCThreads=14 -Xmx8192m" \
            -I ${BAM_DIR}/${SEQ_ID}.analysisReady.bam \
            -O ${PILEUP_TABLE} \
            -V ${VCF_GNOMAD} \
            -L ${TARGET_INTERVAL} \
            -R ${GENOME_FASTA}
    fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | [ Pileup Summary ] FINISHED." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | [ Contamination-Calculation ] START." >> ${LOGFILE}
#==================================================================================================# 

#---| Calculate Contamination |--------------------------------------------------------------------#
    ${GATK4} CalculateContamination \
        --java-options "-XX:ParallelGCThreads=14 -Xmx16384m" \
        -I ${QC_DIR}/${SEQ_ID}.targeted_sequencing.table \
        -O ${QC_DIR}/${SEQ_ID}.contamination.table
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | [ Contamination-Calculation ]  FINISHED." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | [ FILTER-MUTECT2 ] START." >> ${LOGFILE}
#==================================================================================================# 

#---| Run FilterMutectCalls |----------------------------------------------------------------------#
    if [ "${NORMAL_CALL}" = "true" ]; then
        ${GATK4} FilterMutectCalls \
            --java-options "-XX:ParallelGCThreads=14 -Xmx16384m" \
            -V ${VCF_DIR}/${SEQ_ID}.mutect2.vcf \
            -L ${TARGET_INTERVAL} \
            --reference ${GENOME_FASTA} \
            --contamination-table ${QC_DIR}/${SEQ_ID}.contamination.table \
            --ob-priors ${QC_DIR}/${SEQ_ID}.read-orientation-model.tar.gz \
            -O ${VCF_DIR}/${SEQ_ID}.mutect2.filtered.vcf
    fi
    #----------------------------------------------------------------------------------------------#
    if [ "${KEEP_GERMLINE}" = "true" ]; then
        ${GATK4} FilterMutectCalls \
            --java-options "-XX:ParallelGCThreads=14 -Xmx16384m" \
            -V ${VCF_DIR}/${SEQ_ID}.mutect2.keep.germline.vcf \
            -L ${TARGET_INTERVAL} \
            --reference ${GENOME_FASTA} \
            --contamination-table ${QC_DIR}/${SEQ_ID}.contamination.table \
            --ob-priors ${QC_DIR}/${SEQ_ID}.keep.germline.read-orientation-model.tar.gz \
            -O ${VCF_DIR}/${SEQ_ID}.mutect2.keep.germline.filtered.vcf
    fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | [ FILTER-MUTECT2 ]  FINISHED." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | [ SOBEdetector ] START." >> ${LOGFILE}
#==================================================================================================# 

#---| Run SOBDector |------------------------------------------------------------------------------#
    if [ "${NORMAL_CALL}" = "true" ]; then
        ${SOBDetector} --input-type VCF \
            --input-variants ${VCF_DIR}/${SEQ_ID}.mutect2.filtered.vcf \
            --input-bam ${BAM_DIR}/${SEQ_ID}.analysisReady.bam \
            --output-variants ${VCF_DIR}/${SEQ_ID}.mutect2.bias.filtered.vcf \
            --minBaseQuality 20 \
            --minMappingQuality 20 \
            --only-passed false
    fi
    #----------------------------------------------------------------------------------------------#
    if [ "${KEEP_GERMLINE}" = "true" ]; then
        ${SOBDetector} --input-type VCF \
            --input-variants ${VCF_DIR}/${SEQ_ID}.mutect2.keep.germline.filtered.vcf \
            --input-bam ${BAM_DIR}/${SEQ_ID}.analysisReady.bam \
            --output-variants ${VCF_DIR}/${SEQ_ID}.mutect2.keep.germline.bias.filtered.vcf \
            --minBaseQuality 20 \
            --minMappingQuality 20 \
            --only-passed false
    fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | [ SOBEdetector ]  FINISHED." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | VARIANT-CALLING FINISHED." >> ${LOGFILE}
#==================================================================================================# 

#==================================================================================================# 
    echo " " >> ${LOGFILE}
    echo "WES SomaticVariantCalling Mutect2 TumorOnlyMode DONE." >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " " >> ${LOGFILE}
#==================================================================================================# 

