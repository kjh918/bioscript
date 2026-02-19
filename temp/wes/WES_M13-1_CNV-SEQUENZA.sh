
#! /bin/bash 

#---| Input Options |------------------------------------------------------------------------------#
    usage() {
        echo "Usage: WES_run.CNV_Sequenza.sh [ ARGUMENTS ]...
        [ --baseDir     ] : base work dir 
        [ --seqFolder   ] : seq folder dir 
        [ --tumorSeqID  ] : tumoe sample seqID
        [ --normalSeqID ] : normal sample seqID
        [ --assembly    ] : reference genome. default = 'hg19'
        [ --wesLibKit   ] : WES library kit. default = 'twist.exome.2.0'
        [ --threads     ] : N-cpu threads. default = 15
        [ --memLimit    ] : memory usage limit (GB). default = 25 
        [ --binSize     ] : binning size. default = 50 
        [ -h | --help   ] : Print this usage 
        "
    }
#--------------------------------------------------------------------------------------------------#

#---| Set Defualt Values |-------------------------------------------------------------------------#
    THREADS=15
    MEM_LIMIT=20
    BIN_SIZE=100
    GENOME_ASSEMBLY="hg19"
    WES_LIB_KIT="twist.exome.2.0"
#--------------------------------------------------------------------------------------------------#

#---| Parse Option Arguments |---------------------------------------------------------------------#
    ARGS=$(getopt -a -o h: --long baseDir:,seqFolder:,tumorSeqID:,normalSeqID:,assembly:,wesLibKit:,threads:,memLimit:,binSize:,help -- "$@" )
    VALID_ARGS=$?
    if [ "$VALID_ARGS" != "0" ]; then 
        usage >&2 
        exit 2
    fi

    eval set -- "$ARGS"
    while :
    do
    case "$1" in
        --baseDir )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = /data/wes " ; shift ;;
                *)  BASE_DIR=$2 ; shift 2 ;;
            esac ;;
        --seqFolder )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = /data/wes " ; shift ;;
                *)  SEQ_FOLDER=$2 ; shift 2 ;;
            esac ;;
        --tumorSeqID )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 10" ; shift ;;
                *)  TUMOR_SEQ_ID=$2 ; shift 2 ;;
            esac ;;
        --normalSeqID )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 20" ; shift ;;
                *)  NORMAL_SEQ_ID=$2 ; shift 2 ;;
            esac ;;
        --assembly )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 10" ; shift ;;
                *)  GENOME_ASSEMBLY=$2 ; shift 2 ;;
            esac ;;
        --wesLibKit )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 10" ; shift ;;
                *)  WES_LIB_KIT=$2 ; shift 2 ;;
            esac ;;
        --threads )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 10" ; shift ;;
                *)  THREADS=$2 ; shift 2 ;;
            esac ;;
        --memLimit )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 10" ; shift ;;
                *)  MEM_LIMIT=$2 ; shift 2 ;;
            esac ;;
        --binSize )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 10" ; shift ;;
                *)  BIN_SIZE=$2 ; shift 2 ;;
            esac ;;
        -h | --help )
        usage >&2 ; exit 2 ;;   
        --) shift ; break ;;
        *)  usage >&2 ; exit 2 ;;
    esac
    done
#--------------------------------------------------------------------------------------------------#

#---| GC WIG : PRE-PROCESS DATA |------------------------------------------------------------------#
    # GENOME_FASTA=/storage/references_and_index/hg19/fasta/human_g1k_v37_decoy.fasta
    # GENOME_FASTA=/storage/references_and_index/hg38/fasta/Homo_sapiens_assembly38.fasta
    # OUT_WIG=/storage/home/kangsm/myDB/cnv_sequenza/hg38_GC.50.base.wig.gz
    # singularity exec -B /storage,/data /storage/images/sequenza.sif sequenza-utils gc_wiggle -f ${GENOME_FASTA} -o ${OUT_WIG} -w 50
#--------------------------------------------------------------------------------------------------#
  
#---| RUN SEQUENZA |-------------------------------------------------------------------------------#
    SEQZ_DIR=${BASE_DIR}/${SEQ_FOLDER}/${TUMOR_SEQ_ID}/cnv/sequenza
    if [ ! -d ${SEQZ_DIR} ]; then mkdir -p ${SEQZ_DIR}; fi
    #--------------------------------------------------------------------------#
    TMP_DIR=${SEQZ_DIR}/tmp
    DB_DIR=${SEQZ_DIR}/database
    if [ ! -d ${TMP_DIR} ]; then mkdir -p ${TMP_DIR}; fi
    if [ ! -d ${DB_DIR} ]; then mkdir -p ${DB_DIR}; fi
    #--------------------------------------------------------------------------#
    cd ${SEQZ_DIR}
    #--------------------------------------------------------------------------#
    TUMOR_BAM=${BASE_DIR}/${SEQ_FOLDER}/${TUMOR_SEQ_ID}/bam/${TUMOR_SEQ_ID}.analysisReady.bam
    TUMOR_BAM_INDEX=${BASE_DIR}/${SEQ_FOLDER}/${TUMOR_SEQ_ID}/bam/${TUMOR_SEQ_ID}.analysisReady.bam.bai
    NORMAL_BAM=${BASE_DIR}/${SEQ_FOLDER}/${NORMAL_SEQ_ID}/bam/${NORMAL_SEQ_ID}.analysisReady.bam
    NORMAL_BAM_INDEX=${BASE_DIR}/${SEQ_FOLDER}/${NORMAL_SEQ_ID}/bam/${NORMAL_SEQ_ID}.analysisReady.bam.bai
    if [ ! -f ${TUMOR_BAM} ]; then
        echo "Tumor BAM File NOT FOUMD. STOP."
        exit 2
    fi
    if [ ! -f ${NORMAL_BAM} ]; then
        echo "Normal BAM File NOT FOUMD. STOP."
        exit 2
    fi
    #--------------------------------------------------------------------------#
    if [ "${GENOME_ASSEMBLY}" = "hg38" ]; then
        GENOME_FASTA_GZ=/storage/references_and_index/hg38/fasta/Homo_sapiens_assembly38.fasta.gz
        GC_WIG=/storage/home/kangsm/myDB/cnv_sequenza/hg38_GC.${BIN_SIZE}.base.wig.gz
        BED_FILE=/storage/references_and_index/hg38/bed/${WES_LIB_KIT}/hg38_${WES_LIB_KIT}.target.bed
    else
        GENOME_FASTA_GZ=/storage/references_and_index/hg19/fasta/human_g1k_v37_decoy.fasta.gz
        GC_WIG=/storage/home/kangsm/myDB/cnv_sequenza/hg19_GC.${BIN_SIZE}.base.wig.gz
        BED_FILE=/storage/references_and_index/hg19/bed/${WES_LIB_KIT}/hg19_${WES_LIB_KIT}.target.bed
    fi
    #--------------------------------------------------------------------------#
    singularity exec -B /storage,/data -B ${DB_DIR}:/databases -B ${TMP_DIR}:/tmp /storage/images/sequenza.sif sequenza-pipeline \
        --sample-id ${TUMOR_SEQ_ID} \
        --normal-bam ${NORMAL_BAM} --normal-bam-index ${NORMAL_BAM_INDEX} \
        --tumor-bam ${TUMOR_BAM} --tumor-bam-index ${TUMOR_BAM_INDEX} \
        --gc_wig ${GC_WIG} --reference-gz ${GENOME_FASTA_GZ} \
        --breaks ${BED_FILE} --bin ${BIN_SIZE} \
        --mem ${MEM_LIMIT} --ncpu ${THREADS} --no_archive
    #--------------------------------------------------------------------------#
    cp ${TMP_DIR}/workdir/seqz/${TUMOR_SEQ_ID}_bin${BIN_SIZE}.seqz.gz ${SEQZ_DIR}/
    cp ${TMP_DIR}/workdir/seqz/${TUMOR_SEQ_ID}_bin${BIN_SIZE}.seqz.gz.tbi ${SEQZ_DIR}/
    #--------------------------------------------------------------------------#
    Rscript /storage/home/kangsm/runScripts/WES_run.CNV_SequenzaAnalysis.R \
        --BASE_DIR ${BASE_DIR} \
        --SEQ_FOLDER ${SEQ_FOLDER} \
        --SEQ_ID ${TUMOR_SEQ_ID} \
        --RESULT_DIR ${SEQZ_DIR} \
        --BIN_SIZE ${BIN_SIZE}
#--------------------------------------------------------------------------------------------------#





    


