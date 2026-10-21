
#! /bin/bash

# script version v0.1

#---| Input Options |--------------------------------------------------------------------------------------------------------------------------------#
    usage() {
        echo "Usage: GCX-WES-TumorOnlyVariantCalling.sh [ options ]
        [ --tumor ]         : input tumor bam file (with absolute path).
        [ --normal ]        : input matched normal bam file (with absolute path, required if runMode = 'paired').
        [ --ponDir ]        : panel of normal bam files folder path (required if runMode = 'pon').
        [ --sampleId ]      : sample ID (used as result folder name).
        [ --runMode ]       : analysis mode 'paired' or 'pon'. no default.
        [ --runDir ]        : run folder. all files except results will be made in this folder.
        [ --resultDir ]     : result folder.
        [ --bed ]           : target region bed. 'sureselect','sureselect_utr','twist'. default = 'twist'
        [ --genomeVersion ] : reference genome version. 'hg19' or 'hg38'.
        [ -h | --help ]     : Print this usage 
        "
    }
    #---| Set Defualt Values |-------------------------------------------------#
        PON_BAM_DIR=""
        NORMAL_BAM=""
        WES_LIB="twistexome2"
        CNV_RUN_DIR=""
#----------------------------------------------------------------------------------------------------------------------------------------------------#

#---| Parse Option Arguments |-----------------------------------------------------------------------------------------------------------------------#

    ARGS=$(getopt -a -o h: --long tumor:,normal:,ponDir:,sampleId:,runMode:,runDir:,resultDir:,bed:,genomeVersion:,help: -- "$@" )
    VALID_ARGS=$?
    if [ "$VALID_ARGS" != "0" ]; then 
        usage >&2 
        exit 2
    fi

    eval set -- "$ARGS"
    while :
    do
    case "$1" in
        --tumor )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  TUMOR_BAM=$2 ; shift 2 ;;
            esac ;;
        --normal )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = NULL" ; shift ;;
                *)  NORMAL_BAM=$2 ; shift 2 ;;
            esac  ;;
        --ponDir )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = NULL" ; shift ;;
                *)  PON_BAM_DIR=$2 ; shift 2 ;;
            esac ;;
        --sampleId )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  SAMPLE_ID=$2 ; shift 2 ;;
            esac ;;
        --runMode )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  RUN_MODE=$2 ; shift 2 ;;
            esac ;;
        --runDir )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = resultDir" ; shift ;;
                *)  CNV_RUN_DIR=$2 ; shift 2 ;;
            esac ;;
        --resultDir )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  CNV_RES_DIR=$2 ; shift 2 ;;
            esac ;;
        --bed )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  WES_LIB=$2 ; shift 2 ;;
            esac ;;
        --genomeVersion )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  REF_VERSION=$2 ; shift 2 ;;
            esac ;;
        -h | --help )
        usage >&2 ; exit 2 ;;   
        --) shift ; break ;;
        *)  usage >&2 ; exit 2 ;;
    esac
    done
#----------------------------------------------------------------------------------------------------------------------------------------------------#

#---| REQUIRED FILED FOR ANALYSIS RUN |--------------------------------------------------------------------------------------------------------------#

CNV_FILES_DIR=/storage/home/kangsm/myDB/cnv
REFERENCE_DIR=/storage/references_and_index

if [ "${REF_VERSION}" = "hg38" ]; then
    GENOME_FASTA=${REFERENCE_DIR}/hg38/fasta/Homo_sapiens_assembly38.fasta
    REF_FLAT=/storage/home/kangsm/myDB/annotation/hg38/hg38.refFlat.txt
    CNVKIT_ACCESS_BED=${CNV_FILES_DIR}/hg38_access.10k_CNVkit.bed

    if [ "${WES_LIB}" = "sureselect" ]; then
        TARGET_BED=${CNV_FILES_DIR}/hg38_sureselect.v6.probe.bed
    elif [ "${WES_LIB}" = "sureselect_utr" ]; then
       TARGET_BED=${CNV_FILES_DIR}/hg38_sureselect.v6.utr.probe.bed
    else
       TARGET_BED=${CNV_FILES_DIR}/hg38_twist.exome.2.0.probe.bed  
    fi

else
    GENOME_FASTA=${REFERENCE_DIR}/hg19/fasta/human_g1k_v37_decoy.fasta
    REF_FLAT=${REFERENCE_DIR}/hg19/extra/human_g1k_v37.refFlat.txt
    CNVKIT_ACCESS_BED=${CNV_FILES_DIR}/GRCh37_access.10k_CNVkit.bed

    if [ "${WES_LIB}" = "sureselect" ]; then
        TARGET_BED=${CNV_FILES_DIR}/hg38_sureselect.v6.probe.bed
    elif [ "${WES_LIB}" = "sureselect_utr" ]; then
       TARGET_BED=${CNV_FILES_DIR}/hg38_sureselect.v6.utr.probe.bed
    else
       TARGET_BED=${CNV_FILES_DIR}/hg38_twist.exome.2.0.probe.bed  
    fi

fi

RUN_SINGULARITY="singularity exec -B /storage,/data /storage/images"
CNVKIT_RUN_CMD="${RUN_SINGULARITY}/cnvkit-0.9.10.sif cnvkit.py"

#----------------------------------------------------------------------------------------------------------------------------------------------------#

#---| ANALYSIS RUN |---------------------------------------------------------------------------------------------------------------------------------#

if [ "${RUN_MODE}" = "paired" ]; then

    #---| PAIRED MODE |--------------------------------------------------------#
    CNVKIT_RES_DIR=${CNV_RES_DIR}/CNVkit_paired
    if [ ! -d ${CNVKIT_RES_DIR} ]; then mkdir -p ${CNVKIT_RES_DIR}; fi
    if [ -z "${CNV_RUN_DIR}" ]; then CNV_RUN_DIR=${CNVKIT_RES_DIR}; fi

    ${CNVKIT_RUN_CMD} batch ${TUMOR_BAM} \
        --normal ${NORMAL_BAM} \
        --targets ${TARGET_BED} \
        --fasta ${GENOME_FASTA} \
        --access ${CNVKIT_ACCESS_BED} \
        --output-reference ${CNVKIT_RES_DIR}/${SAMPLE_ID}/${SAMPLE_ID}_reference.cnn \
        --output-dir ${CNVKIT_RES_DIR}/${SAMPLE_ID} \
        -p 2 \
        --diagram --scatter    
    #--------------------------------------------------------------------------#   
else
    #---| PON MODE |-----------------------------------------------------------#
    CNVKIT_RES_DIR=${CNV_RES_DIR}/CNVkit_pon
    if [ ! -d ${CNVKIT_RES_DIR} ]; then mkdir -p ${CNVKIT_RES_DIR}; fi
    if [ -z "${CNV_RUN_DIR}" ]; then CNV_RUN_DIR=${CNVKIT_RES_DIR}; fi

    ${CNVKIT_RUN_CMD} batch ${TUMOR_BAM} \
        --normal ${PON_BAM_DIR}/*.bam \
        --targets ${TARGET_BED} \
        --fasta ${GENOME_FASTA} \
        --access ${CNVKIT_ACCESS_BED} \
        --output-reference ${CNVKIT_RES_DIR}/${SAMPLE_ID}/${SAMPLE_ID}_reference.cnn \
        --output-dir ${CNVKIT_RES_DIR}/${SAMPLE_ID} \
        -p 2 \
        --diagram --scatter
    #--------------------------------------------------------------------------#  
fi
