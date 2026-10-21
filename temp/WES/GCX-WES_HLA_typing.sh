
#!/bin/bash

#---| Input Options |
    usage() {
        echo "Usage: GCX-WES-HLA_Typing.sh [ options ]
        [ -d | --dataDir ]        : data processing base directory
        [ -s | --sampleID ]       : sample ID (Used as sample's processing product folder name)
        [ -l | --logFile ]        : logFile name (file will be created in '/path/to/dataDir/smapleID/')
        [ -t | --inputType ]      : input NGS type. DNA or RNA. (default = DNA)
        [ -h | --help ]           : Print this usage 
        "
    }
  #---| Set Defualt Values |
    INPUT_TYPE="DNA"
  #---/
#---/
#---| Parse Option Arguments |
    ARGS=$(getopt -a -o d:s:l:t:h: --long dataDir:,sampleID:,logFile:,inputType:,help -- "$@" )
    VALID_ARGS=$?
    if [ "$VALID_ARGS" != "0" ]; then 
        usage >&2 
        exit 2
    fi

    eval set -- "$ARGS"
    while :
    do
    case "$1" in
        -d | --dataDir )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  SEQ_FOLDER=$2 ; shift 2 ;;
            esac ;;
        -s | --sampleID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  SEQ_ID=$2 ; shift 2 ;;
            esac  ;;
        -l | --logFile )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  LOG_FILE=$2 ; shift 2 ;;
            esac ;;
        -t | --inputType )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = DNA" ; shift ;;
                *)  INPUT_TYPE=$2 ; shift 2 ;;
            esac ;;
        -h | --help )
        usage >&2 ; exit 2 ;;   
        --) shift ; break ;;
        *)  usage >&2 ; exit 2 ;;
    esac
    done
#---/

#---| OptiType |
    if [ "${INPUT_TYPE}" == "DNA" ]; then
        OPTITYPE_REF_FASTA=/storage/references_and_index/optitype_references/hla_reference_dna.fasta
    else
        OPTITYPE_REF_FASTA=/storage/references_and_index/optitype_references/hla_reference_rna.fasta
    fi
  #-| FOLDERS
    mkdir -p ${SEQ_FOLDER}/${SEQ_ID}/hla
    mkdir -p ${SEQ_FOLDER}/${SEQ_ID}/hla/optitype  
    OPTITYPE_DIR=${SEQ_FOLDER}/${SEQ_ID}/hla/optitype
    FASTQ_DIR=${SEQ_FOLDER}/${SEQ_ID}/fastq
    READS="R1 R2"
  #-| CONDA ENV ACTIVATION
    source /storage/apps/miniconda3/bin/activate HLA
  #-| FASTQ FILTERING
    parallel -j 2 -k "razers3 -i 95 -m 1 -dr 0 -o ${OPTITYPE_DIR}/${SEQ_ID}.optitype_fished_{}.bam ${OPTITYPE_REF_FASTA} ${FASTQ_DIR}/${SEQ_ID}_{}.fastq.gz" ::: ${READS}
    parallel -j 2 -k "/storage/apps/bin/samtools bam2fq ${OPTITYPE_DIR}/${SEQ_ID}.optitype_fished_{}.bam > ${OPTITYPE_DIR}/${SEQ_ID}.optitype_{}.fastq" ::: ${READS}
  #-| Run OptiType
    python /storage/apps/miniconda3/envs/HLA/bin/OptiTypePipeline.py \
        --input ${OPTITYPE_DIR}/${SEQ_ID}.optitype_R1.fastq ${OPTITYPE_DIR}/${SEQ_ID}.optitype_R2.fastq \
        --dna \
        --outdir ${OPTITYPE_DIR} \
        --prefix ${SEQ_ID}.optitype
  #-/
#---/
#---| HLA-LA |
  #-| FOLDERS 
    mkdir -p ${SEQ_FOLDER}/${SEQ_ID}/hla/hla-la
    HLA_LA_DIR=${SEQ_FOLDER}/${SEQ_ID}/hla/hla-la
    BAM_DIR=${SEQ_FOLDER}/${SEQ_ID}/bam
  #-| Run HLA-LA  
    /storage/apps/HLA-LA-1.0.3/src/HLA-LA.pl \
        --BAM ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.recal.bam \
        --graph PRG_MHC_GRCh38_withIMGT \
        --sampleID ${SEQ_ID} \
        --maxThreads 25 \
        --workingDir ${HLA_LA_DIR}
  #-/
#---/


