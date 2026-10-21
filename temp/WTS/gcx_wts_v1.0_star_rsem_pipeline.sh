#!/bin/bash

#---| Input Options |
    usage() {
        echo "Usage: GCX-WES-TumorOnlyVariantCalling.sh [ options ]
        [ -d | --seqFolder ]   : data processing base directory
        [ -s | --seqID ]       : sample ID (Used as sample's processing product folder name)
        [ -l | --logFile ]     : logFile name
        [ -t | --threads ]     : cpu cores 
        [ --readLength ]       : library read length 100 or 150, default = 150 
        [ --stranded ]         : stranded seq or not 'yes' or 'no', default= 'no'
        [ -h | --help ]        : Print this usage 
        "
    }
  #---| Set Defualt Values |
  THREADS=10
  READ_LENGTH=150
  STRANDED="no"
  #---/
#---/

#---| Parse Option Arguments |
    ARGS=$(getopt -a -o d:s:l:t:h: --long seqFolder:,seqID:,logFile:,threads:,help -- "$@" )
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
        -l | --logFile )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  LOG_NAME=$2 ; shift 2 ;;
            esac ;;
        -t | --threads )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 10" ; exit 2 ;;
                *)  THREADS=$2 ; shift 2 ;;
            esac ;;
        --readLength )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 150" ; exit 2 ;;
                *)  READ_LENGTH=$2 ; shift 2 ;;
            esac ;;        
        -h | --help )
        usage >&2 ; exit 2 ;;   
        --) shift ; break ;;
        *)  usage >&2 ; exit 2 ;;
    esac
    done
#---/

#---| FOLDERS |
    DATA_DIR=/data/wts/${SEQ_FOLDER}
    FASTQ_DIR=${DATA_DIR}/${SEQ_ID}/fastq
    BAM_DIR=${DATA_DIR}/${SEQ_ID}/bam
    QUANT_DIR=${DATA_DIR}/${SEQ_ID}/quant
    QC_DIR=${DATA_DIR}/${SEQ_ID}/qcfiles

    if [ ! -d ${FASTQ_DIR} ]; then mkdir -p ${FASTQ_DIR}; else rm -rf ${FASTQ_DIR}/*trimmed*; fi
    if [ ! -d ${BAM_DIR} ]; then mkdir -p ${BAM_DIR}; fi
    if [ ! -d ${QUANT_DIR} ]; then mkdir -p ${QUANT_DIR}; fi
    if [ ! -d ${QC_DIR} ]; then mkdir -p ${QC_DIR}; fi
#---/

#---| REFERENCES |
    if [ "${READ_LENGTH}" = "150"  ]; then 
        STAR_INDEX=/storage/references_and_index/hg19/rsem-index-star-pe150
    else
        STAR_INDEX=/storage/references_and_index/hg19/rsem-index-star-pe100
    fi
    GTF_FILE=/storage/references_and_index/hg19/Homo_sapiens_assembly19.gtf
    rRNA_INTERVAL=/storage/references_and_index/hg19/extra/human_g1k_v37.rRNA.interval_list
    REF_FLAT=/storage/references_and_index/hg19/extra/human_g1k_v37.refFlat.txt
    GENOME_FASTA=/storage/references_and_index/hg19/fasta/human_g1k_v37_decoy.fasta
#---/

#---| RUN COMMANDS
    RUN_SINGULARITY_IMG="singularity exec -B /storage,/data /storage/images"
    RUN_FASTP="singularity exec -B /storage,/data /storage/images/fastp-0.23.4.sif"
    RUN_STAR="singularity exec -B /storage,/data /storage/images/star-2.7.11.sif"
    RUN_RSEM="singularity exec -B /storage,/data /storage/images/rsem-1.3.3.sif"
#---/

#---| LOGFILE |
    LOGFILE=${DATA_DIR}/${SEQ_ID}/${LOG_NAME}
    touch ${LOGFILE}
#---/

#---\ 
    echo " " >> ${LOGFILE}
    echo "GCX RNAseq processing log" >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " SEQ FOLDER           : ${DATA_DIR}" >> ${LOGFILE}
    echo " SEQ ID               : ${SEQ_ID}" >> ${LOGFILE}
    echo " REFERENCE            : hg19 (GRCh37)" >> ${LOGFILE}
    echo " ALIGNER              : STAR-2.7.11a (TwoPass Mode)" >> ${LOGFILE}
    echo " QUANTIFICATION       : RSEM-1.3.3" >> ${LOGFILE}
    if [ "${STRANDED}" = "yes" ]; then
        echo " STRANDED LIBRARY     : YES" >>  ${LOGFILE}
    else
        echo " STRANDED LIBRARY     : NO" >>  ${LOGFILE}
    fi
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo "$(date) | Start Processing." >> ${LOGFILE}
#---/

#---| FASTQC |
#---\
    echo "$(date) | FastQC Start." >> ${LOGFILE}
#---/
    /storage/apps/bin/fastqc -t ${THREADS} -o ${QC_DIR} ${FASTQ_DIR}/*.fastq.gz
#---\
    echo "$(date) | FastQC Finished." >> ${LOGFILE}
#---/

#---| FASTP FILTERING |
#---\
    echo "$(date) | Fastq Trimming and Filtering Start." >> ${LOGFILE}
    echo "$(date) | trimming - minimun length cut-off = 105" >> ${LOGFILE}
    echo "$(date) | trimming - base q-score   cut-off = 20" >> ${LOGFILE}
#---/
    ${RUN_SINGULARITY_IMG}/fastp-0.23.4.sif fastp \
        --in1 ${FASTQ_DIR}/${SEQ_ID}_R1.fastq.gz \
        --in2 ${FASTQ_DIR}/${SEQ_ID}_R2.fastq.gz \
        --json ${DATA_DIR}/${SEQ_ID}/qcfiles/${SEQ_ID}.fastp.json \
        --html ${DATA_DIR}/${SEQ_ID}/qcfiles/${SEQ_ID}.fastp.html \
        --length_required 105 \
        --qualified_quality_phred 20 \
        --thread ${THREADS} \
        --trim_poly_g \
        --detect_adapter_for_pe \
        --out1 ${FASTQ_DIR}/${SEQ_ID}.trimmed_R1.fastq.gz \
        --out2 ${FASTQ_DIR}/${SEQ_ID}.trimmed_R2.fastq.gz
#---\
    echo "$(date) | Fastq Trimming and Filtering Finished." >> ${LOGFILE}
#---/

#---| STAR ALIGNMENT |
#---\
    echo "$(date) | STAR Alignment Start." >> ${LOGFILE}
#---/
    ${RUN_SINGULARITY_IMG}/star-2.7.11.sif STAR \
        --readFilesCommand zcat \
        --quantMode TranscriptomeSAM \
        --twopassMode Basic \
        --genomeLoad NoSharedMemory \
        --runThreadN ${THREADS} \
        --outFilterMatchNminOverLread 0.33 \
        --outFilterScoreMinOverLread 0.33 \
        --readQualityScoreBase 20 \
        --outFilterMismatchNmax 15 \
        --outFilterMismatchNoverLmax 0.1 \
        --outSAMtype BAM Unsorted \
        --outSAMunmapped None\
        --genomeDir ${STAR_INDEX} \
        --readFilesIn ${FASTQ_DIR}/${SEQ_ID}.trimmed_R1.fastq.gz ${FASTQ_DIR}/${SEQ_ID}.trimmed_R2.fastq.gz \
        --outFileNamePrefix ${BAM_DIR}/${SEQ_ID} \
        --sjdbGTFfile ${GTF_FILE}
#---\
    echo "$(date) | Sorting BAM for QC." >> ${LOGFILE}
#---/  
    ${RUN_SINGULARITY_IMG}/samtools-1.18.sif samtools sort -@ ${THREADS} \
        -o ${BAM_DIR}/${SEQ_ID}.Aligned.out.Sorted.bam \
        ${BAM_DIR}/${SEQ_ID}.Aligned.out.bam     
#---\
    echo "$(date) | Sorted BAM Generated." >> ${LOGFILE}
    echo "$(date) | STAR Alignment Finished." >> ${LOGFILE}
#---/

#---| QUANTIFICATION |
#---\
    echo "$(date) | RSEM Quantification Start." >> ${LOGFILE}
#---/
    if [ "${STRANDED}" = "yes" ]; then
        ${RUN_SINGULARITY_IMG}/rsem-1.3.3.sif rsem-calculate-expression \
            --paired-end --no-bam-output --estimate-rspd \
            --strandedness reverse \
            --num-threads ${THREADS} \
            --alignments ${BAM_DIR}/${SEQ_ID}.Aligned.toTranscriptome.out.bam \
            ${STAR_INDEX}/human_g1k_v37 \
            ${QUANT_DIR}/${SEQ_ID}
    else
        ${RUN_SINGULARITY_IMG}/rsem-1.3.3.sif rsem-calculate-expression \
            --paired-end --no-bam-output --estimate-rspd \
            --num-threads ${THREADS} \
            --alignments ${BAM_DIR}/${SEQ_ID}.Aligned.toTranscriptome.out.bam \
            ${STAR_INDEX}/human_g1k_v37 \
            ${QUANT_DIR}/${SEQ_ID}
    fi
#---\
    echo "$(date) | RSEM Quantification Finished." >> ${LOGFILE}
#---/

#---| RNAseq ALIGNMENT QC |
#---\
    echo "$(date) | RNAseq Alignment QC Start." >> ${LOGFILE}
    echo "$(date) |     QC.1 [ GATK CollectRnaSeqMetrics ] Start." >> ${LOGFILE}
#---/
    if [ "${STRANDED}" = "yes" ]; then STRAND_OPT=SECOND_READ_TRANSCRIPTION_STRAND; else STRAND_OPT=NONE; fi
    ${RUN_SINGULARITY_IMG}/gatk-4.4.0.0.sif gatk CollectRnaSeqMetrics \
        --java-options "-XX:ConcGCThreads=10 -Xmx32768m" \
        --INPUT ${BAM_DIR}/${SEQ_ID}.Aligned.out.Sorted.bam \
        --OUTPUT ${QC_DIR}/${SEQ_ID}.gatk.rnaseq.metrics \
        --REF_FLAT ${REF_FLAT} \
        --STRAND_SPECIFICITY ${STRAND_OPT} \
        --MINIMUM_LENGTH 500 \
        --CHART_OUTPUT ${QC_DIR}/${SEQ_ID}.gatk.rnaseq.metric.pdf \
        --RIBOSOMAL_INTERVALS ${rRNA_INTERVAL} \
        --REFERENCE_SEQUENCE ${GENOME_FASTA}
#---\
    echo "$(date) |     QC.1 [ GATK CollectRnaSeqMetrics ] Finished." >> ${LOGFILE}
    echo "$(date) |     QC.2 [ RNA-SeQC                  ] Start." >> ${LOGFILE}
#---/
    ${RUN_SINGULARITY_IMG}/rna-seqc-2.4.2.sif rnaseqc \
        --sample=$SEQ_ID \
        --fasta=$GENOME_FASTA \
        --base-mismatch=15 \
        --gene-length=500 \
        --coverage \
        --detection-threshold=3 \
        ${GTF_FILE} \
        ${BAM_DIR}/${SEQ_ID}.Aligned.out.Sorted.bam \
        ${QC_DIR}
#---\
    echo "$(date) |     QC.2 [ RNA-SeQC                  ] Finished." >> ${LOGFILE}
    echo "$(date) | RNAseq Alignment QC Finished." >> ${LOGFILE}
#---/
    mkdir -p ${QC_DIR}/qc_res
    SCRIPT_DIR=/storage/home/kangsm/runScripts
    cp ${QUANT_DIR}/${SEQ_ID}.stat/* ${QC_DIR}/
    ##
    ${RUN_SINGULARITY_IMG}/multiqc-1.16.sif multiqc --force \
        --filename ${SEQ_ID}_qc \
        --outdir ${QC_DIR}/qc_res \
        --data-dir ${QC_DIR}





#---\
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo "$(date) | Processing Finished." >> ${LOGFILE}
#---/
