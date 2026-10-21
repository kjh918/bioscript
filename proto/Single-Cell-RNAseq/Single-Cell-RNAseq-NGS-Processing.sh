
#### Script Usage
    usage() { echo "  Usage: Single-Cell-RNAseq-NGS-Processing.sh [ARGUMENTS]...

        --BASE_DIR          [ REQUIRED ] Base folder for processing and analysis. 

        --BATCH_ID          [ REQUIRED ] NGS data batch-id. 
                            Batch-id also used as folder name of data batch.

        --SEQ_PLATFORM      [ REQUIRED ] RNAseq Type. 
                            available values : BulkRANseq, SmartSeq2, scRNAseq
                            
        --SEQ_ID            [ REQUIRED ] Sample's seq-id. 
                            if RNAseq Type = scRNAseq, SeqID will be used as SeqBatchID
                            
        --THREADS           The number of cores to use
                            default = 8

        --GENOME_ASSEMBLY   [ REQUIRED ] Reference genome assembly version.

        --PAIRED            [ REQUIRED ] PE sequencing type (True) or not (False)

        -h | --help         Print this usage 
        "
    }
####
#### Option Args Parsing
    ARGS=$(getopt -a -o h: --long BASE_DIR:,BATCH_ID:,SEQ_PLATFORM:,SEQ_ID:,THREADS:,GENOME_ASSEMBLY:,PAIRED:,help -- "$@" )
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
        --SEQ_PLATFORM )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  SeqPlatform=$2 ; shift 2 ;;
            esac  ;;
        --SEQ_ID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  SeqID=$2 ; shift 2 ;;
            esac  ;;
        --THREADS )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 8" ; exit 2 ;;
                *)  Threads=$2 ; shift 2 ;;
            esac ;;
        --GENOME_ASSEMBLY )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  GenomeAssembly=$2 ; shift 2 ;;
            esac ;; 
        --PAIRED )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  Paired=$2 ; shift 2 ;;
            esac ;; 
        -h | --help )
        usage >&2 ; exit 2 ;;   
        --) shift ; break ;;
        *)  usage >&2 ; exit 2 ;;
    esac
    done
####

#### Manual Params
    # BaseDir="/data/scRNAseq/02.Public-Datasets"
    # BatchID="004.GSE151334.PRJNA635600.cell-line-smartseq2"
    # SeqPlatform="SmartSeq2" # BulkRNAseq scRNAseq SmartSeq2   
    # SeqID="CCL0539"
    # Threads=16
    # GenomeAssembly="hg38_ercc.plus"
    # Paired="True"
####
#### Param Check
    if [[ -z "${BaseDir}" ||  -z "${BatchID}" || -z "${SeqPlatform}" || -z "${SeqID}" ]]; then
        echo "Some or All Required Params NOT Found."
        exit 2
    fi
    if [[ -z "${Threads}" ]]; then Threads=8; fi
####
#### Swtools Run Commands
    RunFastqc="singularity exec -B/storage,/data /storage/images/fastqc-0.12.1.sif fastqc"
    RunFastp="singularity exec -B/storage,/data /storage/images/fastp-0.24.0.sif fastp"
    RunStar="singularity exec -B /storage,/data /storage/images/star-2.7.11.sif  STAR"
    RunInferStrand="python3 /storage/apps/RSeQC-5.0.1/scripts/infer_experiment.py"
    InferStrandConvertScript="/storage/apps/RSeQC-5.0.1/scripts/infer_result_to_config.R"
    RunFeatureCounts="/storage/apps/subread-2.1.1/bin/featureCounts"
    RunSalmon="/storage/apps/salmon-1.10.0/bin/salmon"
####
#### Resources
    ResourceDir="/storage/references_and_index/RNAseq"
    STAR_INDEX="${ResourceDir}/${GenomeAssembly}/star"
    GENE_GTF="${ResourceDir}/${GenomeAssembly}/genes/genes.gtf"
    TRANSCRIPT_FASTA="${ResourceDir}/${GenomeAssembly}/fasta/transcripts.fa"
    GENE_MODEL_BED="${ResourceDir}/${GenomeAssembly}/genes/gene.model.strandness.infer.bed"

    if [[ "${GenomeAssembly}" = "mm39" ]]; then 
        CELLRANGER_REFERENCE_DIR="/storage/references_and_index/cellranger/refdata-gex-GRCm39-2024-A"
    elif [[ "${GenomeAssembly}" = "hg38" ]]; then 
        CELLRANGER_REFERENCE_DIR="/storage/references_and_index/cellranger/refdata-gex-GRCh38-2024-A"
    elif [[ "${GenomeAssembly}" = "hg38_mm39" ]]; then 
        CELLRANGER_REFERENCE_DIR="/storage/references_and_index/cellranger/refdata-gex-GRCh38_and_GRCm39-2024-A"
    fi

####


#### Data Folders
    RawFastqDir="${BaseDir}/${BatchID}/fastq/${SeqPlatform}"
    ResDir="${BaseDir}/${BatchID}/${SeqPlatform}/${SeqID}"
    mkdir -p "${ResDir}"
    if [[ "${SeqPlatform}" = "BulkRNAseq" || "${SeqPlatform}" = "SmartSeq2" ]]; then
        TrimFastqDir="${ResDir}/fastq"
        BamDir="${ResDir}/bam"
        QuantDir="${ResDir}/quant"
        QcDir="${ResDir}/qcfiles"
        mkdir -p "${TrimFastqDir}"
        mkdir -p "${BamDir}"
        mkdir -p "${QuantDir}"
        mkdir -p "${QcDir}"
    fi
####

#### FASTQ-QC

    if [[ "${SeqPlatform}" = "BulkRNAseq" || "${SeqPlatform}" = "SmartSeq2" ]]; then
        # fastqc    
        if [[ ! -f "${QcDir}/fastqc.done" ]]; then
            export LC_ALL=C
            export LANG=C
            # run fastqc
            FASTQC_INPUT=()
            if [[ "${Paired}" = "True" ]]; then
                FASTQC_INPUT+=("${RawFastqDir}/${SeqID}_R1.fastq.gz")
                FASTQC_INPUT+=("${RawFastqDir}/${SeqID}_R2.fastq.gz")
            else
                FASTQC_INPUT+=("${RawFastqDir}/${SeqID}.fastq.gz")
            fi
            ${RunFastqc} --extract --threads 2 --outdir ${QcDir} ${FASTQC_INPUT[@]}
            # prevent re-run
            touch "${QcDir}/fastqc.done"
        fi

        # fastp
        if [[ ! -f "${QcDir}/fastp.done" ]]; then 
            # fastp args
            FASTP_ARGS=()
            if [[ "${Paired}" = "True" ]]; then
                FASTP_ARGS+=("--in1 ${RawFastqDir}/${SeqID}_R1.fastq.gz")
                FASTP_ARGS+=("--in2 ${RawFastqDir}/${SeqID}_R2.fastq.gz")
                FASTP_ARGS+=("--out1 ${TrimFastqDir}/${SeqID}.trimmed_R1.fastq.gz")
                FASTP_ARGS+=("--out2 ${TrimFastqDir}/${SeqID}.trimmed_R2.fastq.gz")
            else
                FASTP_ARGS+=("--in1 ${RawFastqDir}/${SeqID}.fastq.gz")
                FASTP_ARGS+=("--out1 ${TrimFastqDir}/${SeqID}.trimmed.fastq.gz")
            fi
            FASTP_ARGS+=("--trim_poly_g")
            FASTP_ARGS+=("--length_required 20")
            FASTP_ARGS+=("--thread 6")           
            FASTP_ARGS+=("--json ${QcDir}/${SeqID}.fastp.json")
            FASTP_ARGS+=("--html ${QcDir}/${SeqID}.fastp.html")
            # run fastp
            ${RunFastp} ${FASTP_ARGS[@]}
            # prevent re-run
            touch "${QcDir}/fastp.done"
        fi
    fi
####

#### BulkRNAseq and SmartSeq2 

    # STAR Alignment : BulkRNAseq 
    if [[ "${SeqPlatform}" = "BulkRNAseq" ]]; then
        if [[ ! -f "${QcDir}/STAR.BulkRNAseq.done" ]]; then
            # STAR args
            STAR_ARGS=()
            if [[ "${Paired}" = "True" ]]; then
                STAR_ARGS+=("--readFilesIn ${TrimFastqDir}/${SeqID}.trimmed_R1.fastq.gz ${TrimFastqDir}/${SeqID}.trimmed_R2.fastq.gz")
            else
                STAR_ARGS+=("--readFilesIn ${TrimFastqDir}/${SeqID}.trimmed.fastq.gz")
            fi
            STAR_ARGS+=("--genomeDir ${STAR_INDEX}")
            STAR_ARGS+=("--sjdbGTFfile ${GENE_GTF}")
            STAR_ARGS+=("--runThreadN ${Threads}")
            STAR_ARGS+=("--readFilesCommand zcat")
            STAR_ARGS+=("--quantMode TranscriptomeSAM")
            STAR_ARGS+=("--genomeLoad NoSharedMemory")
            STAR_ARGS+=("--twopassMode Basic")
            STAR_ARGS+=("--outSAMtype BAM SortedByCoordinate")
            STAR_ARGS+=("--outFileNamePrefix ${BamDir}/${SeqID}.")
            # run STAR
            ${RunStar} ${STAR_ARGS[@]}
            # prevent re-run
            touch "${QcDir}/STAR.BulkRNAseq.done"
        fi
        # infer strandness from bam
        export PYTHONPATH=$PYTHONPATH:/storage/apps/RSeQC-5.0.1/lib
        export PATH=$PATH:/storage/apps/RSeQC-5.0.1/scripts
        ${RunInferStrand} -i "${BamDir}/${SeqID}.Aligned.sortedByCoord.out.bam" -r "${GENE_MODEL_BED}" > ${QcDir}/rseqc_strand_infer_result.txt
        Rscript ${InferStrandConvertScript} ${QcDir}/rseqc_strand_infer_result.txt
        StrandInfo=$(cat ${QcDir}/stranded_config.txt)      
        FeatureCountStrandParam=$(echo ${StrandInfo} | awk '{print $1}')
        SalmonLibraryType=$(echo ${StrandInfo} | awk '{print $2}')
    fi

    # STAR Alignment : SmartSeq2 (FACS)
    if [[ "${SeqPlatform}" = "SmartSeq2" ]]; then
        if [[ ! -f "${QcDir}/STAR.SmartSeq2.done" ]]; then
            # STAR args
            STAR_ARGS=()
            if [[ "${Paired}" = "True" ]]; then
                STAR_ARGS+=("--readFilesIn ${TrimFastqDir}/${SeqID}.trimmed_R1.fastq.gz ${TrimFastqDir}/${SeqID}.trimmed_R2.fastq.gz")
            else
                STAR_ARGS+=("--readFilesIn ${TrimFastqDir}/${SeqID}.trimmed.fastq.gz")
            fi
            STAR_ARGS+=("--genomeDir ${STAR_INDEX}")
            STAR_ARGS+=("--sjdbGTFfile ${GENE_GTF}")
            STAR_ARGS+=("--runThreadN ${Threads}")
            STAR_ARGS+=("--readFilesCommand zcat")
            STAR_ARGS+=("--quantMode TranscriptomeSAM")
            STAR_ARGS+=("--genomeLoad NoSharedMemory")
            STAR_ARGS+=("--outSAMtype BAM SortedByCoordinate")
            STAR_ARGS+=("--outFileNamePrefix ${BamDir}/${SeqID}.")
            # run STAR
            ${RunStar} ${STAR_ARGS[@]}
            # prevent re-run
            touch "${QcDir}/STAR.SmartSeq2.done" 
        fi
        # assign strandness
        FeatureCountStrandParam="0"
        if [[ "${Paired}" = "True" ]]; then
            SalmonLibraryType="IU"
        else
            SalmonLibraryType="U"
        fi
    fi

    # Quantification 
    if [[ "${SeqPlatform}" = "BulkRNAseq" || "${SeqPlatform}" = "SmartSeq2" ]]; then
    
        # featureCounts
        OutputDirFeatureCounts=${QuantDir}/featureCounts
        mkdir -p ${OutputDirFeatureCounts}
        FEATURECOUNTS_ARGS=()
        FEATURECOUNTS_ARGS+=("-T ${Threads}")
        FEATURECOUNTS_ARGS+=("-F GTF")
        FEATURECOUNTS_ARGS+=("-s ${FeatureCountStrandParam}")
        FEATURECOUNTS_ARGS+=("-a ${GENE_GTF}")
        if [[ "${Paired}" = "True" ]]; then
            FEATURECOUNTS_ARGS+=("-p")
        fi
        ${RunFeatureCounts} ${FEATURECOUNTS_ARGS[@]} \
            -o ${OutputDirFeatureCounts}/featureCounts.tsv  \
            ${BamDir}/${SeqID}.Aligned.sortedByCoord.out.bam
        rm ${OutputDirFeatureCounts}/*.tmp

        # salmon (align-based)
        OutputDirSalmon=${QuantDir}/salmon
        mkdir -p ${OutputDirSalmon}
        ${RunSalmon} quant -p ${Threads} -l ${SalmonLibraryType} \
            -t ${TRANSCRIPT_FASTA} \
            -a ${BamDir}/${SeqID}.Aligned.toTranscriptome.out.bam \
            -o ${OutputDirSalmon} \
            --quiet
    fi
####

#### scRNAseq (cellranger)

    if [[ "${SeqPlatform}" = "scRNAseq" ]]; then
        
        # input fastq folder
        scRNAseqFastqDir="${BaseDir}/${BatchID}/fastq/${SeqPlatform}/${SeqID}"
        AnalysisRunID="${SeqID}_Force3000"
        ResDir="${BaseDir}/${BatchID}/${SeqPlatform}/${AnalysisRunID}"
        mkdir -p ${ResDir}

        # run CellRanger
        cellranger count \
            --id ${AnalysisRunID} \
            --transcriptome ${CELLRANGER_REFERENCE_DIR} \
            --fastqs ${scRNAseqFastqDir} \
            --sample ${SeqID} \
            --localcores ${Threads} \
            --localmem 64 \
            --create-bam true \
            --output-dir ${ResDir}
            
    fi

####


