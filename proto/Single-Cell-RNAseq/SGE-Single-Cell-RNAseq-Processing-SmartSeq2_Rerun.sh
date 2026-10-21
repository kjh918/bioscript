
#! /bin/bash

#$ -N SC-SmartSeq2
#$ -S /bin/bash
#$ -pe mpi 4
#$ -t 390-1027
#$ -tc 25
#$ -e /data/scRNAseq/02.Public-Datasets/004.GSE151334.PRJNA635600.cell-line-smartseq2/meta/runlog/sge.smartseq2.processing2.err
#$ -o /data/scRNAseq/02.Public-Datasets/004.GSE151334.PRJNA635600.cell-line-smartseq2/meta/runlog/sge.smartseq2.processing2.out

BaseDir="/data/scRNAseq/02.Public-Datasets"
BatchID="004.GSE151334.PRJNA635600.cell-line-smartseq2"
Threads="4"

#### processing config file (sample list file)
    # BASE_DIR BATCH_ID SEQ_PLATFORM SEQ_ID  GENOME_ASSEMBLY Threads
    ProcessingListFile=${BaseDir}/${BatchID}/meta/processing.SmartSeq2.list
####
#### sample info
    RunInfo=$(tail -n +2 ${ProcessingListFile} | sed -n -e "${SGE_TASK_ID} p" - )
    SeqPlatform=$(echo ${RunInfo} | awk '{print $1}')
    SeqID=$(echo ${RunInfo} | awk '{print $2}')
    GenomeAssembly=$(echo ${RunInfo} | awk '{print $3}')
    Paired=$(echo ${RunInfo} | awk '{print $4}')
####
#### manual params
    # BaseDir="/data/scRNAseq/02.Public-Datasets"
    # BatchID="002.Tabula-Muris-Senis"
    # SeqPlatform="BulkRNAseq" # BulkRNAseq scRNAseq SmartSeq2   
    # SeqPlatform="scRNAseq"   # BulkRNAseq scRNAseq SmartSeq2   
    # SeqPlatform="SmartSeq2"  # BulkRNAseq scRNAseq SmartSeq2   
    # SeqID=
    # GenomeAssembly=
    # Threads=
####
#### Run Processing

#### Swtools Run Commands
    RunFeatureCounts="/storage/apps/subread-2.1.1/bin/featureCounts"
    RunSalmon="/storage/apps/salmon-1.10.0/bin/salmon"
####
#### Resources
    ResourceDir="/storage/references_and_index/RNAseq"
    GENE_GTF="${ResourceDir}/${GenomeAssembly}/genes/genes.gtf"
    TRANSCRIPT_FASTA="${ResourceDir}/${GenomeAssembly}/fasta/transcripts.fa"
    GENE_MODEL_BED="${ResourceDir}/${GenomeAssembly}/genes/gene.model.strandness.infer.bed"

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

####

    FeatureCountStrandParam="0"
    if [[ "${Paired}" = "True" ]]; then
        SalmonLibraryType="IU"
    else
        SalmonLibraryType="U"
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
