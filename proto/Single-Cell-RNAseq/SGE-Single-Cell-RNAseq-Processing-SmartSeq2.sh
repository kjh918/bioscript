
#! /bin/bash

#$ -N SmartSeq2
#$ -S /bin/bash
#$ -pe mpi 8
#$ -t 2-1027
#$ -tc 8
#$ -e /data/scRNAseq/02.Public-Datasets/004.GSE151334.PRJNA635600.cell-line-smartseq2/meta/runlog/sge.smartseq2.processing.err
#$ -o /data/scRNAseq/02.Public-Datasets/004.GSE151334.PRJNA635600.cell-line-smartseq2/meta/runlog/sge.smartseq2.processing.out

BaseDir="/data/scRNAseq/02.Public-Datasets"
BatchID="004.GSE151334.PRJNA635600.cell-line-smartseq2"
Threads="8"

#### processing config file (sample list file)
    # BASE_DIR BATCH_ID SEQ_PLATFORM SEQ_ID  GENOME_ASSEMBLY Threads
    ProcessingListFile=${BaseDir}/${BatchID}/meta/processing.SmartSeq2.list
####
#### sample info
    RunInfo=$(tail -n +2 ${ProcessingListFile} | sed -n -e "${SGE_TASK_ID} p" - )
    # BaseDir=$(echo ${RunInfo} | awk '{print $1}')
    # BatchID=$(echo ${RunInfo} | awk '{print $2}')
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
    # Paired=
####
#### Run Processing

    #if [[ -z "${Threads}" ]]; then Threads="8"; fi

    /storage/home/kangsm/runScripts/00_Pipeline_Development/Single-Cell-RNAseq/Single-Cell-RNAseq-NGS-Processing.sh \
        --BASE_DIR ${BaseDir} \
        --BATCH_ID ${BatchID} \
        --SEQ_PLATFORM ${SeqPlatform} \
        --SEQ_ID ${SeqID} \
        --THREADS ${Threads} \
        --GENOME_ASSEMBLY ${GenomeAssembly} \
        --PAIRED ${Paired}



