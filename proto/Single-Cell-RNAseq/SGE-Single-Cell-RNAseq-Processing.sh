
#! /bin/bash

#$ -N SC-RNAseq
#$ -S /bin/bash
#$ -pe mpi THREADS_HERE
#$ -t 1-TOTAL_TASKS_HERE
#$ -tc PARALLEL_JOBS_HERE
#$ -e ERROR_LOG_HERE
#$ -o OUT_LOG_HERE

BaseDir="BASE_DIR_HERE"
BatchID="BATCH_ID_HERE"
Threads="THREADS_HERE"

#### processing config file (sample list file)
    # BASE_DIR BATCH_ID SEQ_PLATFORM SEQ_ID  GENOME_ASSEMBLY Threads
    ProcessingListFile=${BaseDir}/${BatchID}/meta/processing.samples.list
####
#### sample info
    RunInfo=$(tail -n +2 ${ProcessingListFile} | sed -n -e "${SGE_TASK_ID} p" - )
    # BaseDir=$(echo ${RunInfo} | awk '{print $1}')
    # BatchID=$(echo ${RunInfo} | awk '{print $2}')
    SeqPlatform=$(echo ${RunInfo} | awk '{print $3}')
    SeqID=$(echo ${RunInfo} | awk '{print $4}')
    GenomeAssembly=$(echo ${RunInfo} | awk '{print $5}')
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

    if [[ -z "${Threads}" ]]; then Threads="8"; fi

    /storage/home/kangsm/runScripts/00_Pipeline_Development/Single-Cell-RNAseq/Single-Cell-RNAseq-NGS-Processing.sh \
        --BASE_DIR ${BaseDir} \
        --BATCH_ID ${BatchID} \
        --SEQ_PLATFORM ${SeqPlatform} \
        --SEQ_ID ${SeqID} \
        --THREADS ${Threads} \
        --GENOME_ASSEMBLY ${GenomeAssembly}



