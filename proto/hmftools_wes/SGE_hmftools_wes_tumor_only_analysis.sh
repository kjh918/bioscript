

#! /bin/bash

#$ -N WES_T_VAR_CALL
#$ -S /bin/bash
#$ -pe mpi 12
#$ -t 1-18
#$ -tc 6
#$ -e /storage/home/kangsm/analysis/gcx.pipeline.use.hmftools.in.wes/01.preprocess/runlog/tumor.only.analysis.$TASK_ID.err
#$ -o /storage/home/kangsm/analysis/gcx.pipeline.use.hmftools.in.wes/01.preprocess/runlog/tumor.only.analysis.$TASK_ID.out

#---/ Sample List File
    SampleListFile="/storage/home/kangsm/analysis/gcx.pipeline.use.hmftools.in.wes/01.preprocess/analysis.list"
#-------------------------------------------------------------------------------

#---/ Sample-Info from Sample-List-File |---------------------------------------
    # SeqID RGID BatchID WesLibrary
    RunInfo=$(sed -n -e "${SGE_TASK_ID} p" ${SampleListFile})
  
    BatchID=$(echo ${RunInfo} | awk '{print $3}')
    TumorSeqID=$(echo ${RunInfo} | awk '{print $1}')
    NormalSeqID=$(echo ${RunInfo} | awk '{print $2}')
    WesLibrary=$(echo ${RunInfo} | awk '{print $4}')

#-------------------------------------------------------------------------------

#---/ Manual Params |-----------------------------------------------------------
    BaseDir="/data/wes"
    #BatchID="WES-panel.of.normals-hg19-sureselect.v6"
    #SampleListFile="/data/pon/WES-panel.of.normals-hg19-sureselect.v6/meta/sample.list"
    Threads=12
    GenomeAssembly="hg19"
    #WesLibrary="sureselect.v6"
#-------------------------------------------------------------------------------

    /storage/home/kangsm/runScripts/WES_v2/WES_M02-1_VARIANT-CALLING-Tumor-Only.sh \
        --BASE_DIR ${BaseDir} \
        --BATCH_ID ${BatchID} \
        --TUMOR_SEQ_ID ${TumorSeqID} \
        --WES_LIB ${WesLibrary} \
        --GENOME_ASSEMBLY ${GenomeAssembly} \
        --THREADS ${Threads}

#-------------------------------------------------------------------------------