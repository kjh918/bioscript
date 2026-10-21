

#! /bin/bash

#$ -N hmftools_wes
#$ -S /bin/bash
#$ -pe mpi 10
#$ -t 1-27
#$ -tc 10
#$ -e /storage/home/kangsm/analysis/gcx.pipeline.use.hmftools.in.wes/01.preprocess/runlog/preprocessing.$TASK_ID.err
#$ -o /storage/home/kangsm/analysis/gcx.pipeline.use.hmftools.in.wes/01.preprocess/runlog/preprocessing.$TASK_ID.out

#---/ Sample List File
    SampleListFile="/storage/home/kangsm/analysis/gcx.pipeline.use.hmftools.in.wes/01.preprocess/sample.list"
#-------------------------------------------------------------------------------

#---/ Sample-Info from Sample-List-File |---------------------------------------
    # SeqID RGID BatchID WesLibrary
    RunInfo=$(sed -n -e "${SGE_TASK_ID} p" ${SampleListFile})
    
    BatchID=$(echo ${RunInfo} | awk '{print $3}')
    SeqID=$(echo ${RunInfo} | awk '{print $1}')
    RGID=$(echo ${RunInfo} | awk '{print $2}')
    RGSM=$(echo ${RunInfo} | awk '{print $1}')
    WesLibrary=$(echo ${RunInfo} | awk '{print $4}')

#-------------------------------------------------------------------------------

#---/ Manual Params |-----------------------------------------------------------
    BaseDir="/data/wes"
    #BatchID="WES-panel.of.normals-hg19-sureselect.v6"
    #SampleListFile="/data/pon/WES-panel.of.normals-hg19-sureselect.v6/meta/sample.list"
    Threads=16
    GenomeAssembly="hg19"
    #WesLibrary="sureselect.v6"
#-------------------------------------------------------------------------------

#---| HMFTOOLS |----------------------------------------------------------------

    RUN_SIF="singularity exec -B /storage,/data /storage/images"
    RUN_BWA_MEM2="${RUN_SIF}/bwa-mem2-v2.2.1.sif bwa-mem2"
    HMF_RUN_REDUX="${RUN_SIF}/hmftools-redux-v1.2.3.sif redux"

    HMF_GENOME_FASTA="/storage/references_and_index/hg19/hmftools/Homo_sapiens.GRCh37.GATK.illumina.fasta"
    HMF_BWAMEM2_INDEX="/storage/references_and_index/hg19/hmftools/bwa-mem2_index-2.2.1/Homo_sapiens.GRCh37.GATK.illumina.fasta"
    GenomeVersion="37"

    HMF_RESOURCE="/storage/references_and_index/${GenomeAssembly}/hmftools"
    HMF_UNMAP_REGION=${HMF_RESOURCE}/resources_links/UNMAP-REGIONS.tsv
    HMF_MSI_JITTER_SITES=${HMF_RESOURCE}/resources_links/MSI-JITTER-SITES.tsv.gz
    HMF_HIGH_CONFIDENCE_BED=${HMF_RESOURCE}/resources_links/HIGH-CONFIDENCE.bed.gz
    HMF_ENSEMBL_DIR=${HMF_RESOURCE}/resources_links/ENSEMBL-DATA/
    
#-------------------------------------------------------------------------------

    # SortThreads=4
    # BwaThreads=12

    FastqDir=${BaseDir}/${BatchID}/${SeqID}/fastq
    BamDir=${BaseDir}/${BatchID}/${SeqID}/bam
    HmftoolsDir=${BaseDir}/${BatchID}/${SeqID}/hmftools
    mkdir -p ${HmftoolsDir}
    
    # bwa mem2    
    # ${RUN_BWA_MEM2} mem -M -Y -q -K 100000000 -t ${BwaThreads} \
    #     -R "@RG\tID:${RGID}\tSM:${RGSM}\tLB:TWISTEXOME2\tPL:ILLUMINA\tCN:GCX" \
    #     ${HMF_BWAMEM2_INDEX} \
    #     ${FastqDir}/${SeqID}.trimmed_R1.fastq.gz \
    #     ${FastqDir}/${SeqID}.trimmed_R2.fastq.gz \
    #     | samtools sort -@ ${SortThreads} -o ${BamDir}/${SeqID}.hmftools.bwa.mem2.bam - 

    # samtools index -@ ${Threads} ${BamDir}/${SeqID}.hmftools.bwa.mem2.bam 

    # redux
    ReduxOutputDir=${HmftoolsDir}/redux
    mkdir -p ${ReduxOutputDir}

    ${HMF_RUN_REDUX} \
        -sample ${SeqID} \
        -input_bam ${BamDir}/${SeqID}.hmftools.bwa.mem2.bam \
        -output_bam ${ReduxOutputDir}/${SeqID}.redux.bam \
        -output_dir ${ReduxOutputDir} \
        -ref_genome_version ${GenomeVersion} -ref_genome ${HMF_GENOME_FASTA} \
        -unmap_regions ${HMF_UNMAP_REGION} -ref_genome_msi_file ${HMF_MSI_JITTER_SITES} \
        -form_consensus -write_stats -bamtool /opt/conda/bin/samtools \
        -threads ${Threads}

#-------------------------------------------------------------------------------