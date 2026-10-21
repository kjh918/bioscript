
#! /bin/bash

#$ -N hmftools_wes_t_cohort
#$ -S /bin/bash
#$ -pe mpi 16
#$ -t 1-71
#$ -tc 5
#$ -e /data/wgs/GCX_Resources/WGS.panel.of.normal_1000.genome_10.samples/hg19/meta/sge.pre.merge.$JOB_ID.$TASK_ID.err
#$ -o /data/wgs/GCX_Resources/WGS.panel.of.normal_1000.genome_10.samples/hg19/meta/sge.pre.merge.$JOB_ID.$TASK_ID.err

#---/ Manual Params |-----------------------------------------------------------
    BaseDir="/storage/home/kangsm/analysis/gcx.pipeline.use.hmftools.in.wes/02.data_processing"
    BatchID="hg19"
    SampleListFile="/storage/home/kangsm/analysis/gcx.pipeline.use.hmftools.in.wes/01.preprocess/sample.list"
    Threads=16
    GenomeAssembly="hg19"
#-------------------------------------------------------------------------------

#---/ Sample-Info from Sample-List-File |---------------------------------------
    # SeqID  ReadGroupID SampleName
    
    RunInfo=$(sed -n -e "${SGE_TASK_ID} p" ${SampleListFile})
    
    SeqID=$(echo ${RunInfo} | awk '{print $2}')
    ReadGroupID=$(echo ${RunInfo} | awk '{print $3}')
    SampleName=$(echo ${RunInfo} | awk '{print $2}')

#-------------------------------------------------------------------------------

    RunSingularity="singularity exec -B /storage,/data /storage/images"
    RunBwaMem2="${RunSingularity}/bwa-mem2-v2.2.1.sif bwa-mem2"
    HMF_RUN_REDUX="${RunSingularity}/hmftools-redux-v1.2.3.sif redux"
    HMF_GenomeFasta="/storage/references_and_index/hg19/hmftools/Homo_sapiens.GRCh37.GATK.illumina.fasta"
    HMF_BwaMem2Index="/storage/references_and_index/hg19/hmftools/bwa-mem2_index-2.2.1/Homo_sapiens.GRCh37.GATK.illumina.fasta"
    SortThreads=4
    BwaThreads=12
    FastqDir="/storage/home/kangsm/analysis/gcx.pipeline.use.hmftools.in.wes/02.data_processing/fastq_trimmed"
    BamDir="/storage/home/kangsm/analysis/gcx.pipeline.use.hmftools.in.wes/02.data_processing/${SeqID}/bam"
    mkdir -p ${BamDir}
    HMF_RESOURCE="/storage/references_and_index/hg19/hmftools/resources_links"
    HMF_UNMAP_REGION=${HMF_RESOURCE}/UNMAP-REGIONS.tsv
    HMF_MSI_JITTER_SITES=${HMF_RESOURCE}/MSI-JITTER-SITES.tsv.gz
    
    ${RunBwaMem2} mem -M -Y -q -K 100000000 -t ${BwaThreads} \
        -R "@RG\tID:${ReadGroupID}\tSM:${SampleName}\tLB:TWISTEXOME2\tPL:ILLUMINA\tCN:GCX" \
        ${HMF_BwaMem2Index} \
        ${FastqDir}/${SeqID}.trimmed_R1.fastq.gz \
        ${FastqDir}/${SeqID}.trimmed_R2.fastq.gz \
        | samtools sort -@ ${SortThreads} -o ${BamDir}/${SeqID}.hmftools.bwa.mem2.bam - 

    samtools index -@ ${Threads} ${BamDir}/${SeqID}.hmftools.bwa.mem2.bam 


    ${HMF_RUN_REDUX} \
        -sample ${SeqID} \
        -input_bam ${BamDir}/${SeqID}.hmftools.bwa.mem2.bam \
        -output_bam ${BamDir}/${SeqID}.redux.bam \
        -output_dir ${BamDir} \
        -ref_genome_version 37 -ref_genome ${HMF_GenomeFasta} \
        -unmap_regions ${HMF_UNMAP_REGION} -ref_genome_msi_file ${HMF_MSI_JITTER_SITES} \
        -form_consensus -write_stats -bamtool /opt/conda/bin/samtools \
        -threads ${Threads}

#-------------------------------------------------------------------------------




