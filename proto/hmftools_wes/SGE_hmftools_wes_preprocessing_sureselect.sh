

#! /bin/bash

#$ -N hmftools_wes_pon
#$ -S /bin/bash
#$ -pe mpi 16
#$ -t 1-20
#$ -tc 5
#$ -e /data/pon/WES-panel.of.normals-hg19-sureselect.v6/meta/log/sge.hmftools.wes.pon.preprocess.$TASK_ID.err
#$ -o /data/pon/WES-panel.of.normals-hg19-sureselect.v6/meta/log/sge.hmftools.wes.pon.preprocess.$TASK_ID.err

#---/ Manual Params |-----------------------------------------------------------
    BaseDir="/data/pon"
    BatchID="WES-panel.of.normals-hg19-sureselect.v6"
    SampleListFile="/data/pon/WES-panel.of.normals-hg19-sureselect.v6/meta/sample.list"
    Threads=16
    GenomeAssembly="hg19"
    WesLibrary="sureselect.v6"
#-------------------------------------------------------------------------------

#---/ Sample-Info from Sample-List-File |---------------------------------------
    # SeqID  ReadGroupID SampleName
    
    RunInfo=$(sed -n -e "${SGE_TASK_ID} p" ${SampleListFile})
    
    SeqID=$(echo ${RunInfo} | awk '{print $1}')
    RGID=$(echo ${RunInfo} | awk '{print $2}')
    RGSM=$(echo ${RunInfo} | awk '{print $1}')

#-------------------------------------------------------------------------------

#---| HMFTOOLS |----------------------------------------------------------------

    RUN_SIF="singularity exec -B /storage,/data /storage/images"
    RUN_BWA_MEM2="${RUN_SIF}/bwa-mem2-v2.2.1.sif bwa-mem2"
    HMF_RUN_REDUX="${RUN_SIF}/hmftools-redux-v1.2.3.sif redux"
    HMF_RUN_SAGE="${RUN_SIF}/hmftools-sage-v4.2.sif sage"

    if [ "${GenomeAssembly}" = "hg19" ]; then
        HMF_GENOME_FASTA="/storage/references_and_index/hg19/hmftools/Homo_sapiens.GRCh37.GATK.illumina.fasta"
        HMF_BWAMEM2_INDEX="/storage/references_and_index/hg19/hmftools/bwa-mem2_index-2.2.1/Homo_sapiens.GRCh37.GATK.illumina.fasta"
        GenomeVersion="37"
    else
        HMF_GENOME_FASTA="/storage/references_and_index/hg38/hmftools/GRCh38_masked_exclusions_alts_hlas.fasta"
        HMF_BWAMEM2_INDEX="/storage/references_and_index/hg38/hmftools/bwa-mem2_index-2.2.1/GRCh38_masked_exclusions_alts_hlas.fasta"
        GenomeVersion="38"
    fi

    HMF_RESOURCE="/storage/references_and_index/${GenomeAssembly}/hmftools/resources_links"
    HMF_UNMAP_REGION=${HMF_RESOURCE}/UNMAP-REGIONS.tsv
    HMF_MSI_JITTER_SITES=${HMF_RESOURCE}/MSI-JITTER-SITES.tsv.gz
    HMF_HIGH_CONFIDENCE_BED=${HMF_RESOURCE}/HIGH-CONFIDENCE.bed.gz
    HMF_ENSEMBL_DIR=${HMF_RESOURCE}/ENSEMBL-DATA/
    #HMF_PON=${HMF_RESOURCE}/PANEL-OF-NORMALS.tsv.gz

    TARGET_BED=/storage/references_and_index/hg19/hmftools/resources_wes_${WesLibrary}/${GenomeAssembly}.interval.merged.target.bed

#-------------------------------------------------------------------------------

    SortThreads=4
    BwaThreads=12

    FastqDir=${BaseDir}/${BatchID}/${SeqID}/fastq
    
    #WES_BATCH=$(echo ${SeqID} | sed 's/_0[1-5]$//')
    #FastqDir="/data/wes/${WES_BATCH}/${SeqID}/fastq"

    BamDir=${BaseDir}/${BatchID}/${SeqID}/bam
    HmftoolsDir=${BaseDir}/${BatchID}/${SeqID}/hmftools
    mkdir -p ${HmftoolsDir}
    
    # bwa mem2    
    ${RUN_BWA_MEM2} mem -M -Y -q -K 100000000 -t ${BwaThreads} \
        -R "@RG\tID:${RGID}\tSM:${RGSM}\tLB:TWISTEXOME2\tPL:ILLUMINA\tCN:GCX" \
        ${HMF_BWAMEM2_INDEX} \
        ${FastqDir}/${SeqID}.trimmed_R1.fastq.gz \
        ${FastqDir}/${SeqID}.trimmed_R2.fastq.gz \
        | samtools sort -@ ${SortThreads} -o ${BamDir}/${SeqID}.hmftools.bwa.mem2.bam - 

    samtools index -@ ${Threads} ${BamDir}/${SeqID}.hmftools.bwa.mem2.bam 

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

    # sage germline mode
    SageOutputDir=${HmftoolsDir}/sage
    mkdir -p ${SageOutputDir}

    ${HMF_RUN_SAGE} \
        -tumor ${SeqID} \
        -tumor_bam ${ReduxOutputDir}/${SeqID}.redux.bam \
        -output_vcf ${SageOutputDir}/${SeqID}.sage.germline.vcf.gz \
        -ref_genome_version ${GenomeVersion} -ref_genome ${HMF_GENOME_FASTA} \
        -high_confidence_bed ${HMF_HIGH_CONFIDENCE_BED} \
        -ensembl_data_dir ${HMF_ENSEMBL_DIR} \
        -germline \
        -panel_bed ${TARGET_BED} \
        -threads ${Threads}

#-------------------------------------------------------------------------------

