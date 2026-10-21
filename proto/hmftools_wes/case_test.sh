
BaseDir="/data/wes"
BatchID="WES_25_09-WES_25_10_revised"
TumorSeqID="WES_25_09_09"
OrganoidSeqID="WES_25_09_21"
NormalSeqID="WES_25_09_32"
Threads=24
GenomeAssembly="hg19"
WesLibrary="twist.exome.2.0"

### RUN HMFTOOLS
    RUN_SIF="singularity exec -B /storage,/data /storage/images"
    RUN_BWA_MEM2="${RUN_SIF}/bwa-mem2-v2.2.1.sif bwa-mem2"
    HMF_RUN_REDUX="${RUN_SIF}/hmftools-redux-v1.2.3.sif redux"
    HMF_RUN_SAGE="${RUN_SIF}/hmftools-sage-v4.2.sif sage"
    HMF_RUN_PAVE="${RUN_SIF}/hmftools-pave-v1.8.2.sif pave"
    HMF_RUN_AMBER="${RUN_SIF}/hmftools-amber-v4.2.sif amber"
    HMF_RUN_COBALT="${RUN_SIF}/hmftools-cobalt-v2.2.sif cobalt"
    HMF_RUN_PURPLE="${RUN_SIF}/hmftools-purple-v4.3.sif purple"
    SAGE_JAR="/storage/apps/hmftools/sage_v4.2.jar"
###
### HMFTOOLS RESOURCES 
    HMF_GENOME_FASTA="/storage/references_and_index/hg19/hmftools/Homo_sapiens.GRCh37.GATK.illumina.fasta"
    HMF_BWAMEM2_INDEX="/storage/references_and_index/hg19/hmftools/bwa-mem2_index-2.2.1/Homo_sapiens.GRCh37.GATK.illumina.fasta"
    HMF_GENOME=37

    HMF_WES_RESOURCE="/storage/references_and_index/${GenomeAssembly}/hmftools/resources_wes_${WesLibrary}"
    HMF_WES_PON=${HMF_WES_RESOURCE}/WES_PANEL_OF_NORMALS.tsv.gz
    HMF_WES_GNOMAD=${HMF_WES_RESOURCE}/WES_GNOMAD.csv.gz
    HMF_WES_AMBER_LOCI=${HMF_WES_RESOURCE}/WES_AMBER_LOCI.tsv.gz
    HMF_WES_DIPLOID_REGION=${HMF_WES_RESOURCE}/WES_DIPLOID_REGION.bed.gz
    HMF_WES_TARGET_BED=${HMF_WES_RESOURCE}/WES_TARGET_REGION.bed

    HMF_RESOURCE="/storage/references_and_index/${GenomeAssembly}/hmftools/resources_links"
    HMF_UNMAP_REGION=${HMF_RESOURCE}/UNMAP-REGIONS.tsv
    HMF_DRIVER_PANEL=${HMF_RESOURCE}/DRIVER-GENE-PANEL.tsv
    HMF_ENSEMBL_DIR=${HMF_RESOURCE}/ENSEMBL-DATA/
    HMF_MSI_JITTER_SITES=${HMF_RESOURCE}/MSI-JITTER-SITES.tsv.gz
    HMF_SOMATIC_HOTSPOT_VCF=${HMF_RESOURCE}/KNOWN-HOTSPOT-SOMATIC.vcf.gz
    HMF_HIGH_CONFIDENCE_BED=${HMF_RESOURCE}/HIGH-CONFIDENCE.bed.gz
    HMF_COVERAGE_CODINGPANEL=${HMF_RESOURCE}/COVERAGE-CODING-PANEL.bed.gz
    HMF_MAPPABILITY=${HMF_RESOURCE}/MAPPABILITY.bed.gz
    HMF_CLINVAR_VCF=${HMF_RESOURCE}/CLINVAR.vcf.gz
    HMF_PANEL_ARTIFACT=${HMF_RESOURCE}/PANEL-ARTIFACTS.tsv
    HMF_GERMLINE_BLACKLIST_BED=${HMF_RESOURCE}/GERMLINE-BLACKLIST-REGION.bed
    HMF_GERMLINE_BLACKLIST_VCF=${HMF_RESOURCE}/GERMLINE-BLACKLIST-VARS.vcf.gz
    HMF_GERMLINE_HOTSPOT_VCF=${HMF_RESOURCE}/KNOWN-HOTSPOT-GERMLINE.vcf.gz

    HMF_GC_PROFILE=${HMF_RESOURCE}/GC-PROFILE.cnp
    HMF_GERMLINE_DEL_FREQ=${HMF_RESOURCE}/COHORT-GERMLINE-DELETION-FREQ.csv

###

### 
    SortThreads=6
    BwaThreads=18
    for SID in "WES_25_09_09" "WES_25_09_21" "WES_25_09_32"
    do
        ${RUN_BWA_MEM2} mem -M -Y -q -K 100000000 -t ${BwaThreads} \
            -R "@RG\tID:WES_HMFTOOLS\tSM:${SID}\tLB:TWISTEXOME2\tPL:ILLUMINA\tCN:GCX" \
            ${HMF_BWAMEM2_INDEX} \
            ${BaseDir}/${BatchID}/${SID}/fastq/${SID}.trimmed_R1.fastq.gz \
            ${BaseDir}/${BatchID}/${SID}/fastq/${SID}.trimmed_R2.fastq.gz \
            | samtools sort -@ ${SortThreads} -o ${BaseDir}/${BatchID}/${SID}/bam/${SID}.bwa.mem2.bam - 
        samtools index -@ ${Threads} ${BaseDir}/${BatchID}/${SID}/bam/${SID}.bwa.mem2.bam 
    done
    
    for SID in "WES_25_09_09" "WES_25_09_21" "WES_25_09_32"
    do
        REDUX_OUTPUT=${BaseDir}/${BatchID}/${SID}/hmftools/redux
        mkdir -p ${REDUX_OUTPUT} 

        ${HMF_RUN_REDUX} \
            -sample ${SID} -input_bam ${BaseDir}/${BatchID}/${SID}/bam/${SID}.bwa.mem2.bam \
            -output_bam ${REDUX_OUTPUT}/${SID}.redux.bam \
            -output_dir ${REDUX_OUTPUT} \
            -ref_genome_version ${HMF_GENOME} -ref_genome ${HMF_GENOME_FASTA} \
            -unmap_regions ${HMF_UNMAP_REGION} -ref_genome_msi_file ${HMF_MSI_JITTER_SITES} \
            -form_consensus -write_stats -bamtool /opt/conda/bin/samtools \
            -threads ${Threads} 
    done

    cp ${BaseDir}/${BatchID}/WES_25_09_32/hmftools/redux/* ${BaseDir}/${BatchID}/WES_25_09_09/hmftools/redux/
    cp ${BaseDir}/${BatchID}/WES_25_09_32/hmftools/redux/* ${BaseDir}/${BatchID}/WES_25_09_21/hmftools/redux/



    TumorSeqID="WES_25_09_21"
    NormalSeqID="WES_25_09_32"
        
    REDUX_OUTPUT=${BaseDir}/${BatchID}/${TumorSeqID}/hmftools/redux

### MATCHED-NORMAL    
    # SAGE
    SAGE_OUTPUT_NT=${BaseDir}/${BatchID}/${TumorSeqID}/hmftools/sage_matched_normal
    mkdir -p ${SAGE_OUTPUT_NT}
        
    ${HMF_RUN_SAGE} \
        -tumor ${TumorSeqID} -reference ${NormalSeqID} \
        -tumor_bam ${REDUX_OUTPUT}/${TumorSeqID}.redux.bam \
        -reference_bam ${REDUX_OUTPUT}/${NormalSeqID}.redux.bam \
        -output_vcf ${SAGE_OUTPUT_NT}/${TumorSeqID}.sage.NT.vcf.gz \
        -jitter_param_dir ${REDUX_OUTPUT} \
        -ref_genome_version ${HMF_GENOME} -ref_genome ${HMF_GENOME_FASTA} \
        -hotspots ${HMF_SOMATIC_HOTSPOT_VCF} -high_confidence_bed ${HMF_HIGH_CONFIDENCE_BED} \
        -ensembl_data_dir ${HMF_ENSEMBL_DIR} -gnomad_freq_file ${HMF_WES_GNOMAD} \
        -driver_gene_panel ${HMF_DRIVER_PANEL} -pon_file ${HMF_WES_PON} \
        -panel_bed ${HMF_WES_TARGET_BED} \
        -threads ${Threads} 

    # PAVE
    PAVE_OUTPUT_NT=${BaseDir}/${BatchID}/${TumorSeqID}/hmftools/pave_matched_normal
    mkdir -p ${PAVE_OUTPUT_NT}

    SAGE_VCF=${SAGE_OUTPUT_NT}/${TumorSeqID}.sage.NT.vcf.gz
    PAVE_VCF=${PAVE_OUTPUT_NT}/${TumorSeqID}.pave.NT.vcf.gz

    ${HMF_RUN_PAVE} \
        -sample ${TumorSeqID} -input_vcf ${SAGE_VCF} \
        -output_dir ${PAVE_OUTPUT_NT} -output_vcf ${PAVE_VCF} \
        -ref_genome_version ${HMF_GENOME} -ref_genome ${HMF_GENOME_FASTA} \
        -clinvar_vcf ${HMF_CLINVAR_VCF} -driver_gene_panel ${HMF_DRIVER_PANEL} \
        -ensembl_data_dir ${HMF_ENSEMBL_DIR} \
        -gnomad_freq_file ${HMF_WES_GNOMAD} \
        -mappability_bed ${HMF_MAPPABILITY} \
        -pon_artefact_file ${HMF_PANEL_ARTIFACT} \
        -pon_file ${HMF_WES_PON} \
        -blacklist_bed ${HMF_GERMLINE_BLACKLIST_BED} -blacklist_vcf ${HMF_GERMLINE_BLACKLIST_VCF} \
        -only_canonical -write_detailed -write_transcript_data \
        -threads ${Threads} 

    # AMBER
    AMBER_OUTPUT_NT=${BaseDir}/${BatchID}/${TumorSeqID}/hmftools/amber_matched_normal
    mkdir -p ${AMBER_OUTPUT_NT}

    ${HMF_RUN_AMBER} \
        -tumor ${TumorSeqID} \
        -tumor_bam ${REDUX_OUTPUT}/${TumorSeqID}.redux.bam \
        -reference ${NormalSeqID} \
        -reference_bam ${REDUX_OUTPUT}/${NormalSeqID}.redux.bam \
        -output_dir ${AMBER_OUTPUT_NT} \
        -ref_genome_version ${HMF_GENOME} -ref_genome ${HMF_GENOME_FASTA} \
        -loci ${HMF_WES_AMBER_LOCI} -threads ${Threads} 

    # COBALT
    COBALT_OUTPUT_NT=${BaseDir}/${BatchID}/${TumorSeqID}/hmftools/cobalt_matched_normal
    mkdir -p ${COBALT_OUTPUT_NT}

    ${HMF_RUN_COBALT} \
        -tumor ${TumorSeqID} -tumor_bam ${REDUX_OUTPUT}/${TumorSeqID}.redux.bam \
        -reference ${NormalSeqID} -reference_bam ${REDUX_OUTPUT}/${NormalSeqID}.redux.bam \
        -output_dir ${COBALT_OUTPUT_NT} \
        -gc_profile ${HMF_GC_PROFILE} \
        -ref_genome_version ${HMF_GENOME} -ref_genome ${HMF_GENOME_FASTA} \
        -threads ${Threads} 

    # PURPLE
    PURPLE_OUTPUT_NT=${BaseDir}/${BatchID}/${TumorSeqID}/hmftools/purple_matched_normal
    mkdir -p ${PURPLE_OUTPUT_NT}

    ${HMF_RUN_PURPLE} \
        -tumor ${TumorSeqID} -reference ${NormalSeqID} \
        -amber ${AMBER_OUTPUT_NT} \
        -cobalt ${COBALT_OUTPUT_NT} \
        -somatic_vcf ${PAVE_OUTPUT_NT}/${TumorSeqID}.pave.NT.vcf.gz \
        -ref_genome_version ${HMF_GENOME} -ref_genome ${HMF_GENOME_FASTA} \
        -gc_profile ${HMF_GC_PROFILE} \
        -ensembl_data_dir ${HMF_ENSEMBL_DIR} \
        -driver_gene_panel ${HMF_DRIVER_PANEL} \
        -somatic_hotspots ${HMF_SOMATIC_HOTSPOT_VCF} \
        -germline_hotspots ${HMF_GERMLINE_HOTSPOT_VCF} \
        -germline_del_freq_file ${HMF_GERMLINE_DEL_FREQ} \
        -circos "/storage/apps/circos-0.69-10/bin/circos" \
        -target_regions_bed ${HMF_WES_TARGET_BED} \
        -output_dir ${PURPLE_OUTPUT_NT} \
        -threads ${Threads} 
###
#---------------------------------------------------------------------------------------------------
### TUMOR-ONLY
    # SAGE
    SAGE_OUTPUT_TONLY=${BaseDir}/${BatchID}/${TumorSeqID}/hmftools/sage_tumor_only
    mkdir -p ${SAGE_OUTPUT_TONLY}
        
    ${HMF_RUN_SAGE} \
        -tumor ${TumorSeqID} \
        -tumor_bam ${REDUX_OUTPUT}/${TumorSeqID}.redux.bam \
        -output_vcf ${SAGE_OUTPUT_TONLY}/${TumorSeqID}.sage.Tonly.vcf.gz \
        -jitter_param_dir ${REDUX_OUTPUT} \
        -ref_genome_version ${HMF_GENOME} -ref_genome ${HMF_GENOME_FASTA} \
        -hotspots ${HMF_SOMATIC_HOTSPOT_VCF} -high_confidence_bed ${HMF_HIGH_CONFIDENCE_BED} \
        -ensembl_data_dir ${HMF_ENSEMBL_DIR} -gnomad_freq_file ${HMF_WES_GNOMAD} \
        -driver_gene_panel ${HMF_DRIVER_PANEL} -pon_file ${HMF_WES_PON} \
        -panel_bed ${HMF_WES_TARGET_BED} \
        -threads ${Threads} 

    # PAVE
    PAVE_OUTPUT_TONLY=${BaseDir}/${BatchID}/${TumorSeqID}/hmftools/pave_tumor_only
    mkdir -p ${PAVE_OUTPUT_TONLY}

    SAGE_VCF=${SAGE_OUTPUT_TONLY}/${TumorSeqID}.sage.Tonly.vcf.gz
    PAVE_VCF=${PAVE_OUTPUT_TONLY}/${TumorSeqID}.pave.Tonly.vcf.gz

    ${HMF_RUN_PAVE} \
        -sample ${TumorSeqID} -input_vcf ${SAGE_VCF} \
        -output_dir ${PAVE_OUTPUT_TONLY} -output_vcf ${PAVE_VCF} \
        -ref_genome_version ${HMF_GENOME} -ref_genome ${HMF_GENOME_FASTA} \
        -clinvar_vcf ${HMF_CLINVAR_VCF} -driver_gene_panel ${HMF_DRIVER_PANEL} \
        -ensembl_data_dir ${HMF_ENSEMBL_DIR} \
        -gnomad_freq_file ${HMF_WES_GNOMAD} \
        -mappability_bed ${HMF_MAPPABILITY} \
        -pon_artefact_file ${HMF_PANEL_ARTIFACT} \
        -pon_file ${HMF_WES_PON} \
        -blacklist_bed ${HMF_GERMLINE_BLACKLIST_BED} -blacklist_vcf ${HMF_GERMLINE_BLACKLIST_VCF} \
        -only_canonical -write_detailed -write_transcript_data \
        -threads ${Threads} 

    # AMBER
    AMBER_OUTPUT_TONLY=${BaseDir}/${BatchID}/${TumorSeqID}/hmftools/amber_tumor_only
    mkdir -p ${AMBER_OUTPUT_TONLY}

    ${HMF_RUN_AMBER} \
        -tumor ${TumorSeqID} \
        -tumor_bam ${REDUX_OUTPUT}/${TumorSeqID}.redux.bam \
        -output_dir ${AMBER_OUTPUT_TONLY} \
        -ref_genome_version ${HMF_GENOME} -ref_genome ${HMF_GENOME_FASTA} \
        -loci ${HMF_WES_AMBER_LOCI} -threads ${Threads} 

    # COBALT
    COBALT_OUTPUT_TONLY=${BaseDir}/${BatchID}/${TumorSeqID}/hmftools/cobalt_tumor_only
    mkdir -p ${COBALT_OUTPUT_TONLY}

    ${HMF_RUN_COBALT} \
        -tumor ${TumorSeqID} -tumor_bam ${REDUX_OUTPUT}/${TumorSeqID}.redux.bam \
        -output_dir ${COBALT_OUTPUT_TONLY} \
        -gc_profile ${HMF_GC_PROFILE} \
        -tumor_only_diploid_bed ${HMF_WES_DIPLOID_REGION} \
        -ref_genome_version ${HMF_GENOME} -ref_genome ${HMF_GENOME_FASTA} \
        -threads ${Threads} 

    # PURPLE
    PURPLE_OUTPUT_TONLY=${BaseDir}/${BatchID}/${TumorSeqID}/hmftools/purple_tumor_only
    mkdir -p ${PURPLE_OUTPUT_TONLY}

    ${HMF_RUN_PURPLE} \
        -tumor ${TumorSeqID} \
        -amber ${AMBER_OUTPUT_TONLY} \
        -cobalt ${COBALT_OUTPUT_TONLY} \
        -somatic_vcf ${PAVE_OUTPUT_TONLY}/${TumorSeqID}.pave.Tonly.vcf.gz \
        -ref_genome_version ${HMF_GENOME} -ref_genome ${HMF_GENOME_FASTA} \
        -gc_profile ${HMF_GC_PROFILE} \
        -ensembl_data_dir ${HMF_ENSEMBL_DIR} \
        -driver_gene_panel ${HMF_DRIVER_PANEL} \
        -somatic_hotspots ${HMF_SOMATIC_HOTSPOT_VCF} \
        -germline_hotspots ${HMF_GERMLINE_HOTSPOT_VCF} \
        -germline_del_freq_file ${HMF_GERMLINE_DEL_FREQ} \
        -circos "/storage/apps/circos-0.69-10/bin/circos" \
        -target_regions_bed ${HMF_WES_TARGET_BED} \
        -output_dir ${PURPLE_OUTPUT_TONLY} \
        -threads ${Threads} 

###




