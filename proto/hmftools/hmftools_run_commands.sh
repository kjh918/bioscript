


### RUN HMFTOOLS

    HMF_RUN_REDUX="singularity exec -B /storage,/data /storage/images/hmftools-redux-v1.2.3.sif redux"
    HMF_RUN_SAGE="singularity exec -B /storage,/data /storage/images/hmftools-sage-v4.2.sif sage"
    HMF_RUN_PAVE="singularity exec -B /storage,/data /storage/images/hmftools-pave-v1.8.2.sif pave"
    HMF_RUN_AMBER="singularity exec -B /storage,/data /storage/images/hmftools-amber-v4.2.sif amber"
    HMF_RUN_COBALT="singularity exec -B /storage,/data /storage/images/hmftools-cobalt-v2.2.sif cobalt"
    HMF_RUN_ESVEE="singularity exec -B /storage,/data /storage/images/hmftools-esvee-v1.2.sif esvee"
    HMF_RUN_PURPLE="singularity exec -B /storage,/data /storage/images/hmftools-purple-v4.3.sif purple"
    HMF_RUN_LINX="singularity exec -B /storage,/data /storage/images/hmftools-linx-v2.2.sif linx"


### RESOURCES
    HMF_GENOME_FASTA="/storage/references_and_index/hg38/hmftools/GRCh38_masked_exclusions_alts_hlas.fasta"
    HMF_RESOURCE="/storage/references_and_index/hg38/hmftools/hmf_pipeline_resources.38_v2.3.0--2"

    HMF_UNMAP_REGION=${HMF_RESOURCE}/common/unmap_regions.38.tsv
    HMF_MSI_JITTER_SITES=${HMF_RESOURCE}/dna/variants/msi_jitter_sites.38.tsv.gz
    HMF_SOMATIC_HOTSPOT_VCF=${HMF_RESOURCE}/dna/variants/KnownHotspots.somatic.38.vcf.gz
    HMF_HIGH_CONFIDENCE_BED=${HMF_RESOURCE}/dna/variants/HighConfidence.bed.gz
    HMF_DRIVER_PANEL=${HMF_RESOURCE}/common/DriverGenePanel.38.tsv
    HMF_PON=${HMF_RESOURCE}/dna/variants/hmf_wgs_sage_pon_1000.38.tsv.gz
    HMF_ENSEMBL_DIR=${HMF_RESOURCE}/common/ensembl_data
    HMF_GNOMAD_DIR=${HMF_RESOURCE}/dna/variants/gnomad
    HMF_MAPPABILITY=${HMF_RESOURCE}/dna/variants/mappability_150.38.bed.gz
    HMF_CLINVAR_VCF=${HMF_RESOURCE}/dna/variants/clinvar.38.vcf.gz
    HMF_PANEL_ARTIFACT=${HMF_RESOURCE}/dna/variants/PanelArtefacts.38.tsv
    HMF_GERMLINE_BLACKLIST_BED=${HMF_RESOURCE}/dna/variants/KnownBlacklist.germline.38.bed
    HMF_GERMLINE_BLACKLIST_VCF=${HMF_RESOURCE}/dna/variants/KnownBlacklist.germline.38.vcf.gz
    HMF_AMBER_LOCI=${HMF_RESOURCE}/dna/copy_number/AmberGermlineSites.38.tsv.gz
    HMF_GC_PROFILE=${HMF_RESOURCE}/dna/copy_number/GC_profile.1000bp.38.cnp
    HMF_DIPLOID_BED=${HMF_RESOURCE}/dna/copy_number/DiploidRegions.38.bed.gz
    HMF_PON_SV=${HMF_RESOURCE}/dna/sv/sv_pon.38.bedpe.gz
    HMF_PON_SGL=${HMF_RESOURCE}/dna/sv/sgl_pon.38.bed.gz
    HMF_REPEAT_MASK=${HMF_RESOURCE}/dna/sv/repeat_mask_data.38.fa.gz
    HMF_SV_BLACKLIST_BED=${HMF_RESOURCE}/dna/sv/sv_prep_blacklist.38.bed
    HMF_KNOWN_FUSION_BED=${HMF_RESOURCE}/known_fusions.38.bedpe
    HMF_GERMLINE_DEL_FREQ=${HMF_RESOURCE}/dna/copy_number/cohort_germline_del_freq.38.csv
    HMF_GERMLINE_HOTSPOT_VCF=${HMF_RESOURCE}/dna/variants/KnownHotspots.germline.38.vcf.gz
    HMF_LILAC_RESOURCE_DIR=${HMF_RESOURCE}/misc/lilac

#---------------------------------------------------------------------------------------------------
#  REDUX 
#---------------------------------------------------------------------------------------------------
    # SIF : /storage/images/hmftools-redux-v1.2.3.sif
    HMF_RUN_REDUX="singularity exec -B /storage,/data /storage/images/hmftools-redux-v1.2.3.sif redux"

    TumorSeqID=
    TumorBamFile=
    OutputDir=
    OutputBAM=
    Threads=

    ${HMF_RUN_REDUX} \
        -sample ${TumorSeqID} \
        -input_bam ${TumorBamFile} \
        -output_dir ${OutputDir} \
        -output_bam ${OutputBAM} \
        -ref_genome_version 38 -ref_genome ${HMF_GENOME_FASTA} \
        -unmap_regions ${HMF_UNMAP_REGION} -ref_genome_msi_file ${HMF_MSI_JITTER_SITES} \
        -form_consensus -write_stats -bamtool /opt/conda/bin/samtools \
        -threads ${Threads} 
#---------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------
#  SAGE 
#---------------------------------------------------------------------------------------------------
    # SIF : /storage/images/hmftools-redux-v1.2.3.sif
    HMF_RUN_SAGE="singularity exec -B /storage,/data /storage/images/hmftools-sage-v4.2.sif sage"

    TumorSeqID=
    TumorBamFile=
    NormalSeqID=
    NormalBamFile=
    OutputVCF=
    Threads=
    ReduxJitterParamDir=

    ${HMF_RUN_SAGE} \
        -tumor ${TumorSeqID} -tumor_bam ${TumorBamFile} \
        -reference ${NormalSeqID} -reference_bam ${NormalBamFile} \
        -output_vcf ${OutputVCF} \
        -jitter_param_dir ${ReduxJitterParamDir} \
        -ref_genome_version 38 -ref_genome ${HMF_GENOME_FASTA} \
        -hotspots ${HMF_SOMATIC_HOTSPOT_VCF} -high_confidence_bed ${HMF_HIGH_CONFIDENCE_BED} \
        -ensembl_data_dir ${HMF_ENSEMBL_DIR} -gnomad_freq_dir ${HMF_GNOMAD_DIR} \
        -driver_gene_panel ${HMF_DRIVER_PANEL} -pon_file ${HMF_PON} \
        -threads ${Threads}
#---------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------
#  PAVE
#---------------------------------------------------------------------------------------------------
    # SIF : /storage/images/hmftools-redux-v1.2.3.sif
    HMF_RUN_PAVE="singularity exec -B /storage,/data /storage/images/hmftools-pave-v1.8.2.sif pave"

    TumorSeqID=
    InputVCF=
    OutputDir=
    OutputVCF=
    Threads=

    ${HMF_RUN_PAVE} \
        -sample ${TumorSeqID} \
        -input_vcf ${InputVCF} \
        -output_dir ${OutputDir} \
        -output_vcf ${OutputVCF} \
        -ref_genome_version 38 -ref_genome ${HMF_GENOME_FASTA} \
        -clinvar_vcf ${HMF_CLINVAR_VCF} -driver_gene_panel ${HMF_DRIVER_PANEL} \
        -ensembl_data_dir ${HMF_ENSEMBL_DIR} -gnomad_freq_dir ${HMF_GNOMAD_DIR} \
        -mappability_bed ${HMF_MAPPABILITY} -pon_artefact_file ${HMF_PANEL_ARTIFACT} -pon_file ${HMF_PON} \
        -blacklist_bed ${HMF_GERMLINE_BLACKLIST_BED} -blacklist_vcf ${HMF_GERMLINE_BLACKLIST_VCF} \
        -only_canonical -threads ${Threads} 
#---------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------
#  AMBER
#---------------------------------------------------------------------------------------------------
    # SIF : /storage/images/hmftools-amber-v4.2.sif
    HMF_RUN_AMBER="singularity exec -B /storage,/data /storage/images/hmftools-amber-v4.2.sif amber"

    TumorSeqID=
    TumorBamFile=
    NormalSeqID=
    NormalBamFile=
    OutputDir=
    Threads=

    ${HMF_RUN_AMBER} \
        -tumor ${TumorSeqID} \
        -tumor_bam ${TumorBamFile} \
        -reference ${NormalSeqID} \
        -reference_bam ${NormalBamFile} \
        -output_dir ${OutputDir} \
        -ref_genome_version 38 -ref_genome ${HMF_GENOME_FASTA} \
        -loci ${HMF_AMBER_LOCI} -threads ${Threads}
#---------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------
#  COBALT
#---------------------------------------------------------------------------------------------------
    # SIF : /storage/images/hmftools-cobalt-v2.2.sif
    HMF_RUN_COBALT="singularity exec -B /storage,/data /storage/images/hmftools-cobalt-v2.2.sif cobalt"
    
    TumorSeqID=
    TumorBamFile=
    NormalSeqID=
    NormalBamFile=
    OutputDir=
    Threads=

    ${HMF_RUN_COBALT} \
        -tumor ${TumorSeqID} \
        -tumor_bam ${TumorBamFile} \
        -reference ${NormalSeqID} \
        -reference_bam ${NormalBamFile} \
        -output_dir ${OutputDir} \
        -gc_profile ${HMF_GC_PROFILE} \
        -ref_genome_version 38 -ref_genome ${HMF_GENOME_FASTA} \
#---------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------
#  ESVEE
#---------------------------------------------------------------------------------------------------
    # SIF : /storage/images/hmftools-esvee-v1.2.sif
    HMF_RUN_ESVEE="singularity exec -B /storage,/data /storage/images/hmftools-esvee-v1.2.sif esvee"
    TumorSeqID=
    TumorBamFile=
    NormalSeqID=
    NormalBamFile=
    OutputDir=
    Threads=

    ${HMF_RUN_ESVEE} \
        -tumor ${TumorSeqID} \
        -tumor_bam ${TumorBamFile} \
        -reference ${NormalSeqID} \
        -reference_bam ${NormalBamFile} \
        -output_dir ${OutputDir} \
        -ref_genome_version 38 -ref_genome ${HMF_GENOME_FASTA} \
        -known_hotspot_file ${HMF_KNOWN_FUSION_BED} \
        -pon_sgl_file ${HMF_PON_SGL} \
        -pon_sv_file ${HMF_PON_SV} \
        -repeat_mask_file ${HMF_REPEAT_MASK} \
        -blacklist_bed ${HMF_SV_BLACKLIST_BED} \
        -unmap_regions ${HMF_UNMAP_REGION} \
        -threads ${Threads}
#---------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------
#  PURPLE
#---------------------------------------------------------------------------------------------------
    # SIF : /storage/images/hmftools-purple-v4.3.sif
    HMF_RUN_PURPLE="singularity exec -B /storage,/data /storage/images/hmftools-purple-v4.3.sif purple"

    AmberResultDir=
    CobaltResultDir=
    PaveResultVCF=
    EsveeResultVCF=

    TumorSeqID=
    TumorBamFile=
    NormalSeqID=
    NormalBamFile=
    OutputDir=
    Threads=

    ${HMF_RUN_PURPLE} \
        -tumor ${TumorSeqID} \
        -reference ${NormalSeqID} \
        -output_dir ${OUTDIR} \
        -ref_genome_version 38 -ref_genome ${HMF_GENOME_FASTA} \
        -amber ${AmberResultDir} \
        -cobalt ${CobaltResultDir} \
        -somatic_vcf ${PaveResultVCF} \
        -somatic_sv_vcf ${EsveeResultVCF} \
        -ensembl_data_dir ${HMF_ENSEMBL_DIR} \
        -driver_gene_panel ${HMF_DRIVER_PANEL} \
        -gc_profile ${HMF_GC_PROFILE} \
        -germline_del_freq_file ${HMF_GERMLINE_DEL_FREQ} \
        -somatic_hotspots ${HMF_SOMATIC_HOTSPOT_VCF} \
        -germline_hotspots ${HMF_GERMLINE_HOTSPOT_VCF} \
        -threads ${Threads}
#---------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------
#  LINX
#---------------------------------------------------------------------------------------------------
    # SIF : /storage/images/hmftools-linx-v2.2.sif
    HMF_RUN_LINX="singularity exec -B /storage,/data /storage/images/hmftools-linx-v2.2.sif linx"

    TumorSeqID=
    TumorBamFile=
    NormalSeqID=
    NormalBamFile=
    OutputDir=
    Threads=
    PurpleDir=

    ${HMF_RUN_LINX} \
        -sample ${TumorSeqID} \
        -purple_dir ${PurpleDir} \
        -output_dir ${OutputDir} \
        -ref_genome_version 38 -ref_genome ${HMF_GENOME_FASTA} \
        -ensembl_data_dir ${HMF_ENSEMBL_DIR} \
        -driver_gene_panel ${HMF_DRIVER_PANEL} \
        -known_fusion_file ${HMF_KNOWN_FUSION_BED} \
        -threads ${Threads}
#---------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------
#  CHORD
#---------------------------------------------------------------------------------------------------
    # SIF : /storage/images/hmftools-chord-v2.1.2.sif
    HMF_RUN_CHORD="singularity exec -B /storage,/data /storage/images/hmftools-chord-v2.1.2.sif chord"

    TumorSeqID=
    OutputDir=
    Threads=
    PurpleDir=

    ${HMF_RUN_CHORD} \
        -sample ${TumorSeqID} \
        -output_dir ${OutputDir} \
        -purple_dir ${PurpleDir} \
        -ref_genome ${HMF_GENOME_FASTA} \
        -include_non_pass -write_detailed_files \
        -threads ${Threads}
#---------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------
#  LILAC
#---------------------------------------------------------------------------------------------------
    # SIF : /storage/images/hmftools-lilac-v1.7.3.sif
    HMF_RUN_LILAC="singularity exec -B /storage,/data /storage/images/hmftools-lilac-v1.7.3.sif lilac"

    TumorSeqID=
    TumorBamFile=
    NormalSeqID=
    NormalBamFile=
    OutputDir=
    Threads=
    PurpleDir=

    ${HMF_RUN_LILAC} \
        -sample ${TumorSeqID} \
        -tumor_bam ${TumorBamFile} \
        -reference_bam ${NormalBamFile} \
        -output_dir ${OutputDir} \
        -ref_genome_version 38 -ref_genome ${HMF_GENOME_FASTA} \
        -purple_dir ${PURPLE_DIR} \
        -resource_dir ${HMF_LILAC_RESOURCE_DIR} \
        -threads 8 
        #-somatic_vcf ${SOMATIC_VCF} \
        #-gene_copy_number ${GENE_CN_FILE} \
        
        