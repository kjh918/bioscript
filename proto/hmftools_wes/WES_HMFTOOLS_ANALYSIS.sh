


#! /bin/bash

#---| Input Options |------------------------------------------------------------------------------#
    usage() {
        echo "
        Usage: WES_M02-1_MUTECT2-Tumor-Only.sh [ ARGUMENTS ]...
        
        --BASE_DIR          analysis base folder by NGS data type
                            default = /data/wgs
        
        --BATCH_ID          [ REQUIRED ] dataset batch-id
        
        --TUMOR_SEQ_ID      [ REQURIED ] tumor sample-id 

        --NORMAL_SEQ_ID     matched-normal smaple-id
                            Required when 'RUN_MODE' = 'NT'
                
        --THREADS           parallel processing cpu cores 
                            default = 16
        
        --GENOME_ASSEMBLY   Reference genome assembly version
                            default = hghg19

        --WES_LIB           WES capture library 
                            default = twist.exome.2.0

        --RUN_MODE          analysis run mode. 
                            'NT' or 'Tonly'
                            default = 'NT'

        -h | --help ]       Print this usage 
        "
    }
#--------------------------------------------------------------------------------------------------#

#---| Parse Option Arguments |---------------------------------------------------------------------#
    ARGS=$(getopt -a -o h: --long BATCH_ID:,TUMOR_SEQ_ID:,NORMAL_SEQ_ID:,WES_LIB:,BASE_DIR:,THREADS:,GENOME_ASSEMBLY:,NORMAL_TYPE:,RUN_MODE:,help -- "$@" )
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
                -* | --* | "") echo "no value at $1, default = '/data/wes'" ; shift ;;
                *)  BaseDir=$2 ; shift 2 ;;
            esac ;;
        --BATCH_ID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  BatchID=$2 ; shift 2 ;;
            esac ;;
        --TUMOR_SEQ_ID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  TumorSeqID=$2 ; shift 2 ;;
            esac  ;;
        --NORMAL_SEQ_ID )
            case "$2" in
                -* | --* | "") echo "no value at $1" ; shift ;;
                *)  NormalSeqID=$2 ; shift 2 ;;
            esac  ;;
        --WES_LIB )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 'twist.exome.2.0'" ; shift ;;
                *)  WesLibrary=$2 ; shift 2 ;;
            esac ;;
        --THREADS )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 16" ; shift ;;
                *)  Threads=$2 ; shift 2 ;;
            esac ;;
        --GENOME_ASSEMBLY )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = hg19" ; shift ;;
                *)  GenomeAssembly=$2 ; shift 2 ;;
            esac ;;
        --RUN_MODE )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = NT" ; shift ;;
                *)  RunMode=$2 ; shift 2 ;;
            esac ;;
        -h | --help )
        usage >&2 ; exit 2 ;;   
        --) shift ; break ;;
        *)  usage >&2 ; exit 2 ;;
    esac
    done
#--------------------------------------------------------------------------------------------------#

#---| DEFAULT PARAMS |------------------------------------------------------------------------------
    if [[ -z "${WesLibrary}" ]]; then WesLibrary="twist.exome.2.0"; fi
    if [[ -z "${GenomeAssembly}" ]]; then GenomeAssembly="hg19"; fi
    if [[ -z "${Threads}" ]]; then GenomeAssembly="16"; fi
    if [[ "${RunMode}" = "NT" && -z "${NormalSeqID}" ]]; then exit 2; fi
#---------------------------------------------------------------------------------------------------

#---| Manual Params |-------------------------------------------------------------------------------
    BaseDir="/data/wes"
    BatchID="WES_24_02"
    TumorSeqID="WES_24_02_09"
    NormalSeqID="WES_24_02_10"
    WesLibrary="twist.exome.2.0"
    GenomeAssembly="hg19"  # fixed
    RunMode="Tonly" # or "Tonly"
    Threads=16
#---------------------------------------------------------------------------------------------------

###---| RUN HMFTOOLS |------------------------------------------------------------------------------
    RUN_SIF="singularity exec -B /storage,/data /storage/images"
    HMF_RUN_REDUX="${RUN_SIF}/hmftools-redux-v1.2.3.sif redux"
    HMF_RUN_BAMTOOLS="${RUN_SIF}/hmftools-bamtools-v1.5.sif bamtools"
    HMF_RUN_SAGE="${RUN_SIF}/hmftools-sage-v4.2.sif sage"
    HMF_RUN_PAVE="${RUN_SIF}/hmftools-pave-v1.8.2.sif pave"
    HMF_RUN_PAVE_JAR="${RUN_SIF}/hmftools-pave-v1.8.2.sif java -Xmx32G -jar /usr/share/java/pave_v1.8.2.jar"
    HMF_RUN_AMBER="${RUN_SIF}/hmftools-amber-v4.2.sif amber"
    HMF_RUN_COBALT="${RUN_SIF}/hmftools-cobalt-v2.2.sif cobalt"
    HMF_RUN_PURPLE="${RUN_SIF}/hmftools-purple-v4.3.sif java -Xmx48G -jar /usr/share/java/purple_v4.3.jar"   
#---------------------------------------------------------------------------------------------------

###---| RESOURCES |---------------------------------------------------------------------------------
    
    HMF_RESOURCE_BASE="/storage/references_and_index/${GenomeAssembly}/hmftools"
        
    HMF_WES_REOURCE="${HMF_RESOURCE_BASE}/resources_wes_${WesLibrary}"
    HMF_WES_AMBER_LOCI="${HMF_WES_REOURCE}/WES_AMBER_LOCI.tsv.gz"               # <-- HMF_AMBER_LOCI
    HMF_WES_GNOMAD="${HMF_WES_REOURCE}/WES_GNOMAD.csv.gz"                       # <-- GNOMAD_RESOURCE_ARGS
    HMF_WES_PANEL_OF_NORMALS="${HMF_WES_REOURCE}/WES_PANEL_OF_NORMALS.tsv.gz"   # <-- HMF_PON
    HMF_WES_TARGET_BED="${HMF_WES_REOURCE}/WES_TARGET_REGION.bed"
    #HMF_WES_DIPLOID_REGION="${HMF_WES_REOURCE}/WES_DIPLOID_REGION.bed.gz"
    

    if [ "${GenomeAssembly}" = "hg19" ]; then
        HMF_GENOME_FASTA="/storage/references_and_index/hg19/hmftools/Homo_sapiens.GRCh37.GATK.illumina.fasta"
        HMF_BWAMEM2_INDEX="/storage/references_and_index/hg19/hmftools/bwa-mem2_index-2.2.1/Homo_sapiens.GRCh37.GATK.illumina.fasta"
        HMF_GENOME="37"
    else
        HMF_GENOME_FASTA="/storage/references_and_index/hg38/hmftools/GRCh38_masked_exclusions_alts_hlas.fasta"
        HMF_BWAMEM2_INDEX="/storage/references_and_index/hg38/hmftools/bwa-mem2_index-2.2.1/GRCh38_masked_exclusions_alts_hlas.fasta"
        HMF_GENOME="38"
    fi
    
    HMF_COMMON_RESOURCE="${HMF_RESOURCE_BASE}/resources_links"
    HMF_UNMAP_REGION="${HMF_COMMON_RESOURCE}/UNMAP-REGIONS.tsv"
    HMF_MSI_JITTER_SITES="${HMF_COMMON_RESOURCE}/MSI-JITTER-SITES.tsv.gz"
    HMF_DRIVER_PANEL="${HMF_COMMON_RESOURCE}/DRIVER-GENE-PANEL.tsv"
    HMF_ENSEMBL_DIR="${HMF_COMMON_RESOURCE}/ENSEMBL-DATA/"
    HMF_SOMATIC_HOTSPOT_VCF="${HMF_COMMON_RESOURCE}/KNOWN-HOTSPOT-SOMATIC.vcf.gz"
    HMF_GERMLINE_HOTSPOT_VCF="${HMF_COMMON_RESOURCE}/KNOWN-HOTSPOT-GERMLINE.vcf.gz"
    HMF_HIGH_CONFIDENCE_BED="${HMF_COMMON_RESOURCE}/HIGH-CONFIDENCE.bed.gz"
    HMF_CLINVAR_VCF="${HMF_COMMON_RESOURCE}/CLINVAR.vcf.gz"
    HMF_MAPPABILITY="${HMF_COMMON_RESOURCE}/MAPPABILITY.bed.gz"
    HMF_PANEL_ARTIFACT="${HMF_COMMON_RESOURCE}/PANEL-ARTIFACTS.tsv"
    HMF_GERMLINE_BLACKLIST_BED="${HMF_COMMON_RESOURCE}/GERMLINE-BLACKLIST-REGION.bed"
    HMF_GERMLINE_BLACKLIST_VCF="${HMF_COMMON_RESOURCE}/GERMLINE-BLACKLIST-VARS.vcf.gz"
    HMF_GC_PROFILE="${HMF_COMMON_RESOURCE}/GC-PROFILE.cnp"
    HMF_GERMLINE_DEL_FREQ="${HMF_COMMON_RESOURCE}/COHORT-GERMLINE-DELETION-FREQ.csv"
    HMF_DIPLOID_REGIONS="${HMF_COMMON_RESOURCE}/DIPLOID-REGIONS.bed.gz"

#---------------------------------------------------------------------------------------------------

###---| CHECK REDUX-BAM |---------------------------------------------------------------------------

    if [[ ! -f "${BaseDir}/${BatchID}/${TumorSeqID}/hmftools/redux/${TumorSeqID}.redux.bam" ]]; then
        echo "NO TUMOR REDUX BAM."
        exit 2
    fi

    # hmftools result dir
    HMFTOOLS_RES_BASE="${BaseDir}/${BatchID}/${TumorSeqID}/hmftools"
    # run log
    mkdir -p ${HMFTOOLS_RES_BASE}/run_log
    RunLog=${HMFTOOLS_RES_BASE}/run_log/hmftools.analysis.run.mode.${RunMode}.log
    echo " " > "${RunLog}"
    # redux bam check
    OUTPUT_REDUX="${HMFTOOLS_RES_BASE}/redux"
    if [[ "${RunMode}" = "NT" ]]; then
        if [[ ! -f "${OUTPUT_REDUX}/${NormalSeqID}.redux.bam" ]]; then
            if [[ -f "${BaseDir}/${BatchID}/${NormalSeqID}/hmftools/redux/${NormalSeqID}.redux.bam" ]]; then
                cp ${BaseDir}/${BatchID}/${NormalSeqID}/hmftools/redux/* ${OUTPUT_REDUX}/
            else
                echo "No Matched-Normal Sample's REDUX BAM FILE."
                exit 2
            fi
        fi
    fi
    ### redux bam
    TUMOR_REDUX_BAM=${OUTPUT_REDUX}/${TumorSeqID}.redux.bam
    if [[ "${RunMode}" = "NT" ]]; then NORMAL_REDUX_BAM=${OUTPUT_REDUX}/${NormalSeqID}.redux.bam; fi
#---------------------------------------------------------------------------------------------------

###---| SAGE |--------------------------------------------------------------------------------------

    SAGE_INPUT_ARGS=("-tumor ${TumorSeqID}")
    SAGE_INPUT_ARGS+=("-tumor_bam ${TUMOR_REDUX_BAM}")
    if [[ "${RunMode}" == "NT" ]]; then
        SAGE_INPUT_ARGS+=("-reference ${NormalSeqID}")
        SAGE_INPUT_ARGS+=("-reference_bam ${NORMAL_REDUX_BAM}")
    else
        SAGE_INPUT_ARGS+=("-high_depth_mode")
    fi

    # somatic
    if [ "${RunMode}" = "Tonly" ]; then
        OUTPUT_SAGE="${HMFTOOLS_RES_BASE}/sage_Tonly"
    else    
        OUTPUT_SAGE="${HMFTOOLS_RES_BASE}/sage"
    fi   
    mkdir -p "${OUTPUT_SAGE}"

    ${HMF_RUN_SAGE} \
        ${SAGE_INPUT_ARGS[@]} \
        -output_vcf "${OUTPUT_SAGE}/${TumorSeqID}.sage.vcf.gz" \
        -panel_bed "${HMF_WES_TARGET_BED}" \
        -jitter_param_dir "${OUTPUT_REDUX}" \
        -ref_genome_version "${HMF_GENOME}" \
        -ref_genome "${HMF_GENOME_FASTA}" \
        -hotspots "${HMF_SOMATIC_HOTSPOT_VCF}" \
        -high_confidence_bed "${HMF_HIGH_CONFIDENCE_BED}" \
        -ensembl_data_dir "${HMF_ENSEMBL_DIR}" \
        -gnomad_freq_file "${HMF_WES_GNOMAD}" \
        -driver_gene_panel "${HMF_DRIVER_PANEL}" \
        -pon_file "${HMF_WES_PANEL_OF_NORMALS}" \
        -bqr_write_plot \
        -threads "${Threads}" 2>&1 | tee -a "${RunLog}"
    #
    # germline 
    if [ "${RunMode}" = "Tonly" ]; then
        OUTPUT_SAGE_GERMLINE="${HMFTOOLS_RES_BASE}/sage_germline_Tonly"
    else    
        OUTPUT_SAGE_GERMLINE="${HMFTOOLS_RES_BASE}/sage_germline"
    fi   
    mkdir -p "${OUTPUT_SAGE_GERMLINE}"

    ${HMF_RUN_SAGE} -germline \
        ${SAGE_INPUT_ARGS[@]} \
        -output_vcf "${OUTPUT_SAGE_GERMLINE}/${TumorSeqID}.sage.germline.vcf.gz" \
        -jitter_param_dir "${OUTPUT_REDUX}" \
        -ref_genome_version "${HMF_GENOME}" \
        -ref_genome "${HMF_GENOME_FASTA}" \
        -ensembl_data_dir "${HMF_ENSEMBL_DIR}" \
        -driver_gene_panel "${HMF_DRIVER_PANEL}" \
        -gnomad_freq_file "${HMF_WES_GNOMAD}" \
        -bqr_write_plot \
        -threads "${Threads}" 2>&1 | tee -a "${RunLog}"
    #
#---------------------------------------------------------------------------------------------------

###---| PAVE |--------------------------------------------------------------------------------------

    if [ "${RunMode}" = "Tonly" ]; then
        OUTPUT_PAVE="${HMFTOOLS_RES_BASE}/pave_Tonly"
    else    
        OUTPUT_PAVE="${HMFTOOLS_RES_BASE}/pave"
    fi   
    mkdir -p "${OUTPUT_PAVE}"

    # somatic
    ${HMF_RUN_PAVE} \
        -sample "${TumorSeqID}" \
        -input_vcf "${OUTPUT_SAGE}/${TumorSeqID}.sage.vcf.gz" \
        -output_dir "${OUTPUT_PAVE}" \
        -output_vcf "${OUTPUT_PAVE}/${TumorSeqID}.pave.vcf.gz" \
        -ref_genome_version "${HMF_GENOME}" \
        -ref_genome "${HMF_GENOME_FASTA}" \
        -clinvar_vcf "${HMF_CLINVAR_VCF}" \
        -driver_gene_panel "${HMF_DRIVER_PANEL}" \
        -ensembl_data_dir "${HMF_ENSEMBL_DIR}" \
        -gnomad_freq_file "${HMF_WES_GNOMAD}" \
        -mappability_bed "${HMF_MAPPABILITY}" \
        -pon_artefact_file "${HMF_PANEL_ARTIFACT}" \
        -pon_file "${HMF_WES_PANEL_OF_NORMALS}" \
        -blacklist_bed "${HMF_GERMLINE_BLACKLIST_BED}" \
        -blacklist_vcf "${HMF_GERMLINE_BLACKLIST_VCF}" \
        -only_canonical \
        -write_detailed \
        -write_transcript_data \
        -threads "${Threads}" 2>&1 | tee -a "${RunLog}"
    #
    # germline
    ${HMF_RUN_PAVE} \
        -sample "${TumorSeqID}" \
        -input_vcf "${OUTPUT_SAGE_GERMLINE}/${TumorSeqID}.sage.germline.vcf.gz" \
        -output_dir "${OUTPUT_PAVE}" \
        -output_vcf "${OUTPUT_PAVE}/${TumorSeqID}.pave.germline.vcf.gz" \
        -ref_genome_version "${HMF_GENOME}" \
        -ref_genome "${HMF_GENOME_FASTA}" \
        -clinvar_vcf "${HMF_CLINVAR_VCF}" \
        -driver_gene_panel "${HMF_DRIVER_PANEL}" \
        -ensembl_data_dir "${HMF_ENSEMBL_DIR}" \
        -gnomad_freq_file "${HMF_WES_GNOMAD}" \
        -mappability_bed "${HMF_MAPPABILITY}" \
        -pon_artefact_file "${HMF_PANEL_ARTIFACT}" \
        -pon_file "${HMF_WES_PANEL_OF_NORMALS}" \
        -blacklist_bed "${HMF_GERMLINE_BLACKLIST_BED}" \
        -blacklist_vcf "${HMF_GERMLINE_BLACKLIST_VCF}" \
        -only_canonical \
        -threads "${Threads}" 2>&1 | tee -a "${RunLog}"
    #
#---------------------------------------------------------------------------------------------------

###---| AMBER |-------------------------------------------------------------------------------------

    if [ "${RunMode}" = "Tonly" ]; then
        OUTPUT_AMBER="${HMFTOOLS_RES_BASE}/amber_Tonly"
    else    
        OUTPUT_AMBER="${HMFTOOLS_RES_BASE}/amber"
    fi   
    mkdir -p "${OUTPUT_AMBER}"

    AMBER_INPUT_ARGS=("-tumor ${TumorSeqID}")
    AMBER_INPUT_ARGS+=("-tumor_bam ${TUMOR_REDUX_BAM}")
    if [[ "${RunMode}" == "NT" ]]; then
        AMBER_INPUT_ARGS+=("-reference ${NormalSeqID}")
        AMBER_INPUT_ARGS+=("-reference_bam ${NORMAL_REDUX_BAM}")
    else
        AMBER_INPUT_ARGS+=("-tumor_min_depth 2")
    fi

    ${HMF_RUN_AMBER} \
        ${AMBER_INPUT_ARGS[@]} \
        -output_dir "${OUTPUT_AMBER}" \
        -ref_genome_version "${HMF_GENOME}" \
        -ref_genome "${HMF_GENOME_FASTA}" \
        -loci "${HMF_WES_AMBER_LOCI}" \
        -target_regions_bed "${HMF_WES_TARGET_BED}" \
        -threads "${Threads}" 2>&1 | tee -a "${RunLog}"

#---------------------------------------------------------------------------------------------------

###---| COBALT |------------------------------------------------------------------------------------

    if [ "${RunMode}" = "Tonly" ]; then
        OUTPUT_COBALT="${HMFTOOLS_RES_BASE}/cobalt_Tonly"
    else    
        OUTPUT_COBALT="${HMFTOOLS_RES_BASE}/cobalt"
    fi   
    mkdir -p "${OUTPUT_COBALT}"

    COBLAT_INPUT_ARGS=("-tumor ${TumorSeqID}")
    COBLAT_INPUT_ARGS+=("-tumor_bam ${TUMOR_REDUX_BAM}")
    if [[ "${RunMode}" == "NT" ]]; then
        COBLAT_INPUT_ARGS+=("-reference ${NormalSeqID}")
        COBLAT_INPUT_ARGS+=("-reference_bam ${NORMAL_REDUX_BAM}")
    else
        COBLAT_INPUT_ARGS+=("-tumor_only_diploid_bed ${HMF_DIPLOID_REGIONS}")
        COBLAT_INPUT_ARGS+=("-pcf_gamma 50")
    fi

    ${HMF_RUN_COBALT} \
        ${COBLAT_INPUT_ARGS[@]} \
        -output_dir "${OUTPUT_COBALT}" \
        -gc_profile "${HMF_GC_PROFILE}" \
        -ref_genome_version "${HMF_GENOME}" \
        -ref_genome "${HMF_GENOME_FASTA}" \
        -threads "${Threads}" 2>&1 | tee -a "${RunLog}"

#---------------------------------------------------------------------------------------------------

###---| PURPLE |------------------------------------------------------------------------------------

    if [ "${RunMode}" = "Tonly" ]; then
        OUTPUT_PURPLE="${HMFTOOLS_RES_BASE}/purple_Tonly"
    else    
        OUTPUT_PURPLE="${HMFTOOLS_RES_BASE}/purple"
    fi   
    mkdir -p "${OUTPUT_COBALT}"

    PURPLE_INPUT_ARGS=("-tumor ${TumorSeqID}")
    if [[ "${RunMode}" == "NT" ]]; then
        PURPLE_INPUT_ARGS+=("-reference ${NormalSeqID}")
    fi

    ${HMF_RUN_PURPLE} \
        ${PURPLE_INPUT_ARGS[@]} \
        -amber "${OUTPUT_AMBER}" \
        -cobalt "${OUTPUT_COBALT}" \
        -somatic_vcf "${OUTPUT_PAVE}/${TumorSeqID}.pave.vcf.gz" \
        -germline_vcf "${OUTPUT_PAVE}/${TumorSeqID}.pave.germline.vcf.gz" \
        -ref_genome_version "${HMF_GENOME}" \
        -ref_genome "${HMF_GENOME_FASTA}" \
        -gc_profile "${HMF_GC_PROFILE}" \
        -ensembl_data_dir "${HMF_ENSEMBL_DIR}" \
        -driver_gene_panel "${HMF_DRIVER_PANEL}" \
        -somatic_hotspots "${HMF_SOMATIC_HOTSPOT_VCF}" \
        -germline_hotspots "${HMF_GERMLINE_HOTSPOT_VCF}" \
        -germline_del_freq_file "${HMF_GERMLINE_DEL_FREQ}" \
        -circos "/storage/apps/circos-0.69-10/bin/circos" \
        -target_regions_bed "${HMF_WES_TARGET_BED}" \
        -output_dir "${OUTPUT_PURPLE}" \
        -threads 1 2>&1 | tee -a "${RunLog}" 

#---------------------------------------------------------------------------------------------------


