
#! /bin/bash

#---| Input Options |------------------------------------------------------------------------------#
    usage() {
        echo "
        Usage: WGS_M05-1_RUN-HMFTOOLS-ANALYSIS.sh [ ARGUMENTS ]...

        --BASE_DIR                  [ REQUIRED ] analysis base folder by NGS data type
        
        --BATCH_ID                  [ REQUIRED ] dataset batch-id
        
        --TUMOR_SEQ_ID              [ REQURIED ] tumor sample-id

        --NORMAL_SEQ_ID             [ REQURIED ] normal sample-id
        
        --THREADS                   parallel processing cpu cores 
                                    default = 16

        --GENOME_ASSEMBLY           reference genome aseembly version
                                    'hg19' or 'hg38' enable
                                    default = 'hg38' 

        --WITH_RNASEQ               analysis with RNAseq data
                                    default = false

        --RNASEQ_HMFTOOLS_DATA_DIR  RNAseq hmftools analysis data folder
                                    required when 'WITH_RNASEQ' = true

        --RUN_BWA_MEM2              analysis start from bwa-meme2 mapping
                                    if 'false', bwa-mem2 aligned BAM file should be located in BAM-file folder
                                    default = false

        --PRIMARY_TUMOR_DOID        [ REQUIRED ]primary tumor tissue doid 

        -h | --help         Print this usage 
        "
    }
#--------------------------------------------------------------------------------------------------#

#---| OPTION PARAMETERS PARSING --------------------------------------------------------------------
    ARGS=$(getopt -a -o d:s:h: --long BASE_DIR:,BATCH_ID:,SEQ_ID:,TUMOR_SEQ_ID:,NORMAL_SEQ_ID:,THREADS:,GENOME_ASSEMBLY:,WITH_RNASEQ:,RNASEQ_HMFTOOLS_DATA_DIR:,RUN_BWA_MEM2:,PRIMARY_TUMOR_DOID:,help -- "$@" )
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
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  BaseDir=$2 ; shift 2 ;;
            esac ;;
        -s | --BATCH_ID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  BatchID=$2 ; shift 2 ;;
            esac  ;;
        --TUMOR_SEQ_ID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  TumorSeqID=$2 ; shift 2 ;;
            esac ;;
        --NORMAL_SEQ_ID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  NormalSeqID=$2 ; shift 2 ;;
            esac ;;
        --THREADS )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 24" ; shift ;;
                *)  Threads=$2 ; shift 2 ;;
            esac ;;
        --GENOME_ASSEMBLY )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = hg38" ; shift ;;
                *)  GenomeAssembly=$2 ; shift 2 ;;
            esac ;;
        --WITH_RNASEQ )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = false" ; shift ;;
                *)  WithRnaseqData=$2 ; shift 2 ;;
            esac ;;
        --RNASEQ_HMFTOOLS_DATA_DIR )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = null" ; shift ;;
                *)  RnaseqHmftoolsDataDir=$2 ; shift 2 ;;
            esac ;;
        --RUN_BWA_MEM2 )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = false" ; shift ;;
                *)  RunBwaMem2Mapping=$2 ; shift 2 ;;
            esac ;;
        --PRIMARY_TUMOR_DOID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  PrimaryTumorDoid=$2 ; shift 2 ;;
            esac ;;
        -h | --help )
        usage >&2 ; exit 2 ;;   
        --) shift ; break ;;
        *)  usage >&2 ; exit 2 ;;
    esac
    done
#---------------------------------------------------------------------------------------------------

###---| Check Params |------------------------------------------------------------------------------
    if [ -z "${Threads}" ]; then Threads=24; fi
    if [ -z "${GenomeAssembly}" ]; then GenomeAssembly="hg38"; fi
    if [ -z "${GenomeAssembly}" ]; then GenomeAssembly="hg38"; fi
    if [ -z "${WithRnaseqData}" ]; then WithRnaseqData="false"; fi
    if [ -z "${RunBwaMem2Mapping}" ]; then RunBwaMem2Mapping="false"; fi
    if [ "${WithRnaseqData}" == "true" ]; then
        if [ -z ${RnaseqHmftoolsDataDir} ]; then 
            echo "*** NO RNAseq-HMFTOOLS-DATA-FOLDER. REQUIRED --WITH_RNASEQ=true " 2>&1 | tee -a "${RunLog}"
            exit 2
        fi
    fi
###-------------------------------------------------------------------------------------------------

#---/ Manual Params /-------------------------------------------------------------------------------
    BaseDir="/data/wgs/RUO_Service"
    BatchID="WGS-26-02"
    TumorSeqID="WGS-26-02_03"
    NormalSeqID="WGS-26-02_01"
    Threads=24
    GenomeAssembly="hg19"
    PrimaryTumorDoid="1612"
    RunBwaMem2Mapping="false"
    WithRnaseqData="false"
    # RnaseqHmftoolsDataDir=
#---------------------------------------------------------------------------------------------------

###---| RUN HMFTOOLS |------------------------------------------------------------------------------

    RUN_SIF="singularity exec -B /storage,/data /storage/images"
    RUN_BWA_MEM2="${RUN_SIF}/bwa-mem2-v2.2.1.sif bwa-mem2"
    HMF_RUN_REDUX="${RUN_SIF}/hmftools-redux-v1.2.3.sif redux"
    HMF_RUN_BAMTOOLS="${RUN_SIF}/hmftools-bamtools-v1.5.sif bamtools"
    HMF_RUN_SAGE="${RUN_SIF}/hmftools-sage-v4.2.sif sage"
    HMF_RUN_PAVE="${RUN_SIF}/hmftools-pave-v1.8.2.sif pave"
    HMF_RUN_PAVE_JAR="${RUN_SIF}/hmftools-pave-v1.8.2.sif java -Xmx32G -jar /usr/share/java/pave_v1.8.2.jar"
    HMF_RUN_AMBER="${RUN_SIF}/hmftools-amber-v4.2.sif amber"
    HMF_RUN_COBALT="${RUN_SIF}/hmftools-cobalt-v2.2.sif cobalt"
    HMF_RUN_ESVEE="${RUN_SIF}/hmftools-esvee-v1.2.sif esvee"
    HMF_RUN_PURPLE="${RUN_SIF}/hmftools-purple-v4.3.sif java -Xmx48G -jar /usr/share/java/purple_v4.3.jar"   
    HMF_RUN_LINX="${RUN_SIF}/hmftools-linx-v2.2.sif linx"
    HMF_RUN_LILAC="${RUN_SIF}/hmftools-lilac-v1.7.3.sif lilac"
    HMF_RUN_CHORD="${RUN_SIF}/hmftools-chord-v2.1.2.sif chord"
    HMF_RUN_CUPPA="${RUN_SIF}/hmftools-cuppa-v2.4.sif cuppa"
    HMF_RUN_NEO="${RUN_SIF}/hmftools-neo-v1.2.1.sif neo"
    HMF_RUN_TEAL="${RUN_SIF}/hmftools-teal-v1.3.6.sif java -Xmx16G -cp /usr/share/java/teal_v1.3.6.jar com.hartwig.hmftools.teal.TealPipelineApp"
    HMF_RUN_ORANGE="${RUN_SIF}/hmftools-orange-v4.1.3.sif orange"
    HMF_RUN_SIGS="java -jar /storage/apps/hmftools/sigs_v1.2.jar"
    HMF_RUN_PEACH="java -jar /storage/apps/hmftools/peach_v2.0.0.jar"
    HMF_RUN_VIRUS_INTERPRETATER="java -jar /storage/apps/hmftools/virus-interpreter_v1.5.0.jar"
    HMF_RUN_WISP="java -jar /storage/apps/hmftools/wisp_v1.2.jar"
    HMF_RUN_VIRUS_BREAKEND="/storage/apps/gridss-v2.13.2/virusbreakend"
    HMF_RUN_GRIDSS="/storage/apps/gridss-v2.13.2/gridss"
    HMF_GRIDSS_JAR="/storage/apps/gridss-v2.13.2/gridss-2.13.2-gridss-jar-with-dependencies.jar"

###-------------------------------------------------------------------------------------------------

###---| HMFTOOLS RESOURCES |------------------------------------------------------------------------
    if [ "${GenomeAssembly}" = "hg19" ]; then
        HMF_GENOME_FASTA="/storage/references_and_index/hg19/hmftools/Homo_sapiens.GRCh37.GATK.illumina.fasta"
        HMF_BWAMEM2_INDEX="/storage/references_and_index/hg19/hmftools/bwa-mem2_index-2.2.1/Homo_sapiens.GRCh37.GATK.illumina.fasta"
        HMF_GENOME="37"
    else
        HMF_GENOME_FASTA="/storage/references_and_index/hg38/hmftools/GRCh38_masked_exclusions_alts_hlas.fasta"
        HMF_BWAMEM2_INDEX="/storage/references_and_index/hg38/hmftools/bwa-mem2_index-2.2.1/GRCh38_masked_exclusions_alts_hlas.fasta"
        HMF_GENOME="38"
    fi

    HMF_RESOURCE="/storage/references_and_index/${GenomeAssembly}/hmftools/resources_links"
    
    HMF_UNMAP_REGION="${HMF_RESOURCE}/UNMAP-REGIONS.tsv"
    HMF_DRIVER_PANEL="${HMF_RESOURCE}/DRIVER-GENE-PANEL.tsv"
    HMF_ENSEMBL_DIR="${HMF_RESOURCE}/ENSEMBL-DATA/"
    HMF_MSI_JITTER_SITES="${HMF_RESOURCE}/MSI-JITTER-SITES.tsv.gz"
    HMF_SOMATIC_HOTSPOT_VCF="${HMF_RESOURCE}/KNOWN-HOTSPOT-SOMATIC.vcf.gz"
    HMF_HIGH_CONFIDENCE_BED="${HMF_RESOURCE}/HIGH-CONFIDENCE.bed.gz"
    HMF_GNOMAD_DIR="${HMF_RESOURCE}/GNOMAD-DATA/"
    HMF_GNOMAD_VARS="${HMF_RESOURCE}/GNOMAD-VARIANTS.csv.gz"
    HMF_PON="${HMF_RESOURCE}/PANEL-OF-NORMALS.tsv.gz"
    HMF_ACTIONAL_CODING_PANEL="${HMF_RESOURCE}/ACTIONABLE-CODING-PANEL.bed.gz"
    HMF_COVERAGE_CODINGPANEL="${HMF_RESOURCE}/COVERAGE-CODING-PANEL.bed.gz"
    HMF_MAPPABILITY="${HMF_RESOURCE}/MAPPABILITY.bed.gz"
    HMF_CLINVAR_VCF="${HMF_RESOURCE}/CLINVAR.vcf.gz"
    HMF_PANEL_ARTIFACT="${HMF_RESOURCE}/PANEL-ARTIFACTS.tsv"
    HMF_GERMLINE_BLACKLIST_BED="${HMF_RESOURCE}/GERMLINE-BLACKLIST-REGION.bed"
    HMF_GERMLINE_BLACKLIST_VCF="${HMF_RESOURCE}/GERMLINE-BLACKLIST-VARS.vcf.gz"
    HMF_GERMLINE_HOTSPOT_VCF="${HMF_RESOURCE}/KNOWN-HOTSPOT-GERMLINE.vcf.gz"
    HMF_PON_SV="${HMF_RESOURCE}/SV-PON.bedpe.gz"
    HMF_PON_SGL="${HMF_RESOURCE}/SGL-PON.bed.gz"
    HMF_REPEAT_MASK="${HMF_RESOURCE}/REPEAT-MASK.fa.gz"
    HMF_SV_BLACKLIST_BED="${HMF_RESOURCE}/SV-BLACKLIST.bed"
    HMF_KNOWN_FUSION_BED="${HMF_RESOURCE}/KNOWN-FUSIONS.bedpe"
    HMF_KNOWN_FUSION_CSV="${HMF_RESOURCE}/KNOWN-FUSION-DATA.csv"
    HMF_GC_PROFILE="${HMF_RESOURCE}/GC-PROFILE.cnp"
    HMF_AMBER_LOCI="${HMF_RESOURCE}/AMBER-GERMLINE-SITES.tsv.gz"
    HMF_DIPLOID_BED="${HMF_RESOURCE}/DIPLOID-REGIONS.bed.gz"
    HMF_GERMLINE_DEL_FREQ="${HMF_RESOURCE}/COHORT-GERMLINE-DELETION-FREQ.csv"
    HMF_LILAC_RESOURCE_DIR="${HMF_RESOURCE}/LILAC-RESOURCES"
    HMF_SIGS_SIGNATURE="${HMF_RESOURCE}/COSMIC-SNV-SIGNATURES.csv"
    HMF_SIGATURE_ETIOLOGY="${HMF_RESOURCE}/SIGNATURE_ETIOLOGY.tsv"
    HMF_PEACH_HAPLOTYPE="${HMF_RESOURCE}/HAPLOTYPES.tsv"
    HMF_HAPLOTYPE_FUNCTION="${HMF_RESOURCE}/HAPLOTYPE_FUNCTION.tsv"
    HMF_PEACH_DRUGS="${HMF_RESOURCE}/PEACH_DRUGS.tsv"
    HMF_CUPPA_PICKLE="${HMF_RESOURCE}/CUPPA-CLASSFIER.pickle.gz"
    HMF_CUPPA_SJ_SITES="${HMF_RESOURCE}/CUPPA-ALT-SJ-LOCI.tsv.gz"
    HMF_VIRUSBREAKEND_DB="${HMF_RESOURCE}/VIRUSBREAKEND_DB/"
    HMF_VIRUS_TAX_DB="${HMF_RESOURCE}/TAXONOMY_DB.tsv"
    HMF_VIRUS_BLACKLIST_DB="${HMF_RESOURCE}/VIRUS_BLACKLIST_DB.tsv"
    HMF_VIRUS_REPORT_DB="${HMF_RESOURCE}/VIRUS_REPORTING_DB.tsv"
    HMF_COHORT_MAPPING="${HMF_RESOURCE}/COHORT_MAPPING.tsv"
    HMF_COHORT_PERCENTILES="${HMF_RESOURCE}/COHORT_PERCENTILES.tsv"
    HMF_DOID_JSON="${HMF_RESOURCE}/DOID.json"


    if [[ "${HMF_GENOME}" = "37" ]]; then
        GNOMAD_RESOURCE_ARGS=("-gnomad_freq_file ${HMF_GNOMAD_VARS}")
    else
        GNOMAD_RESOURCE_ARGS=("-gnomad_freq_dir ${HMF_GNOMAD_DIR}")
    fi

###-------------------------------------------------------------------------------------------------

###---| VEP RESOURCES |-----------------------------------------------------------------------------
    # Run VEP
    RunVep="singularity exec -B /storage,/data /storage/images/vep-110.1.sif /opt/vep/src/ensembl-vep/vep"
    RunVcf2Maf="perl /storage/apps/vcf2maf-1.6.22/vcf2maf.pl"

    ## VEP Resources
    VepCacheDir=/storage/references_and_index/${GenomeAssembly}/vep
    VepPluginDir=${VepCacheDir}/Plugins/110
    VepBam=${VepCacheDir}/bam/${GenomeAssembly}.vep.bam
    VepTopLevelFasta=${VepCacheDir}/toplevel_fasta/${GenomeAssembly}.toplevel.fa.gz
    VepData_Clinvar=/storage/references_and_index/${GenomeAssembly}/annotation/clinvar/clinvar.vcf.gz
    VepData_AlphaMissense=/storage/references_and_index/${GenomeAssembly}/annotation/AlphaMissense/alpha-missense.tsv.gz
    VepData_Cosmic=/storage/references_and_index/${GenomeAssembly}/annotation/cosmic/Cosmic_GenomeScreensMutant_v103.vcf.gz
###-------------------------------------------------------------------------------------------------

###---| INPUT/OUTPUT DATA & FOLDERS |---------------------------------------------------------------

    # tumor/normal data folders
    TumorDataDir="${BaseDir}/${BatchID}/${TumorSeqID}"
    NormalDataDir="${BaseDir}/${BatchID}/${NormalSeqID}"

    # hmftools output base folder
    HMF_RES_DIR="${BaseDir}/${BatchID}/${TumorSeqID}/hmftools"
    mkdir -p "${HMF_RES_DIR}/runlog"

    # run log 
    RunLog="${HMF_RES_DIR}/runlog/hmftools.run.$(date '+%Y-%m-%d').log"
    echo " " > "${RunLog}"
    echo "# HMFTOOLS ANALYSIS RUN LOG" >> "${RunLog}"
    echo " " >> "${RunLog}"
    echo " TUMOR            ${TumorSeqID}" >> "${RunLog}"
    echo " MATCHED-NORMAL   ${NormalSeqID}" >> "${RunLog}"
    echo " " >> "${RunLog}"

###-------------------------------------------------------------------------------------------------

###---| BWA-MEM2 ALIGN |----------------------------------------------------------------------------
    if [[ "${RunBwaMem2Mapping}" = "true" ]]; then

        echo "###-----/ BWA-MEM2 /----------------------------------------------------------------------###" >> "${RunLog}"

        case "${Threads}" in
            24) SortThreads=6 ;;
            32) SortThreads=8 ;;
            *) SortThreads=4 ;; 
        esac
        BwaThreads=$(( Threads - SortThreads ))

        # tumor
        echo "--- TUMOR BWA-MEM2 BAM ---" >> "${RunLog}"
        mkdir -p "${BaseDir}/${BatchID}/${TumorSeqID}/bam"
        "${RUN_BWA_MEM2}" mem -M -Y -q -K 100000000 -t "${BwaThreads}" \
            -R "@RG\tID:HMFTOOLS\tSM:${TumorSeqID}\tLB:WGS\tPL:ILLUMINA\tCN:GCX" \
            "${HMF_BWAMEM2_INDEX}" \
            "${BaseDir}/${BatchID}/${TumorSeqID}/fastq/${TumorSeqID}.trimmed_R1.fastq.gz" \
            "${BaseDir}/${BatchID}/${TumorSeqID}/fastq/${TumorSeqID}.trimmed_R2.fastq.gz" \
            | samtools sort -@ "${SortThreads}" -o "${BaseDir}/${BatchID}/${TumorSeqID}/bam/${TumorSeqID}.bwa.mem2.bam" - 
        samtools index -@ "${Threads}" "${BaseDir}/${BatchID}/${TumorSeqID}/bam/${TumorSeqID}.bwa.mem2.bam" 
        echo " " >> "${RunLog}"

        # normal 
        echo "--- NORMAL BWA-MEM2 BAM ---" >> "${RunLog}"
        mkdir -p "${BaseDir}/${BatchID}/${NormalSeqID}/bam"
        "${RUN_BWA_MEM2}" mem -M -Y -q -K 100000000 -t "${BwaThreads}" \
            -R "@RG\tID:HMFTOOLS\tSM:${NormalSeqID}\tLB:WGS\tPL:ILLUMINA\tCN:GCX" \
            "${HMF_BWAMEM2_INDEX}" \
            "${BaseDir}/${BatchID}/${NormalSeqID}/fastq/${NormalSeqID}.trimmed_R1.fastq.gz" \
            "${BaseDir}/${BatchID}/${NormalSeqID}/fastq/${NormalSeqID}.trimmed_R2.fastq.gz" \
            | samtools sort -@ "${SortThreads}" -o "${BaseDir}/${BatchID}/${NormalSeqID}/bam/${NormalSeqID}.bwa.mem2.bam" - 
        
        samtools index -@ "${Threads}" "${BaseDir}/${BatchID}/${NormalSeqID}/bam/${NormalSeqID}.bwa.mem2.bam" 

        echo " " >> "${RunLog}"
    else
        TumorHmftoolsBam="${BaseDir}/${BatchID}/${TumorSeqID}/bam/${TumorSeqID}.bwa.mem2.bam" 
        NormalHmftoolsBam="${BaseDir}/${BatchID}/${NormalSeqID}/bam/${NormalSeqID}.bwa.mem2.bam"
        if [[ ! -f "${TumorHmftoolsBam}" || ! -f "${NormalHmftoolsBam}" ]]; then
            echo "NO HMFTOOLS ANALYSIS-READY BAM FILE FOUND." 2>&1 | tee -a "${RunLog}" 
            echo "Please Do Alignment using bwa-mem2 first OR Set '--RUN_BWA_MEM2' to 'true'" 2>&1 | tee -a "${RunLog}"
            exit 2
        fi
    fi
#---------------------------------------------------------------------------------------------------

###---| HMFTOOLS : REDUX |--------------------------------------------------------------------------

    echo "###-----/ REDUX /-------------------------------------------------------------------------###" >> "${RunLog}"
    
    REDUX_OUTPUT="${HMF_RES_DIR}/redux"
    mkdir -p "${REDUX_OUTPUT}"

    # tumor
    echo "--- TUMOR BAM ---" >> "${RunLog}"
    ${HMF_RUN_REDUX} \
        -sample "${TumorSeqID}" \
        -input_bam "${BaseDir}/${BatchID}/${TumorSeqID}/bam/${TumorSeqID}.bwa.mem2.bam" \
        -output_bam "${REDUX_OUTPUT}/${TumorSeqID}.redux.bam" \
        -output_dir "${REDUX_OUTPUT}" \
        -ref_genome_version "${HMF_GENOME}" \
        -ref_genome "${HMF_GENOME_FASTA}" \
        -unmap_regions "${HMF_UNMAP_REGION}" \
        -ref_genome_msi_file "${HMF_MSI_JITTER_SITES}" \
        -form_consensus \
        -write_stats \
        -bamtool "/opt/conda/bin/samtools" \
        -threads "${Threads}" 2>&1 | tee -a "${RunLog}"

    echo " " >> "${RunLog}"

    # normal
    echo "--- MATCHED-NORMAL BAM ---" >> "${RunLog}"
    REDUX_NORMAL_OUTPUT="${BaseDir}/${BatchID}/${NormalSeqID}/hmftools/redux"
    mkdir -p ${REDUX_NORMAL_OUTPUT}

    if [[ ! -f "${REDUX_NORMAL_OUTPUT}/${NormalSeqID}.redux.bam" ]]; then

        ${HMF_RUN_REDUX} \
            -sample "${NormalSeqID}" \
            -input_bam "${BaseDir}/${BatchID}/${NormalSeqID}/bam/${NormalSeqID}.bwa.mem2.bam" \
            -output_bam "${REDUX_NORMAL_OUTPUT}/${NormalSeqID}.redux.bam" \
            -output_dir "${REDUX_NORMAL_OUTPUT}" \
            -ref_genome_version "${HMF_GENOME}" \
            -ref_genome "${HMF_GENOME_FASTA}" \
            -unmap_regions "${HMF_UNMAP_REGION}" \
            -ref_genome_msi_file "${HMF_MSI_JITTER_SITES}" \
            -form_consensus \
            -write_stats \
            -bamtool "/opt/conda/bin/samtools" \
            -threads "${Threads}" 2>&1 | tee -a "${RunLog}"

        cp ${REDUX_NORMAL_OUTPUT}/${NormalSeqID}.redux.bam ${REDUX_OUTPUT}/
        cp ${REDUX_NORMAL_OUTPUT}/${NormalSeqID}.redux.bam.bai ${REDUX_OUTPUT}/
        cp ${REDUX_NORMAL_OUTPUT}/${NormalSeqID}.jitter_params.tsv ${REDUX_OUTPUT}/
        cp ${REDUX_NORMAL_OUTPUT}/${NormalSeqID}.ms_table.tsv.gz ${REDUX_OUTPUT}/
        cp ${REDUX_NORMAL_OUTPUT}/${NormalSeqID}.redux.duplicate_freq.tsv ${REDUX_OUTPUT}/

    fi
    
    echo " " >> "${RunLog}"

#---------------------------------------------------------------------------------------------------

###---| HMFTOOLS ANALYSIS INPUT BAM |---------------------------------------------------------------

    TUMOR_REDUX_BAM="${REDUX_OUTPUT}/${TumorSeqID}.redux.bam"
    NORMAL_REDUX_BAM="${REDUX_OUTPUT}/${NormalSeqID}.redux.bam"

###-------------------------------------------------------------------------------------------------

###---| BAM-MATRICS |-------------------------------------------------------------------------------

    echo "###-----/ BamTools - BAM Metrics /--------------------------------------------------------###" >> "${RunLog}"
    TUMOR_BAMTOOLS_OUTPUT="${HMF_RES_DIR}/tumor_bam_metric"
    NORMAL_BAMTOOLS_OUTPUT="${HMF_RES_DIR}/normal_bam_metric"
    mkdir -p "${TUMOR_BAMTOOLS_OUTPUT}"
    mkdir -p "${NORMAL_BAMTOOLS_OUTPUT}"

    # tumor bam
    ${HMF_RUN_BAMTOOLS} \
        -sample "${TumorSeqID}" \
        -bam_file "${TUMOR_REDUX_BAM}" \
        -output_dir "${TUMOR_BAMTOOLS_OUTPUT}" \
        -driver_gene_panel "${HMF_DRIVER_PANEL}" \
        -ensembl_data_dir "${HMF_ENSEMBL_DIR}" \
        -ref_genome_version "${HMF_GENOME}" \
        -ref_genome "${HMF_GENOME_FASTA}" \
        -threads "${Threads}" 2>&1 | tee -a "${RunLog}"

    echo " " >> "${RunLog}"

    # normal bam
    ${HMF_RUN_BAMTOOLS} \
        -sample "${NormalSeqID}" \
        -bam_file "${NORMAL_REDUX_BAM}" \
        -output_dir "${NORMAL_BAMTOOLS_OUTPUT}" \
        -driver_gene_panel "${HMF_DRIVER_PANEL}" \
        -ensembl_data_dir "${HMF_ENSEMBL_DIR}" \
        -ref_genome_version "${HMF_GENOME}" \
        -ref_genome "${HMF_GENOME_FASTA}" \
        -threads "${Threads}" 2>&1 | tee -a "${RunLog}"
 
    echo " " >> "${RunLog}"
###-------------------------------------------------------------------------------------------------

###---| GATK : WGS-MATRICS |------------------------------------------------------------------------

    echo "###-----/ GATK - WGS Metrics /------------------------------------------------------------###" >> "${RunLog}"

    # GATK Run
    GATK_VERSION="4.6.2.0"
    RUN_GATK4="singularity exec -B /storage,/data /storage/images/gatk-${GATK_VERSION}.sif gatk"

    ${RUN_GATK4} --java-options "-Xmx16g" CollectWgsMetrics \
        --INPUT "${REDUX_OUTPUT}/${TumorSeqID}.redux.bam" \
        --OUTPUT "${REDUX_OUTPUT}/${TumorSeqID}.redux.wgs.metrics.txt" \
        --REFERENCE_SEQUENCE "${HMF_GENOME_FASTA}" \
        --VALIDATION_STRINGENCY LENIENT 2>&1 | tee -a "${RunLog}"
    

    if [[ ! -f "${REDUX_NORMAL_OUTPUT}/${NormalSeqID}.redux.wgs.metrics.txt" ]]; then

        ${RUN_GATK4} --java-options "-Xmx16g" CollectWgsMetrics \
            --INPUT "${REDUX_NORMAL_OUTPUT}/${NormalSeqID}.redux.bam" \
            --OUTPUT "${REDUX_NORMAL_OUTPUT}/${NormalSeqID}.redux.wgs.metrics.txt" \
            --REFERENCE_SEQUENCE "${HMF_GENOME_FASTA}" \
            --VALIDATION_STRINGENCY LENIENT 2>&1 | tee -a "${RunLog}"

        cp ${REDUX_NORMAL_OUTPUT}/${NormalSeqID}.redux.wgs.metrics.txt ${REDUX_OUTPUT}/

    fi

    echo " " >> "${RunLog}"

###-------------------------------------------------------------------------------------------------

###---| HMFTOOLS : SAGE |---------------------------------------------------------------------------

    echo "###-----/ SAGE /--------------------------------------------------------------------------###" >> "${RunLog}"

    # somatic
    echo "--- SAGE SOMATIC VARIANTS ---" >> "${RunLog}"

    SAGE_OUTPUT="${HMF_RES_DIR}/sage"
    mkdir -p "${SAGE_OUTPUT}"
   
    ${HMF_RUN_SAGE} \
        -tumor "${TumorSeqID}" \
        -tumor_bam "${TUMOR_REDUX_BAM}" \
        -reference "${NormalSeqID}" \
        -reference_bam "${NORMAL_REDUX_BAM}" \
        -output_vcf "${SAGE_OUTPUT}/${TumorSeqID}.sage.vcf.gz" \
        -jitter_param_dir "${REDUX_OUTPUT}" \
        -ref_genome_version "${HMF_GENOME}" \
        -ref_genome "${HMF_GENOME_FASTA}" \
        -hotspots "${HMF_SOMATIC_HOTSPOT_VCF}" \
        -high_confidence_bed "${HMF_HIGH_CONFIDENCE_BED}" \
        -ensembl_data_dir "${HMF_ENSEMBL_DIR}" \
        ${GNOMAD_RESOURCE_ARGS[@]} \
        -driver_gene_panel "${HMF_DRIVER_PANEL}" \
        -pon_file "${HMF_PON}" \
        -bqr_write_plot \
        -threads "${Threads}" 2>&1 | tee -a "${RunLog}"

    echo " " >> "${RunLog}"

    # germline
    echo "--- SAGE GERMLINE VARIANTS ---" >> "${RunLog}"

    SAGE_GERMLINE_OUTPUT="${HMF_RES_DIR}/sage_germline"
    mkdir -p "${SAGE_GERMLINE_OUTPUT}"
   
    # actional-coding-panel only
    ${HMF_RUN_SAGE} \
        -germline \
        -tumor "${NormalSeqID}" \
        -tumor_bam "${NORMAL_REDUX_BAM}" \
        -reference "${TumorSeqID}"\
        -reference_bam "${TUMOR_REDUX_BAM}" \
        -output_vcf "${SAGE_GERMLINE_OUTPUT}/${TumorSeqID}.sage.germline.vcf.gz" \
        -jitter_param_dir "${REDUX_OUTPUT}" \
        -ref_genome_version "${HMF_GENOME}" \
        -ref_genome "${HMF_GENOME_FASTA}" \
        -ref_sample_count 0 \
        -ensembl_data_dir "${HMF_ENSEMBL_DIR}" \
        -hotspots "${HMF_GERMLINE_HOTSPOT_VCF}" \
        ${GNOMAD_RESOURCE_ARGS[@]} \
        -pon_file "${HMF_PON}" \
        -panel_only \
        -panel_bed "${HMF_ACTIONAL_CODING_PANEL}" \
        -bqr_write_plot \
        -threads "${Threads}" 2>&1 | tee -a "${RunLog}"
    #
    # all genome sites
    ${HMF_RUN_SAGE} \
        -germline \
        -tumor "${NormalSeqID}" \
        -tumor_bam "${NORMAL_REDUX_BAM}" \
        -reference "${TumorSeqID}"\
        -reference_bam "${TUMOR_REDUX_BAM}" \
        -output_vcf "${SAGE_GERMLINE_OUTPUT}/${TumorSeqID}.sage.germline.all.sites.vcf.gz" \
        -jitter_param_dir "${REDUX_OUTPUT}" \
        -ref_genome_version "${HMF_GENOME}" \
        -ref_genome "${HMF_GENOME_FASTA}" \
        -ref_sample_count 0 \
        -ensembl_data_dir "${HMF_ENSEMBL_DIR}" \
        -hotspots "${HMF_GERMLINE_HOTSPOT_VCF}" \
        ${GNOMAD_RESOURCE_ARGS[@]} \
        -pon_file "${HMF_PON}" \
        -bqr_write_plot \
        -threads "${Threads}" 2>&1 | tee -a "${RunLog}"
    #
    echo " " >> "${RunLog}"

###-------------------------------------------------------------------------------------------------

###---| HMFTOOLS : PAVE |---------------------------------------------------------------------------
    
    echo "###-----/ PAVE /--------------------------------------------------------------------------###" >> "${RunLog}"

    PAVE_OUTPUT="${HMF_RES_DIR}/pave"
    mkdir -p "${PAVE_OUTPUT}"

    # somatic
    echo "--- PAVE SOMATIC VARIANTS ---" >> "${RunLog}"

    SAGE_VCF="${SAGE_OUTPUT}/${TumorSeqID}.sage.vcf.gz"
    PAVE_VCF="${PAVE_OUTPUT}/${TumorSeqID}.pave.vcf.gz"

    ${HMF_RUN_PAVE} \
        -sample "${TumorSeqID}" \
        -input_vcf "${SAGE_OUTPUT}/${TumorSeqID}.sage.vcf.gz" \
        -output_dir "${PAVE_OUTPUT}" \
        -output_vcf "${PAVE_OUTPUT}/${TumorSeqID}.pave.vcf.gz" \
        -ref_genome_version "${HMF_GENOME}" \
        -ref_genome "${HMF_GENOME_FASTA}" \
        -clinvar_vcf "${HMF_CLINVAR_VCF}" \
        -driver_gene_panel "${HMF_DRIVER_PANEL}" \
        -ensembl_data_dir "${HMF_ENSEMBL_DIR}" \
        ${GNOMAD_RESOURCE_ARGS[@]} \
        -mappability_bed "${HMF_MAPPABILITY}" \
        -pon_artefact_file "${HMF_PANEL_ARTIFACT}" \
        -pon_file "${HMF_PON}" \
        -blacklist_bed "${HMF_GERMLINE_BLACKLIST_BED}" \
        -blacklist_vcf "${HMF_GERMLINE_BLACKLIST_VCF}" \
        -only_canonical \
        -write_detailed \
        -write_transcript_data \
        -threads "${Threads}" 2>&1 | tee -a "${RunLog}"
    
    echo " " >> "${RunLog}"

    # germline
    echo "--- PAVE GERMLINE VARIANTS ---" >> "${RunLog}"

    SAGE_GERMLINE_VCF="${SAGE_GERMLINE_OUTPUT}/${TumorSeqID}.sage.germline.vcf.gz"
    PAVE_GERMLINE_VCF="${PAVE_OUTPUT}/${TumorSeqID}.pave.germline.vcf.gz"

    ${HMF_RUN_PAVE_JAR} \
        -sample "${TumorSeqID}" \
        -input_vcf "${SAGE_GERMLINE_OUTPUT}/${TumorSeqID}.sage.germline.vcf.gz" \
        -output_dir "${PAVE_OUTPUT}" \
        -output_vcf "${PAVE_OUTPUT}/${TumorSeqID}.pave.germline.vcf.gz" \
        -ref_genome_version "${HMF_GENOME}" \
        -ref_genome "${HMF_GENOME_FASTA}" \
        -clinvar_vcf "${HMF_CLINVAR_VCF}" \
        -driver_gene_panel "${HMF_DRIVER_PANEL}" \
        -ensembl_data_dir "${HMF_ENSEMBL_DIR}" \
        ${GNOMAD_RESOURCE_ARGS[@]} \
        -mappability_bed "${HMF_MAPPABILITY}" \
        -pon_artefact_file "${HMF_PANEL_ARTIFACT}" \
        -pon_file "${HMF_PON}" \
        -blacklist_bed "${HMF_GERMLINE_BLACKLIST_BED}" \
        -blacklist_vcf "${HMF_GERMLINE_BLACKLIST_VCF}" \
        -write_pass_only \
        -threads 1 2>&1 | tee -a "${RunLog}"
    
    echo " " >> "${RunLog}"

###-------------------------------------------------------------------------------------------------

###---| HMFTOOLS : AMBER |--------------------------------------------------------------------------

    echo "###-----/ AMBER /-------------------------------------------------------------------------###" >> "${RunLog}"

    AMBER_OUTPUT="${HMF_RES_DIR}/amber"
    mkdir -p "${AMBER_OUTPUT}"

    ${HMF_RUN_AMBER} \
        -tumor "${TumorSeqID}" \
        -tumor_bam "${TUMOR_REDUX_BAM}" \
        -reference "${NormalSeqID}" \
        -reference_bam "${NORMAL_REDUX_BAM}" \
        -output_dir "${AMBER_OUTPUT}" \
        -ref_genome_version "${HMF_GENOME}" \
        -ref_genome "${HMF_GENOME_FASTA}" \
        -loci "${HMF_AMBER_LOCI}" \
        -threads "${Threads}" 2>&1 | tee -a "${RunLog}"

    echo " " >> "${RunLog}"

###-------------------------------------------------------------------------------------------------

###---| HMFTOOLS : COBALT |-------------------------------------------------------------------------

    echo "###-----/ COBALT /-------------------------------------------------------------------------###" >> "${RunLog}"

    COBALT_OUTPUT="${HMF_RES_DIR}/cobalt"
    mkdir -p "${COBALT_OUTPUT}"

    ${HMF_RUN_COBALT} \
        -tumor "${TumorSeqID}" \
        -tumor_bam "${TUMOR_REDUX_BAM}" \
        -reference "${NormalSeqID}" \
        -reference_bam "${NORMAL_REDUX_BAM}" \
        -output_dir "${COBALT_OUTPUT}" \
        -gc_profile "${HMF_GC_PROFILE}" \
        -ref_genome_version "${HMF_GENOME}" \
        -ref_genome "${HMF_GENOME_FASTA}" \
        -threads "${Threads}" 2>&1 | tee -a "${RunLog}"

    echo " " >> "${RunLog}"

###-------------------------------------------------------------------------------------------------

###---| HMFTOOLS : ESVEE |--------------------------------------------------------------------------

    echo "###-----/ ESVEE /-------------------------------------------------------------------------###" >> "${RunLog}"

    ESVEE_OUTPUT="${HMF_RES_DIR}/esvee"
    mkdir -p "${ESVEE_OUTPUT}"

    ${HMF_RUN_ESVEE} \
        -tumor "${TumorSeqID}" \
        -tumor_bam "${TUMOR_REDUX_BAM}" \
        -reference "${NormalSeqID}" \
        -reference_bam "${NORMAL_REDUX_BAM}" \
        -output_dir "${ESVEE_OUTPUT}" \
        -ref_genome_version "${HMF_GENOME}" \
        -ref_genome "${HMF_GENOME_FASTA}" \
        -known_hotspot_file "${HMF_KNOWN_FUSION_BED}" \
        -pon_sgl_file "${HMF_PON_SGL}" \
        -pon_sv_file "${HMF_PON_SV}" \
        -target_regions_bed "${HMF_COVERAGE_CODINGPANEL}" \
        -repeat_mask_file "${HMF_REPEAT_MASK}" \
        -blacklist_bed "${HMF_SV_BLACKLIST_BED}" \
        -unmap_regions "${HMF_UNMAP_REGION}" \
        -write_types "PREP_STANDARD;ASSEMBLY_STANDARD" \
        -bamtool "/opt/conda/bin/sambamba" \
        -write_breakend_tsv \
        -threads "${Threads}" 2>&1 | tee -a "${RunLog}"

    echo " " >> "${RunLog}"

###-------------------------------------------------------------------------------------------------

###---| HMFTOOLS : PURPLE |-------------------------------------------------------------------------

    echo "###-----/ PURPLE /------------------------------------------------------------------------###" >> "${RunLog}"

    PURPLE_OUTPUT="${HMF_RES_DIR}/purple"
    mkdir -p "${PURPLE_OUTPUT}"

    ${HMF_RUN_PURPLE} \
        -tumor "${TumorSeqID}" \
        -reference "${NormalSeqID}" \
        -amber "${AMBER_OUTPUT}" \
        -cobalt "${COBALT_OUTPUT}" \
        -somatic_vcf "${PAVE_OUTPUT}/${TumorSeqID}.pave.vcf.gz" \
        -germline_vcf "${PAVE_OUTPUT}/${TumorSeqID}.pave.germline.vcf.gz" \
        -somatic_sv_vcf "${ESVEE_OUTPUT}/${TumorSeqID}.esvee.somatic.vcf.gz"  \
        -germline_sv_vcf "${ESVEE_OUTPUT}/${TumorSeqID}.esvee.germline.vcf.gz"  \
        -ref_genome_version "${HMF_GENOME}" \
        -ref_genome "${HMF_GENOME_FASTA}" \
        -gc_profile "${HMF_GC_PROFILE}" \
        -ensembl_data_dir "${HMF_ENSEMBL_DIR}" \
        -driver_gene_panel "${HMF_DRIVER_PANEL}" \
        -somatic_hotspots "${HMF_SOMATIC_HOTSPOT_VCF}" \
        -germline_hotspots "${HMF_GERMLINE_HOTSPOT_VCF}" \
        -germline_del_freq_file "${HMF_GERMLINE_DEL_FREQ}" \
        -circos "/storage/apps/circos-0.69-10/bin/circos" \
        -output_dir "${PURPLE_OUTPUT}" \
        -threads 1 2>&1 | tee -a "${RunLog}" 

    echo " " >> "${RunLog}"

###-------------------------------------------------------------------------------------------------
       
###---| HMFTOOLS : LINX |---------------------------------------------------------------------------

    echo "###-----/ LINX /--------------------------------------------------------------------------###" >> "${RunLog}"

    LINX_OUTPUT="${HMF_RES_DIR}/linx"
    mkdir -p "${LINX_OUTPUT}"
    
    # somatic 
    ${HMF_RUN_LINX} \
        -sample "${TumorSeqID}" \
        -ref_genome_version "${HMF_GENOME}" \
        -purple_dir "${PURPLE_OUTPUT}" \
        -sv_vcf "${PURPLE_OUTPUT}/${TumorSeqID}.purple.sv.vcf.gz" \
        -output_dir "${LINX_OUTPUT}" \
        -ensembl_data_dir "${HMF_ENSEMBL_DIR}" \
        -known_fusion_file "${HMF_KNOWN_FUSION_CSV}" \
        -driver_gene_panel "${HMF_DRIVER_PANEL}" \
        -write_neo_epitopes \
        -write_all \
        -write_type "SV_DATA;CLUSTER;LINK;ECDNA;CLUSTERING;VIS_DATA" \
        -threads ${Threads} 2>&1 | tee -a "${RunLog}" 

    echo " " >> "${RunLog}"

    # germline
    LINX_GERMLINE_OUTPUT="${HMF_RES_DIR}/linx_germline"
    mkdir -p "${LINX_GERMLINE_OUTPUT}"

    ${HMF_RUN_LINX} \
        -sample "${TumorSeqID}" \
        -germline \
        -ref_genome_version "${HMF_GENOME}" \
        -sv_vcf "${PURPLE_OUTPUT}/${TumorSeqID}.purple.sv.germline.vcf.gz" \
        -output_dir "${LINX_GERMLINE_OUTPUT}" \
        -ensembl_data_dir "${HMF_ENSEMBL_DIR}" \
        -driver_gene_panel "${HMF_DRIVER_PANEL}" \
        -threads ${Threads} 2>&1 | tee -a "${RunLog}" 

    echo " " >> "${RunLog}"


###-------------------------------------------------------------------------------------------------

###---| HMFTOOLS : LILAC |--------------------------------------------------------------------------

    echo "###-----/ LILAC /-------------------------------------------------------------------------###" >> "${RunLog}"

    LILAC_NT_OUTPUT="${HMF_RES_DIR}/lilac/matched_normal"
    LILAC_TONLY_OUTPUT="${HMF_RES_DIR}/lilac/tumor_only"

    mkdir -p "${LILAC_NT_OUTPUT}"
    mkdir -p "${LILAC_TONLY_OUTPUT}"
    
    ### matched-normal
    echo "--- Matched-Normal Mode ---" >> "${RunLog}"
    if [ "${WithRnaseqData}" = "true" ]; then

        ${HMF_RUN_LILAC} \
            -sample "${TumorSeqID}" \
            -tumor_bam "${TUMOR_REDUX_BAM}" \
            -rna_bam "${RnaseqHmftoolsDataDir}/bam/${TumorSeqID}.hmftools.rnaseq.bam" \
            -reference_bam "${NORMAL_REDUX_BAM}" \
            -ref_genome_version "${HMF_GENOME}" \
            -ref_genome "${HMF_GENOME_FASTA}" \
            -resource_dir "${HMF_LILAC_RESOURCE_DIR}" \
            -output_dir "${LILAC_NT_OUTPUT}" \
            -somatic_vcf "${PURPLE_OUTPUT}/${TumorSeqID}.purple.somatic.vcf.gz" \
            -gene_copy_number "${PURPLE_OUTPUT}/${TumorSeqID}.purple.cnv.gene.tsv" \
            -write_types "ALL" \
            -threads ${Threads}  2>&1 | tee -a "${RunLog}"

    else
    
        ${HMF_RUN_LILAC} \
            -sample "${TumorSeqID}" \
            -tumor_bam "${TUMOR_REDUX_BAM}" \
            -reference_bam "${NORMAL_REDUX_BAM}" \
            -ref_genome_version "${HMF_GENOME}" \
            -ref_genome "${HMF_GENOME_FASTA}" \
            -resource_dir "${HMF_LILAC_RESOURCE_DIR}" \
            -output_dir "${LILAC_NT_OUTPUT}" \
            -purple_dir "${PURPLE_OUTPUT}" \
            -somatic_vcf "${PURPLE_OUTPUT}/${TumorSeqID}.purple.somatic.vcf.gz" \
            -gene_copy_number "${PURPLE_OUTPUT}/${TumorSeqID}.purple.cnv.gene.tsv" \
            -write_types "ALL" \
            -threads ${Threads}  2>&1 | tee -a "${RunLog}"

    fi
    ###
    echo " " >> "${RunLog}"
    echo "--- Tumor-Only Mode ---" >> "${RunLog}"
    if [ "${WithRnaseqData}" = "true" ]; then

        ${HMF_RUN_LILAC} \
            -sample "${TumorSeqID}" \
            -tumor_bam "${TUMOR_REDUX_BAM}" \
            -rna_bam "${RnaseqHmftoolsDataDir}/bam/${TumorSeqID}.hmftools.rnaseq.bam" \
            -ref_genome_version "${HMF_GENOME}" \
            -ref_genome "${HMF_GENOME_FASTA}" \
            -resource_dir "${HMF_LILAC_RESOURCE_DIR}" \
            -output_dir "${LILAC_TONLY_OUTPUT}" \
            -somatic_vcf "${PURPLE_OUTPUT}/${TumorSeqID}.purple.somatic.vcf.gz" \
            -gene_copy_number "${PURPLE_OUTPUT}/${TumorSeqID}.purple.cnv.gene.tsv" \
            -write_types "ALL" \
            -threads ${Threads} 2>&1 | tee -a "${RunLog}"

    else

        ${HMF_RUN_LILAC} \
            -sample "${TumorSeqID}" \
            -tumor_bam "${TUMOR_REDUX_BAM}" \
            -ref_genome_version "${HMF_GENOME}" \
            -ref_genome "${HMF_GENOME_FASTA}" \
            -resource_dir "${HMF_LILAC_RESOURCE_DIR}" \
            -output_dir "${LILAC_TONLY_OUTPUT}" \
            -purple_dir "${PURPLE_OUTPUT}" \
            -somatic_vcf "${PURPLE_OUTPUT}/${TumorSeqID}.purple.somatic.vcf.gz" \
            -gene_copy_number "${PURPLE_OUTPUT}/${TumorSeqID}.purple.cnv.gene.tsv" \
            -write_types "ALL" \
            -threads ${Threads} 2>&1 | tee -a "${RunLog}"

    fi

    echo " " >> "${RunLog}"

    if [ -s "${LILAC_TONLY_OUTPUT}/${TumorSeqID}.lilac.reads.tsv" ]; then 
        touch "${HMF_RES_DIR}/lilac/Tumor-Only-Mode_Analyzed"
    else
        touch "${HMF_RES_DIR}/lilac/Tumor-Only-Mode_NOT_Analyzed"
    fi

    if [ -s "${LILAC_NT_OUTPUT}/${TumorSeqID}.lilac.reads.tsv" ]; then 
        touch "${HMF_RES_DIR}/lilac/Matched-Normal-Mode_Analyzed"
    else
        touch "${HMF_RES_DIR}/lilac/Matched-Normal-Mode_NOT_Analyzed"
    fi


###-------------------------------------------------------------------------------------------------

###---| HMFTOOLS : CHORD |--------------------------------------------------------------------------

    echo "###-----/ CHORD /-------------------------------------------------------------------------###" >> "${RunLog}"

    CHORD_OUTPUT="${HMF_RES_DIR}/chord"
    mkdir -p "${CHORD_OUTPUT}"

    ${HMF_RUN_CHORD} \
        -sample "${TumorSeqID}" \
        -output_dir "${CHORD_OUTPUT}" \
        -purple_dir "${PURPLE_OUTPUT}" \
        -ref_genome "${HMF_GENOME_FASTA}" \
        -write_detailed_files \
        -include_non_pass \
        -threads ${Threads} 2>&1 | tee -a "${RunLog}"
       
    echo " " >> "${RunLog}"
 
###-------------------------------------------------------------------------------------------------

###---| HMFTOOLS : SIGS |---------------------------------------------------------------------------

    echo "###-----/ SIGS /--------------------------------------------------------------------------###" >> "${RunLog}"

    SIGS_OUTPUT="${HMF_RES_DIR}/sigs"
    mkdir -p "${SIGS_OUTPUT}"

    ${HMF_RUN_SIGS} \
       -sample "${TumorSeqID}" \
       -output_dir "${SIGS_OUTPUT}" \
       -signatures_file "${HMF_SIGS_SIGNATURE}" \
       -somatic_vcf_file "${PURPLE_OUTPUT}/${TumorSeqID}.purple.somatic.vcf.gz" \
       -write_pos_freq_coords 2>&1 | tee -a "${RunLog}"
           
    echo " " >> "${RunLog}"
 
###-------------------------------------------------------------------------------------------------

###---| HMFTOOLS : PEACH |--------------------------------------------------------------------------

    echo "###-----/ PEACH /-------------------------------------------------------------------------###" >> "${RunLog}"

    PEACH_OUTPUT="${HMF_RES_DIR}/peach"
    mkdir -p "${PEACH_OUTPUT}"

    ${HMF_RUN_PEACH} \
       -sample_name "${NormalSeqID}" \
       -haplotypes_file "${HMF_PEACH_HAPLOTYPE}" \
       -vcf_file "${PURPLE_OUTPUT}/${TumorSeqID}.purple.germline.vcf.gz" \
       -function_file "${HMF_HAPLOTYPE_FUNCTION}" \
       -drugs_file "${HMF_PEACH_DRUGS}" \
       -output_dir "${PEACH_OUTPUT}" 2>&1 | tee -a "${RunLog}"
            
    echo " " >> "${RunLog}"
 
###-------------------------------------------------------------------------------------------------

###---| HMFTOOLS : VIRUS-BREAKEND |-----------------------------------------------------------------

    echo "###-----/ VIRUS-BREAKEND /----------------------------------------------------------------###" >> "${RunLog}"

    VIRUSBREAKEND_OUTPUT="${HMF_RES_DIR}/virus-breakend"
    mkdir -p "${VIRUSBREAKEND_OUTPUT}"
    
    cd "${VIRUSBREAKEND_OUTPUT}"
    export PATH=$PATH:/storage/apps/gridss-v2.13.2:/storage/apps/kraken2:/storage/apps/RepeatMasker-4.1.8:/storage/apps/rmblast-2.14.0/bin

    ${HMF_RUN_VIRUS_BREAKEND} \
        --reference "${HMF_GENOME_FASTA}" \
        --jar "${HMF_GRIDSS_JAR}" \
        --output "${VIRUSBREAKEND_OUTPUT}/${TumorSeqID}.virus.breakend.vcf" \
        --db "${HMF_VIRUSBREAKEND_DB}" \
        --threads "${Threads}" \
        --workingdir "${VIRUSBREAKEND_OUTPUT}" \
        "${TUMOR_REDUX_BAM}" 2>&1 | tee -a "${RunLog}"

    echo " " >> "${RunLog}"
    echo "--- Cross-Check with Matched-Normal ---" >> "${RunLog}"

    VIRUSBREAKEND_OUTPUT_NORMAL="${HMF_RES_DIR}/virus-breakend_matched_normal"
    mkdir -p "${VIRUSBREAKEND_OUTPUT_NORMAL}"
    cd "${VIRUSBREAKEND_OUTPUT_NORMAL}"

    ${HMF_RUN_VIRUS_BREAKEND} \
        --reference "${HMF_GENOME_FASTA}" \
        --jar "${HMF_GRIDSS_JAR}" \
        --output "${VIRUSBREAKEND_OUTPUT_NORMAL}/${NormalSeqID}.virus.breakend.vcf" \
        --db "${HMF_VIRUSBREAKEND_DB}" \
        --threads "${Threads}" \
        --workingdir "${VIRUSBREAKEND_OUTPUT_NORMAL}" \
        "${NORMAL_REDUX_BAM}" 2>&1 | tee -a "${RunLog}"

###-------------------------------------------------------------------------------------------------

###---| HMFTOOLS : VIRUS-INTEPRETER |---------------------------------------------------------------

    echo "###-----/ VIRUS-INTERPRETER /--------------------------------------------------------------###" >> "${RunLog}"

    VIRUSINTERPRETER_OUTPUT="${HMF_RES_DIR}/virus-interpreter"
    mkdir -p "${VIRUSINTERPRETER_OUTPUT}"

    ${HMF_RUN_VIRUS_INTERPRETATER} \
        -sample "${TumorSeqID}" \
        -output_dir "${VIRUSINTERPRETER_OUTPUT}" \
        -purple_dir "${PURPLE_OUTPUT}" \
        -taxonomy_db_tsv "${HMF_VIRUS_TAX_DB}" \
        -tumor_sample_wgs_metrics_file "${REDUX_OUTPUT}/${TumorSeqID}.redux.wgs.metrics.txt" \
        -virus_blacklisting_db_tsv "${HMF_VIRUS_BLACKLIST_DB}" \
        -virus_breakend_tsv "${VIRUSBREAKEND_OUTPUT}/${TumorSeqID}.virus.breakend.vcf.summary.tsv"  \
        -virus_reporting_db_tsv "${HMF_VIRUS_REPORT_DB}" 2>&1 | tee -a "${RunLog}"

    echo " " >> "${RunLog}"


###-------------------------------------------------------------------------------------------------

###---| HMFTOOLS : CUPPA |--------------------------------------------------------------------------

    echo "###-----/ CUPPA /-------------------------------------------------------------------------###" >> "${RunLog}"

    CUPPA_OUTPUT="${HMF_RES_DIR}/cuppa"
    mkdir -p "${CUPPA_OUTPUT}"

    if [ "${WithRnaseqData}" = "true" ]; then

        ${HMF_RUN_CUPPA} \
            -sample "${TumorSeqID}" \
            -output_dir "${CUPPA_OUTPUT}" \
            -categories ALL \
            -classifier_path "${HMF_CUPPA_PICKLE}" \
            -linx_dir "${LINX_OUTPUT}" \
            -purple_dir "${PURPLE_OUTPUT}" \
            -isofox_dir "${RnaseqHmftoolsDataDir}/isofox" \
            -rna_sample "${TumorSeqID}" \
            -python_path "/opt/conda/bin/python" \
            -ref_alt_sj_sites "${HMF_CUPPA_SJ_SITES}" \
            -ref_genome_version "${HMF_GENOME}" \
            -virus_dir "${VIRUSINTERPRETER_OUTPUT}" \
            -threads "${Threads}"  2>&1 | tee -a "${RunLog}"

    else

        ${HMF_RUN_CUPPA} \
            -sample "${TumorSeqID}" \
            -output_dir "${CUPPA_OUTPUT}" \
            -categories DNA \
            -classifier_path "${HMF_CUPPA_PICKLE}" \
            -linx_dir "${LINX_OUTPUT}" \
            -purple_dir "${PURPLE_OUTPUT}" \
            -python_path "/opt/conda/bin/python" \
            -ref_alt_sj_sites "${HMF_CUPPA_SJ_SITES}" \
            -ref_genome_version "${HMF_GENOME}" \
            -virus_dir "${VIRUSINTERPRETER_OUTPUT}" \
            -threads "${Threads}"  2>&1 | tee -a "${RunLog}"
    
    fi

    echo " " >> "${RunLog}"
 
###-------------------------------------------------------------------------------------------------

###---| NEO |---------------------------------------------------------------------------------------

    echo "###-----/ NEO /---------------------------------------------------------------------------###" >> "${RunLog}"

    NEO_OUTPUT="${HMF_RES_DIR}/neo"
    mkdir -p "${NEO_OUTPUT}"

    ${HMF_RUN_NEO} \
        -sample "${TumorSeqID}" \
        -output_dir "${NEO_OUTPUT}" \
        -ref_genome_version "${HMF_GENOME}" \
        -ref_genome "${HMF_GENOME_FASTA}" \
        -linx_dir "${LINX_OUTPUT}" \
        -somatic_vcf "${PURPLE_OUTPUT}/${TumorSeqID}.purple.somatic.vcf.gz" \
        -ensembl_data_dir "${HMF_ENSEMBL_DIR}" \
        -write_trans_data \
        -threads "${Threads}"  2>&1 | tee -a "${RunLog}"

    echo " " >> "${RunLog}"

###-------------------------------------------------------------------------------------------------

###---| TEAL |--------------------------------------------------------------------------------------

    echo "###-----/ TEAL /--------------------------------------------------------------------------###" >> "${RunLog}"

    TEAL_OUTPUT="${HMF_RES_DIR}/teal"
    mkdir -p "${TEAL_OUTPUT}"

    ${HMF_RUN_TEAL} \
        -tumor "${TumorSeqID}" \
        -tumor_bam "${TUMOR_REDUX_BAM}" \
        -reference "${NormalSeqID}" \
        -reference_bam "${NORMAL_REDUX_BAM}" \
        -output_dir "${TEAL_OUTPUT}" \
        -ref_genome_version "${HMF_GENOME}" \
        -ref_genome "${HMF_GENOME_FASTA}" \
        -reference_wgs_metrics "${NORMAL_BAMTOOLS_OUTPUT}/${NormalSeqID}.bam_metric.summary.tsv" \
        -tumor_wgs_metrics "${TUMOR_BAMTOOLS_OUTPUT}/${TumorSeqID}.bam_metric.summary.tsv" \
        -purple "${PURPLE_OUTPUT}" \
        -cobalt "${COBALT_OUTPUT}" \
        -threads "${Threads}"  2>&1 | tee -a "${RunLog}"

    echo " " >> "${RunLog}"

###-------------------------------------------------------------------------------------------------

###---| ORANGE |------------------------------------------------------------------------------------

    echo "###-----/ ORANGE /--------------------------------------------------------------------------###" >> "${RunLog}"

    ORANGE_OUTPUT="${HMF_RES_DIR}/orange"
    mkdir -p "${ORANGE_OUTPUT}"

    HMFTOOLS_RESULT_DIR_ARGS=("-purple_dir ${PURPLE_OUTPUT}")
    if [[ -d "${PURPLE_OUTPUT}/plot" ]]; then HMFTOOLS_RESULT_DIR_ARGS+=("-purple_plot_dir ${PURPLE_OUTPUT}/plot"); fi
    if [[ -d "${CHORD_OUTPUT}" ]]; then HMFTOOLS_RESULT_DIR_ARGS+=("-chord_dir ${CHORD_OUTPUT}"); fi
    if [[ -d "${CUPPA_OUTPUT}" ]]; then HMFTOOLS_RESULT_DIR_ARGS+=("-cuppa_dir ${CUPPA_OUTPUT}"); fi
    if [[ -f "${HMF_RES_DIR}/lilac/Matched-Normal-Mode_Analyzed" ]]; then
        HMFTOOLS_RESULT_DIR_ARGS+=("-lilac_dir ${HMF_RES_DIR}/lilac/matched_normal")       
    elif [[ -f "${HMF_RES_DIR}/lilac/Tumor-Only-Mode_Analyzed" ]]; then
        HMFTOOLS_RESULT_DIR_ARGS+=("-lilac_dir ${HMF_RES_DIR}/lilac/tumor_only")   
    fi
    if [[ -d "${LINX_OUTPUT}" ]]; then HMFTOOLS_RESULT_DIR_ARGS+=("-linx_dir ${LINX_OUTPUT}"); fi
    if [[ -d "${LINX_GERMLINE_OUTPUT}" ]]; then HMFTOOLS_RESULT_DIR_ARGS+=("-linx_germline_dir ${LINX_GERMLINE_OUTPUT}"); fi
    if [[ -d "${PEACH_OUTPUT}" ]]; then HMFTOOLS_RESULT_DIR_ARGS+=("-peach_dir ${PEACH_OUTPUT}"); fi
    if [[ -d "${SAGE_OUTPUT}" ]]; then HMFTOOLS_RESULT_DIR_ARGS+=("-sage_dir ${SAGE_OUTPUT}"); fi
    if [[ -d "${SAGE_GERMLINE_OUTPUT}" ]]; then HMFTOOLS_RESULT_DIR_ARGS+=("-sage_germline_dir ${SAGE_GERMLINE_OUTPUT}"); fi
    if [[ -d "${SIGS_OUTPUT}" ]]; then HMFTOOLS_RESULT_DIR_ARGS+=("-sigs_dir ${SIGS_OUTPUT}"); fi
    if [[ -d "${VIRUSINTERPRETER_OUTPUT}" ]]; then HMFTOOLS_RESULT_DIR_ARGS+=("-virus_dir ${VIRUSINTERPRETER_OUTPUT}"); fi
    if [[ -d "${TUMOR_BAMTOOLS_OUTPUT}" ]]; then HMFTOOLS_RESULT_DIR_ARGS+=("-tumor_metrics_dir ${TUMOR_BAMTOOLS_OUTPUT}"); fi
    if [[ -d "${NORMAL_BAMTOOLS_OUTPUT}" ]]; then HMFTOOLS_RESULT_DIR_ARGS+=("-ref_metrics_dir ${NORMAL_BAMTOOLS_OUTPUT}"); fi

    ${HMF_RUN_ORANGE} \
        -experiment_type "WGS" \
        -primary_tumor_doids "${PrimaryTumorDoid}" \
        -tumor_sample_id "${TumorSeqID}" \
        -reference_sample_id "${NormalSeqID}" \
        -output_dir "${ORANGE_OUTPUT}" \
        -cohort_mapping_tsv "${HMF_COHORT_MAPPING}" \
        -cohort_percentiles_tsv "${HMF_COHORT_PERCENTILES}" \
        -doid_json "${HMF_DOID_JSON}" \
        -ref_genome_version "${HMF_GENOME}" \
        -driver_gene_panel "${HMF_DRIVER_PANEL}" \
        -ensembl_data_dir "${HMF_ENSEMBL_DIR}" \
        -known_fusion_file "${HMF_KNOWN_FUSION_CSV}" \
        -signatures_etiology_tsv "${HMF_SIGATURE_ETIOLOGY}" \
        ${HMFTOOLS_RESULT_DIR_ARGS[@]}

    echo " " >> "${RunLog}"

###-------------------------------------------------------------------------------------------------

###---| PURPLE VARIANTS VEP ANNOTATION |------------------------------------------------------------

    echo "###-----/ PURPLE SOMATIC VARIANTS ANNOTATION /--------------------------------------------###" >> "${RunLog}"

    # input somatic variants vcf
    InputPurpleSomaticVCF="${PURPLE_OUTPUT}/${TumorSeqID}.purple.somatic.vcf.gz"

    # output files prefix
    OutputVCF_Refseq="${InputPurpleSomaticVCF/vcf.gz/vep.refseq.vcf}"
    OutputVCF_Ensembl="${InputPurpleSomaticVCF/vcf.gz/vep.ensembl.vcf}"
    OutputVCF_Merged="${InputPurpleSomaticVCF/vcf.gz/vep.annot.vcf}"
    OutputVCF_Normalized="${InputPurpleSomaticVCF/vcf.gz/vep.annotation.vcf}"
    OutputMAF_Merged="${InputPurpleSomaticVCF/vcf.gz/vep.annotation.maf}"

    if [ "${GenomeAssembly}" = "hg38" ]; then
        VepGenome="GRCh38"
    else
        VepGenome="GRCh37"
    fi

    # refseq
    ${RunVep} --force_overwrite --offline \
        --input_file "${InputPurpleSomaticVCF}" --vcf \
        --output_file "${OutputVCF_Refseq}" \
        --assembly "${VepGenome}" --species "homo_sapiens" \
        --cache --dir_cache "${VepCacheDir}" --dir_plugins "${VepPluginDir}" \
        --bam "${VepBam}" --fasta "${HMF_GENOME_FASTA}" \
        --stats_text --show_ref_allele --refseq --use_transcript_ref --hgvs --canonical --pick --mane_select \
        --fork "${Threads}" --buffer_size 50000 

    # ensembl 
    ${RunVep} --force_overwrite --offline \
        --input_file "${InputPurpleSomaticVCF}" --vcf \
        --output_file "${OutputVCF_Ensembl}" \
        --assembly "${VepGenome}" --species "homo_sapiens" \
        --cache --dir_cache "${VepCacheDir}" --dir_plugins "${VepPluginDir}" \
        --bam ${VepBam} --fasta ${HMF_GENOME_FASTA} \
        --custom file=${VepData_Clinvar},short_name=ClinVar,format=vcf,fields=CLNSIG%CLNDN%ORIGIN \
        --custom file=${VepData_Cosmic},short_name=COSMIC,format=vcf,fields=HGVSC%LEGACY_ID%TIER \
        --plugin AlphaMissense,file=${VepData_AlphaMissense} \
        --stats_text --variant_class --canonical --hgvs --hgvsg --max_af --af_1kg --af_gnomade \
        --ccds --uniprot --sift b --polyphen b --gene_phenotype --regulatory --biotype --symbol \
        --numbers --domains --exclude_predicted --pick --mane_select \
        --fork ${Threads} --buffer_size 50000 

    bcftools +split-vep "${OutputVCF_Refseq}" \
        -c "Feature,Gene,EXON,INTRON,HGVSc,HGVSp" -p "Refseq_" -Ou \
        | bcftools annotate -x "INFO/CSQ" \
        -Oz -o "${PURPLE_OUTPUT}/refseq.only.annot.vcf.gz"
    bcftools index -t "${PURPLE_OUTPUT}/refseq.only.annot.vcf.gz"

    bgzip -c "${OutputVCF_Ensembl}" > "${OutputVCF_Ensembl}.gz"
    bcftools index -t "${OutputVCF_Ensembl}.gz"
    bcftools annotate \
        -a "${PURPLE_OUTPUT}/refseq.only.annot.vcf.gz" \
        -c "INFO/Refseq_Feature,INFO/Refseq_Gene,INFO/Refseq_EXON,INFO/Refseq_INTRON,INFO/Refseq_HGVSc,INFO/Refseq_HGVSp" \
        -Ov -o "${OutputVCF_Merged}" \
        "${OutputVCF_Ensembl}.gz"

    rm "${PURPLE_OUTPUT}/refseq.only.annot.vcf.gz" "${PURPLE_OUTPUT}/refseq.only.annot.vcf.gz.tbi" "${OutputVCF_Ensembl}.gz" "${OutputVCF_Ensembl}.gz.tbi"

    bgzip -c "${OutputVCF_Merged}" > "${OutputVCF_Merged}.gz"
    tabix -p vcf -f "${OutputVCF_Merged}.gz"
    
    bcftools norm -f "${HMF_GENOME_FASTA}" -m -any "${OutputVCF_Merged}.gz" -Ov -o "${OutputVCF_Normalized}"

    RetainInfoItems="Refseq_Feature,Refseq_Gene,Refseq_EXON,Refseq_INTRON,Refseq_HGVSc,Refseq_HGVSp,DP,GERMQ,MBQ,MFRL,MMQ,MPOS,POPAF,PON,ROQ,TLOD,FILTER,BLACKLIST_BED,BLACKLIST_VCF,HOTSPOT,MSG,TIER"

    ${RunVcf2Maf} --inhibit-vep \
        --input-vcf "${OutputVCF_Normalized}" \
        --output-maf "${OutputMAF_Merged}" \
        --tumor-id "${TumorSeqID}" \
        --normal-id "${NormalSeqID}" \
        --ref-fasta "${HMF_GENOME_FASTA}" \
        --ncbi-build "${VepGenome}" \
        --retain-info "${RetainInfoItems}" \
        --retain-ann "HGVSg,am_class,am_pathogenicity,ClinVar_CLNSIG,ClinVar_CLNDN,ClinVar_ORIGIN,gnomADe_AF,gnomADe_EAS_AF,MAX_AF,MAX_AF_POPS,TRANSCRIPTION_FACTORS,COSMIC,COSMIC_HGVSC,COSMIC_LEGACY_ID,COSMIC_TIER"

###-------------------------------------------------------------------------------------------------


# ###---| HMFTOOLS : WISP |--------------------------------------------------------------------------- <<<< UNDER CONSTRUCTION

#     echo "###-----/ WISP /--------------------------------------------------------------------------###" >> "${RunLog}"

#     WISP_OUTPUT=${HMF_RES_DIR}/wisp
#     mkdir -p ${WISP_OUTPUT}

#     ${HMF_RUN_WISP} \
#         -patient_id PATIENT_ID
#         -tumor_id TUMOR_ID
#         -samples SAMPLE_ID
#         -purple_dir /sample_data/purple/
#         -amber_dir ${AMBER_OUTPUT} \
#         -cobalt_dir ${COBALT_OUTPUT} \
#         -somatic_vcf /sample_data/sage_append/TUMOR_ID.purple.somatic.vcf.gz
#         -bqr_dir ${SAGE_OUTPUT} \
#         -ref_genome_version 38 -ref_genome ${HMF_GENOME_FASTA} \
#         -output_dir /path_to_output/
###-------------------------------------------------------------------------------------------------




