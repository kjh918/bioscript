
#!/bin/bash

# script version v0.5

#---| Input Options |-------------------------------------------------------------------------------
    usage() {
        echo "
        Usage: WES_M04-1_VARIANT-ANNOTATION.sh [ ARGUMENTS ]...

        --BASE_DIR                  analysis base folder by NGS data type
                                    default = /data/wes
        
        --BATCH_ID                  [ REQUIRED ] dataset batch-id
        
        --TUMOR_SEQ_ID              [ REQURIED ] tumor-sample-id
        
        --NORMAL_SEQ_ID             normal-sample-id. 
                                    Required when matched-normal vcf annotation

        --ANNOTATE_TARGET_VCF_TYPE  [ REQUIRED ] annotation target vcf type by variant call method. 
                                    Tonly = tumor only call
                                    NT_TS = matched normal call using normal tissue 
                                    NT_ORG = matched normal call using normal organoid
                                    default = 'Tonly'

        --GENOME_ASSEMBLY           reference genome assembly version
                                    default = hg19

        --THREADS                   parallel processing cpu cores 
                                    default = 16

        -h | --help                 Print this usage 
        "
    }
#---------------------------------------------------------------------------------------------------

#---| Parse Option Arguments |----------------------------------------------------------------------
    ARGS=$(getopt -a -o h: --long BASE_DIR:,BATCH_ID:,TUMOR_SEQ_ID:,NORMAL_SEQ_ID:,ANNOTATE_TARGET_VCF_TYPE:,THREADS:,GENOME_ASSEMBLY:,help -- "$@" )
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
                -* | --* | "") echo "no value at $1, default = null" ; shift ;;
                *)  NormalSeqID=$2 ; shift 2 ;;
            esac  ;;
        --ANNOTATE_TARGET_VCF_TYPE )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  AnnotateTargetVCF=$2 ; shift 2 ;;
            esac ;;
        --THREADS )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 16" ; shift ;;
                *)  Threads=$2 ; shift 2 ;;
            esac ;;
        --GENOME_ASSEMBLY )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 'hg19'" ; shift ;;
                *)  GenomeAssembly=$2 ; shift 2 ;;
            esac ;;
        -h | --help )
        usage >&2 ; exit 2 ;;   
        --) shift ; break ;;
        *)  usage >&2 ; exit 2 ;;
    esac
    done
#---------------------------------------------------------------------------------------------------

#---| Check Params and Set Default |----------------------------------------------------------------
    if [[ -z "${BaseDir}" ]]; then BaseDir="/data/wes"; fi
    if [[ -z "${Threads}" ]]; then Threads=16; fi
    if [[ -z "${GenomeAssembly}" ]]; then GenomeAssembly="hg19"; fi
#---------------------------------------------------------------------------------------------------

#---| Manual Params |-------------------------------------------------------------------------------
    # BaseDir="/data/wgs"
    # BatchID="WGS-26-01"
    # TumorSeqID="WGS-26-01_01"
    # NormalSeqID="WGS-26-01_03"
    # AnnotateTargetVCF="NT_TS"
    # Threads=16
    # GenomeAssembly="hg19"
#---------------------------------------------------------------------------------------------------

#---| FOLDERS |-------------------------------------------------------------------------------------
    LogDir=${BaseDir}/${BatchID}/${TumorSeqID}/log
    VcfDir=${BaseDir}/${BatchID}/${TumorSeqID}/vcf
    if [ ! -d ${LogDir} ]; then mkdir -p ${LogDir}; fi
#---------------------------------------------------------------------------------------------------

#---| SWTOOLS RUN |---------------------------------------------------------------------------------   
    RunVep="singularity exec -B /storage,/data /storage/images/vep-110.1.sif /opt/vep/src/ensembl-vep/vep"
    RunVcf2Maf="perl /storage/apps/vcf2maf-1.6.22/vcf2maf.pl"
#---------------------------------------------------------------------------------------------------

#---| LogFile |-------------------------------------------------------------------------------------
    LogFile=${LogDir}/${TumorSeqID}.vep.variant.annotation.$(date '+%Y%m%d').log
    if [ ! -f ${LogFile} ]; then touch ${LogFile}; fi
#---------------------------------------------------------------------------------------------------

#---| RESOURCES |-----------------------------------------------------------------------------------
    GenomeFasta=/storage/references_and_index/${GenomeAssembly}/fasta/${GenomeAssembly}.primary_decoy_HLA_no.ALT.fa
    VepCacheDir=/storage/references_and_index/${GenomeAssembly}/vep
    VepPluginDir=${VepCacheDir}/Plugins/110
    VepBam=${VepCacheDir}/bam/${GenomeAssembly}.vep.bam
    VepTopLevelFasta=${VepCacheDir}/toplevel_fasta/${GenomeAssembly}.toplevel.fa.gz
    VepData_Clinvar=/storage/references_and_index/${GenomeAssembly}/annotation/clinvar/clinvar.vcf.gz
    VepData_AlphaMissense=/storage/references_and_index/${GenomeAssembly}/annotation/AlphaMissense/alpha-missense.tsv.gz
    VepData_Cosmic=/storage/references_and_index/${GenomeAssembly}/annotation/cosmic/Cosmic_GenomeScreensMutant_v103.vcf.gz
#---------------------------------------------------------------------------------------------------

#---| INPUT VCF TYPE |------------------------------------------------------------------------------
    shopt -s nullglob

    if [[ "${AnnotateTargetVCF}" = "NT_TS" ]]; then 
        INPUT_VCF=${VcfDir}/${TumorSeqID}.NT.VARIANTS.vcf.gz
    elif [[ "${AnnotateTargetVCF}" = "NT_ORG" ]]; then 
        INPUT_VCF=${VcfDir}/${TumorSeqID}.NT.VARIANTS.vcf.gz
    else
        INPUT_VCF=${VcfDir}/${TumorSeqID}.Tonly.VARIANTS.vcf.gz
    fi
  
    if [[ "${AnnotateTargetVCF}" = "NT_TS" || "${AnnotateTargetVCF}" = "NT_ORG" ]]; then
        if [[ -z "${NormalSeqID}" ]]; then 
            echo " 'NORMAL_SEQ_ID' is REQUIRED when 'VARAINT_CALL_METHD' = 'NT' " >&2 ; exit 2
        fi
    fi

#---------------------------------------------------------------------------------------------------

#---| OUTPUT FILENAME |-----------------------------------------------------------------------------
    OutputVCF_Refseq="${INPUT_VCF/.vcf.gz/.refseq.vcf}"
    OutputVCF_Ensembl="${INPUT_VCF/.vcf.gz/.ensembl.vcf}"
    OutputVCF_Merged="${INPUT_VCF/.vcf.gz/.vep.vcf}"
    OutputVCF_Normalized="${INPUT_VCF/.vcf.gz/.VEP.ANNOTATED.vcf}"
    OutputMAF_Merged="${INPUT_VCF/.vcf.gz/.VEP.ANNOTATED.maf}"
#---------------------------------------------------------------------------------------------------

#=================================================================================================== 
    echo " " >> ${LogFile}
    echo "WES VariantAnnotation log" >> ${LogFile}
    echo "----------------------------------------------------------------------" >> ${LogFile}
    echo " RUN DATE             : $(date '+%Y-%m-%d %H:%M:%S')" >> ${LogFile}
    echo " BATCH ID             : ${BatchID}" >> ${LogFile}
    echo " SEQ ID               : ${TumorSeqID}" >> ${LogFile}
    echo " GENOME ASSEMBLY      : hg19" >> ${LogFile}
    echo "----------------------------------------------------------------------" >> ${LogFile}
    echo " " >> ${LogFile}
#=================================================================================================== 

#---| GENOME TAG for VEP |--------------------------------------------------------------------------
    if [ "${GenomeAssembly}" = "hg38" ]; then
        VepGenome="GRCh38"
    else
        VepGenome="GRCh37"
    fi
#---------------------------------------------------------------------------------------------------

#=================================================================================================== 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | VARIANT ANNOTATION START." >> ${LogFile}
#=================================================================================================== 

#---| Run VEP : RefSeq |----------------------------------------------------------------------------
    NumVars=$(zgrep -vc '^#' ${INPUT_VCF})
    
    if (( ${NumVars} > 0 )); then
        
        # refseq
        ${RunVep} --force_overwrite --offline \
            --input_file ${INPUT_VCF} --vcf \
            --output_file ${OutputVCF_Refseq} \
            --stats_text --show_ref_allele --species homo_sapiens --assembly ${VepGenome} \
            --cache --dir_cache ${VepCacheDir} --dir_plugins ${VepPluginDir} \
            --refseq --use_transcript_ref --hgvs --canonical --pick --mane_select \
            --bam ${VepBam} --fasta ${GenomeFasta} \
            --fork ${Threads} --buffer_size 50000 

        # ensembl 
        ${RunVep} --force_overwrite --offline \
            --input_file ${INPUT_VCF} --vcf \
            --output_file ${OutputVCF_Ensembl} \
            --stats_text --species homo_sapiens --assembly ${VepGenome} \
            --cache --dir_cache ${VepCacheDir} --dir_plugins ${VepPluginDir} \
            --variant_class --canonical --hgvs --hgvsg --max_af --af_1kg --af_gnomade \
            --ccds --uniprot --sift b --polyphen b --gene_phenotype --regulatory --biotype --symbol \
            --numbers --domains --exclude_predicted --pick --mane_select \
            --custom file=${VepData_Clinvar},short_name=ClinVar,format=vcf,fields=CLNSIG%CLNDN%ORIGIN \
            --custom file=${VepData_Cosmic},short_name=COSMIC,format=vcf,fields=HGVSC%LEGACY_ID%TIER \
            --plugin AlphaMissense,file=${VepData_AlphaMissense} \
            --bam ${VepBam} --fasta ${GenomeFasta} \
            --fork ${Threads} --buffer_size 50000 

        bcftools +split-vep ${OutputVCF_Refseq} \
            -c Feature,Gene,EXON,INTRON,HGVSc,HGVSp -p Refseq_ -Ou \
            | bcftools annotate -x INFO/CSQ \
            -Oz -o ${VcfDir}/refseq.add.vcf.gz
        bcftools index -t ${VcfDir}/refseq.add.vcf.gz

        bgzip -c ${OutputVCF_Ensembl} > ${OutputVCF_Ensembl}.gz
        bcftools index -t ${OutputVCF_Ensembl}.gz
        
        bcftools annotate \
            -a ${VcfDir}/refseq.add.vcf.gz \
            -c INFO/Refseq_Feature,INFO/Refseq_Gene,INFO/Refseq_EXON,INFO/Refseq_INTRON,INFO/Refseq_HGVSc,INFO/Refseq_HGVSp \
            -Ov -o ${OutputVCF_Merged} \
            ${OutputVCF_Ensembl}.gz

        #rm ${VcfDir}/refseq.add.vcf.gz ${VcfDir}/refseq.add.vcf.gz.tbi ${OutputVCF_Ensembl}.gz ${OutputVCF_Ensembl}.gz.tbi

        bgzip -c ${OutputVCF_Merged} > ${OutputVCF_Merged}.gz
        tabix -p vcf -f ${OutputVCF_Merged}.gz
        bcftools norm -f ${GenomeFasta} -m -any ${OutputVCF_Merged}.gz -Ov -o ${OutputVCF_Normalized}

        RetainInfoItems="Refseq_Feature,Refseq_Gene,Refseq_EXON,Refseq_INTRON,Refseq_HGVSc,Refseq_HGVSp,DP,GERMQ,MBQ,MFRL,MMQ,MPOS,POPAF,PON,ROQ,TLOD,FILTER,GNOMAD_AF,GNOMAD_DB"

        VCF2MAF_INPUT_ARGS=("--tumor-id ${TumorSeqID}")
        if [[ "${AnnotateTargetVCF}" = "NT_TS" || "${AnnotateTargetVCF}" = "NT_ORG" ]]; then
            VCF2MAF_INPUT_ARGS+=("--normal-id ${NormalSeqID}")
        fi

        ${RunVcf2Maf} --inhibit-vep \
            --input-vcf "${OutputVCF_Normalized}" \
            --output-maf "${OutputMAF_Merged}" \
            ${VCF2MAF_INPUT_ARGS[@]} \
            --ref-fasta "${GenomeFasta}" \
            --retain-info "${RetainInfoItems}" \
            --ncbi-build "${VepGenome}" \
            --retain-ann "HGVSg,am_class,am_pathogenicity,ClinVar_CLNSIG,ClinVar_CLNDN,ClinVar_ORIGIN,gnomADe_AF,gnomADe_EAS_AF,MAX_AF,MAX_AF_POPS,TRANSCRIPTION_FACTORS,COSMIC,COSMIC_HGVSC,COSMIC_LEGACY_ID,COSMIC_TIER"



    else
        touch ${OutputVCF_Merged} 
    fi
#---------------------------------------------------------------------------------------------------

#---| VCF Normlization |----------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------


#=================================================================================================== 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | VEP-ANNOTATION FINISHED." >> ${LogFile}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | VCF to MAF START." >> ${LogFile}
#=================================================================================================== 

#---| VCF to MAF |----------------------------------------------------------------------------------
    
    if (( ${NumVars} > 0 )); then

        if [ "${PurecnAddedVcf}" = "true" ]; then
            RetainInfoItems="Refseq_Feature,Refseq_Gene,Refseq_EXON,Refseq_INTRON,Refseq_HGVSc,Refseq_HGVSp,DP,GERMQ,MBQ,MFRL,MMQ,MPOS,POPAF,PON,ROQ,TLOD,FILTER,MLK,NUM_TOOLS,SOMATIC,PureCN.PS,PureCN.ML.SOMATIC,PureCN.ML.C,PureCN.CS,PureCN.ML.LOH,PureCN.GHOMOZYGOUS,PureCN.LR,PureCN.CF,PureCN.ML.M"
        else
            RetainInfoItems="Refseq_Feature,Refseq_Gene,Refseq_EXON,Refseq_INTRON,Refseq_HGVSc,Refseq_HGVSp,DP,GERMQ,MBQ,MFRL,MMQ,MPOS,POPAF,PON,ROQ,TLOD,FILTER"
        fi

        if [[ "${AnnotateTargetVCF}" = "NT_TS" || "${AnnotateTargetVCF}" = "NT_ORG" ]]; then
            
            ${RunVcf2Maf} --inhibit-vep \
                --input-vcf ${OutputVCF_Normalized} \
                --output-maf ${OutputMAF_Merged} \
                --tumor-id ${TumorSeqID} \
                --normal-id ${NormalSeqID} \
                --ref-fasta ${GenomeFasta} \
                --retain-info ${RetainInfoItems} \
                --retain-ann HGVSg,am_class,am_pathogenicity,ClinVar_CLNSIG,ClinVar_CLNDN,ClinVar_ORIGIN,gnomADe_AF,gnomADe_EAS_AF,MAX_AF,MAX_AF_POPS,TRANSCRIPTION_FACTORS,COSMIC,COSMIC_HGVSC,COSMIC_LEGACY_ID,COSMIC_TIER
            
        else

            ${RunVcf2Maf} --inhibit-vep \
                --input-vcf ${OutputVCF_Merged} \
                --output-maf ${OutputMAF_Merged} \
                --tumor-id ${TumorSeqID} \
                --ref-fasta ${GenomeFasta} \
                --retain-info ${RetainInfoItems} \
                --retain-ann HGVSg,am_class,am_pathogenicity,ClinVar_CLNSIG,ClinVar_CLNDN,ClinVar_ORIGIN,gnomADe_AF,gnomADe_EAS_AF,MAX_AF,MAX_AF_POPS,TRANSCRIPTION_FACTORS,COSMIC,COSMIC_HGVSC,COSMIC_LEGACY_ID,COSMIC_TIER
        fi
    else
        touch ${OutputMAF_Merged}
    fi
#---------------------------------------------------------------------------------------------------

#=================================================================================================== 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | VCF to MAF FINISHED." >> ${LogFile}
    echo " " >> ${LogFile}
    echo "WES VariantAnnotation DONE." >> ${LogFile}
    echo "----------------------------------------------------------------------" >> ${LogFile}
    echo " " >> ${LogFile}
#=================================================================================================== 

