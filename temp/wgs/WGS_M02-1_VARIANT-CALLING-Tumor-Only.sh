
#!/bin/bash

# SCRIPT VERSION : v0.1
# GUIDELINE ID   : GMC-NGS-B02-A09
# DATE           : 2024-03-11
# AUTHOR         : kangsm@gencurix.com

#---| Input Options |------------------------------------------------------------------------------#
    usage() {echo "
        Usage: WES_M02-1_MUTECT2-Tumor-Only.sh [ ARGUMENTS ]...
        
        --BASE_DIR          analysis base folder by NGS data type
                            default = /data/wgs
        
        --BATCH_ID          [ REQUIRED ] dataset batch-id
        
        --TUMOR_SEQ_ID      [ REQURIED ] tumor sample-id 
                
        --THREADS           parallel processing cpu cores 
                            default = 16
        
        --GENOME_ASSEMBLY   Reference genome assembly version
                            default = hghg19

        --WES_LIB           WES capture library 
                            default = twist.exome.2.0

        -h | --help ]       Print this usage 
        "
    }
#--------------------------------------------------------------------------------------------------#

#---| Parse Option Arguments |---------------------------------------------------------------------#
    ARGS=$(getopt -a -o h: --long BATCH_ID:,TUMOR_SEQ_ID:,WES_LIB:,BASE_DIR:,THREADS:,GENOME_ASSEMBLY:,NORMAL_TYPE:,help -- "$@" )
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
        -h | --help )
        usage >&2 ; exit 2 ;;   
        --) shift ; break ;;
        *)  usage >&2 ; exit 2 ;;
    esac
    done
#--------------------------------------------------------------------------------------------------#

#---| Check Params and Set Default |----------------------------------------------------------------
    if [ -z "${BaseDir}" ]; then BaseDir="/data/wes"; fi
    if [ -z "${Threads}" ]; then Threads=16; fi
    if [ -z "${GenomeAssembly}" ]; then GenomeAssembly="hg19"; fi
    if [ -z "${WesLibrary}" ]; then WesLibrary="twist.exome.2.0"; fi
#---------------------------------------------------------------------------------------------------

#---| Manual Params |------------------------------------------------------------------------------#
    # BaseDir=/data/wes
    # BatchID=WES_25_09-WES_25_10_revised
    # TumorSeqID=WES_25_09_09   # WES_25_09_21
    # WesLibrary=twist.exome.2.0
    # Threads=16
    # GenomeAssembly=hg19
    # NormalType=TS
#--------------------------------------------------------------------------------------------------#

#---| FOLDERS |------------------------------------------------------------------------------------#
    DataDir=${BaseDir}/${BatchID}

    # tumor sample folders
    TumorBamDir=${DataDir}/${TumorSeqID}/bam
    TumorVcfDir=${DataDir}/${TumorSeqID}/vcf
    TumorQcDir=${DataDir}/${TumorSeqID}/qcfiles
   
    MutectDir_TONLY=${TumorVcfDir}/mutect2_tumor_only
    mkdir -p ${MutectDir_TONLY}
   
    LofreqDir_TONLY=${TumorVcfDir}/lofreq_tumor_only
    mkdir -p ${LofreqDir_TONLY}
   
    PurecnDir=${DataDir}/${TumorSeqID}/cnv/purecn/tumor_only
    mkdir -p ${PurecnDir}
   
    SomaticseqDir=${TumorVcfDir}/somaticseq_tumor_only
    mkdir -p ${SomaticseqDir}

    # additional folders
    TmpDir=${DataDir}/${TumorSeqID}/tmp
    LogDir=${DataDir}/${TumorSeqID}/log
    if [ ! -d ${LogDir} ]; then mkdir -p ${LogDir}; fi
#--------------------------------------------------------------------------------------------------#

#---| SWTOOLS RUN |--------------------------------------------------------------------------------#   
    GATK_VERSION="4.6.2.0"
    RunGatk4="singularity exec -B /storage,/data /storage/images/gatk-${GATK_VERSION}.sif gatk"
    RunLofreq="singularity exec -B /storage,/data /storage/images/lofreq-v2.1.sif lofreq"
    RunSomaticseq="singularity exec -B /storage,/data /storage/images/somaticseq-v3.11.1.sif somaticseq"
    RunPurecn="singularity exec -B /storage,/data /storage/images/purecn-v2.10.0.sif"
    COVERAGE_R=/opt/PureCN/Coverage.R
    PURECN_R=/opt/PureCN/PureCN.R
#--------------------------------------------------------------------------------------------------#

#---| RESOURCES |----------------------------------------------------------------------------------#
    
    # reference fasta
    GenomeFasta=/storage/references_and_index/${GenomeAssembly}/fasta/${GenomeAssembly}.primary_decoy_HLA_no.ALT.fa

    # resource vcf-files
    VCF_gnomad=/storage/references_and_index/${GenomeAssembly}/vcf/${GenomeAssembly}.gnomad.af.only.vcf.gz
    VCF_gnomad_mutect_annot=/storage/references_and_index/${GenomeAssembly}/vcf/${GenomeAssembly}.gnomad.split.gnomad-af.vcf.gz
    VCF_pon=/storage/references_and_index/${GenomeAssembly}/vcf/${GenomeAssembly}.panel.of.normals.vcf.gz
    VCF_common=/storage/references_and_index/${GenomeAssembly}/vcf/${GenomeAssembly}.common.biallelic.snp.vcf.gz 

    # target regions
    TargetRegionIntervals=/storage/references_and_index/${GenomeAssembly}/bed/${WesLibrary}/${GenomeAssembly}_${WesLibrary}.target.interval_list    
    TargetRegionBed=/storage/references_and_index/${GenomeAssembly}/bed/${WesLibrary}/${GenomeAssembly}_${WesLibrary}.target.bed 

    # purecn resources
    PurecnResourceDir=/storage/references_and_index/${GenomeAssembly}/purecn/resources_${WesLibrary}
    PurecnNormalDb=${PurecnResourceDir}/NORMALDB.RDS
    PurecnMappingBias=${PurecnResourceDir}/MAPPING-BIAS.RDS
    PurecnIntervals=${PurecnResourceDir}/INTERVALS.RDS
    PurecnSnpBlacklist=${PurecnResourceDir}/SNP-BLACKLIST.bed

#--------------------------------------------------------------------------------------------------#

#---| RunLog |------------------------------------------------------------------------------------#
    LogFile=${LogDir}/${TumorSeqID}.variant.calling.tumor.only.mode.$(date '+%Y%m%d').log
    echo " " > ${LogFile}
    echo "WES VARIANT CALLING Tumor-Only log" >> ${LogFile}
    echo "----------------------------------------------------------------------" >> ${LogFile}
    echo " RUN DATE             : $(date '+%Y-%m-%d %H:%M:%S')" >> ${LogFile}
    echo " BATCH ID             : ${BatchID}" >> ${LogFile}
    echo " SEQ ID               : ${TumorSeqID}" >> ${LogFile}
    echo " GENOME ASSEMBLY      : ${GenomeAssembly}" >> ${LogFile}
    echo " VARIANT CALLING MODE : Tumor-Only Mode" >> ${LogFile}
    echo "----------------------------------------------------------------------" >> ${LogFile}
    echo " " >> ${LogFile}
#==================================================================================================# 

#---| Output Prefix |------------------------------------------------------------------------------#
    # mutect2
    OutputPrefix_Mutect_VCF=${MutectDir_TONLY}/${TumorSeqID}.Tonly.mutect2
    OutputPrefix_Mutect_QC=${TumorQcDir}/${TumorSeqID}.Tonly.mutect2

    # lofreq
    OutputPrefix_Lofreq_VCF=${LofreqDir_TONLY}/${TumorSeqID}.Tonly.lofreq
    OutputPrefix_Lofreq_QC=${TumorQcDir}/${TumorSeqID}.Tonly.lofreq
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] VARIANT-CALLING MUTECT2 START." >> ${LogFile}
    echo " " >> ${LogFile}
#==================================================================================================# 

#---| MUTECT2 Tumor-Only Call |--------------------------------------------------------------------#

    ## run mutect2
    ${RunGatk4} --java-options "-Xmx32g" Mutect2 \
        -pairHMM AVX_LOGLESS_CACHING_OMP \
        --native-pair-hmm-threads ${Threads} \
        --reference ${GenomeFasta} \
        --intervals ${TargetRegionIntervals} \
        --input ${TumorBamDir}/${TumorSeqID}.analysisReady.bam \
        -tumor ${TumorSeqID} \
        --germline-resource ${VCF_gnomad} \
        --af-of-alleles-not-in-resource 2.5e-06 \
        --panel-of-normals ${VCF_pon} \
        --genotype-germline-sites true \
        --genotype-pon-sites true \
        --f1r2-tar-gz ${OutputPrefix_Mutect_QC}.f1r2.tar.gz \
        --f1r2-max-depth 600 \
        --min-base-quality-score 20 \
        --base-quality-score-threshold 20 \
        --interval-padding 50 \
        --max-mnp-distance 0 \
        --output ${OutputPrefix_Mutect_VCF}.unfiltered.vcf 2>&1 | tee -a ${LogFile}
    
    echo " " >> ${LogFile}

    ## orientation-bias models
    ${RunGatk4} --java-options "-Xmx8g" LearnReadOrientationModel \
        -I ${OutputPrefix_Mutect_QC}.f1r2.tar.gz \
        -O ${OutputPrefix_Mutect_QC}.read-orientation-model.tar.gz 2>&1 | tee -a ${LogFile}

    echo " " >> ${LogFile}

    ## pileup summary
    if [ ! -f "${OutputPrefix_Mutect_QC}.pileups.table" ]; then

        ${RunGatk4} --java-options "-Xmx4g" GetPileupSummaries \
            -R ${GenomeFasta} \
            -I ${TumorBamDir}/${TumorSeqID}.analysisReady.bam \
            -V ${VCF_common} \
            -L ${TargetRegionIntervals} \
            -O ${OutputPrefix_Mutect_QC}.pileups.table 2>&1 | tee -a ${LogFile}
        
        echo " " >> ${LogFile}
    fi
    
    ## calculate contamination
    ${RunGatk4} CalculateContamination \
        -I ${OutputPrefix_Mutect_QC}.pileups.table \
        -O ${OutputPrefix_Mutect_QC}.contamination.table \
        --tumor-segmentation ${OutputPrefix_Mutect_QC}.segments.table 2>&1 | tee -a ${LogFile}
    
    echo " " >> ${LogFile}

    ## filter mutect2 calls
    ${RunGatk4} FilterMutectCalls \
        -R ${GenomeFasta} \
        -V ${OutputPrefix_Mutect_VCF}.unfiltered.vcf \
        --contamination-table ${OutputPrefix_Mutect_QC}.contamination.table \
        --tumor-segmentation ${OutputPrefix_Mutect_QC}.segments.table \
        --ob-priors ${OutputPrefix_Mutect_QC}.read-orientation-model.tar.gz \
        -O ${OutputPrefix_Mutect_VCF}.filtered.vcf 2>&1 | tee -a ${LogFile}

    echo " " >> ${LogFile}

    # add gnomad info to vcf  
    bgzip -@ ${Threads} -c ${OutputPrefix_Mutect_VCF}.filtered.vcf > ${OutputPrefix_Mutect_VCF}.filtered.vcf.gz 
    tabix -p vcf ${OutputPrefix_Mutect_VCF}.filtered.vcf.gz

    bcftools norm -f ${GenomeFasta} -m -both -Oz -o ${OutputPrefix_Mutect_VCF}.filtered.split.vcf.gz ${OutputPrefix_Mutect_VCF}.filtered.vcf.gz
    tabix -p vcf ${OutputPrefix_Mutect_VCF}.filtered.split.vcf.gz

    bcftools annotate -a ${VCF_gnomad_mutect_annot} -c INFO/GNOMAD_AF -Oz ${OutputPrefix_Mutect_VCF}.filtered.split.vcf.gz > ${OutputPrefix_Mutect_VCF}.filtered.gnomad-af.vcf.gz
    tabix -p vcf ${OutputPrefix_Mutect_VCF}.filtered.gnomad-af.vcf.gz

    echo '##INFO=<ID=GNOMAD_DB,Number=0,Type=Flag,Description="Present in gnomAD population database">' > ${OutputPrefix_Mutect_VCF}.gnomad_db.header

    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/GNOMAD_AF\n' ${OutputPrefix_Mutect_VCF}.filtered.gnomad-af.vcf.gz \
        | awk 'BEGIN{OFS="\t"} $5 != "." && $5+0 >= 0.01 {print $1,$2,$3,$4,1}' \
        > ${OutputPrefix_Mutect_VCF}.gnomad_db.annot.tsv

    bgzip -c ${OutputPrefix_Mutect_VCF}.gnomad_db.annot.tsv > ${OutputPrefix_Mutect_VCF}.gnomad_db.annot.tsv.gz
    tabix -s 1 -b 2 -e 2 ${OutputPrefix_Mutect_VCF}.gnomad_db.annot.tsv.gz

    bcftools annotate \
        -h ${OutputPrefix_Mutect_VCF}.gnomad_db.header \
        -a ${OutputPrefix_Mutect_VCF}.gnomad_db.annot.tsv.gz \
        -c CHROM,POS,REF,ALT,INFO/GNOMAD_DB \
        -Oz -o ${OutputPrefix_Mutect_VCF}.filtered.gnomad-af-db.vcf.gz \
        ${OutputPrefix_Mutect_VCF}.filtered.gnomad-af.vcf.gz

    tabix -p vcf ${OutputPrefix_Mutect_VCF}.filtered.gnomad-af-db.vcf.gz

#---------------------------------------------------------------------------------------------------


#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | VARIANT-CALLING MUTECT2 FINISHED." >> ${LogFile}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | PURECN WITH MUTECT2 RESULT START." >> ${LogFile}
#==================================================================================================# 

#---| Mutect2-PureCN tumor-only mode |--------------------------------------------------------------

    export LC_ALL=C
    export LANG=C
    mkdir -p ${PurecnDir}/coverages

    ## coverages
    ${RunPurecn} Rscript ${COVERAGE_R} \
        --bam ${TumorBamDir}/${TumorSeqID}.analysisReady.bam \
        --intervals ${PurecnIntervals} \
        --out-dir ${PurecnDir}/coverages \
        --parallel --cores ${Threads} --force

    ## run purecn
    NumMutectVars=$(zgrep -vc '^#' ${OutputPrefix_Mutect_VCF}.filtered.gnomad-af.vcf.gz)

    if (( ${NumMutectVars} > 0 )); then

        RunPureCN="true"

        ${RunPurecn} Rscript ${PURECN_R} \
            --sampleid ${TumorSeqID} \
            --tumor ${PurecnDir}/coverages/${TumorSeqID}.analysisReady_coverage_loess.txt.gz \
            --vcf ${OutputPrefix_Mutect_VCF}.filtered.gnomad-af-db.vcf.gz \
            --mapping-bias-file ${PurecnMappingBias} \
            --normaldb ${PurecnNormalDb} \
            --intervals ${PurecnIntervals} \
            --snp-blacklist ${PurecnSnpBlacklist} \
            --genome ${GenomeAssembly} --min-base-quality 0 --fun-segmentation PSCBS \
            --popaf-info-field POPAF --db-info-flag GNOMAD_DB \
            --min-purity 0.05 --max-purity 0.95 --max-copy-number 8 --min-af 0.015 --log-ratio-calibration 0.1 \
            --post-optimize \
            --out ${PurecnDir}/${TumorSeqID}.Tonly.mutect2.purecn --out-vcf \
            --cores ${Threads} --seed 123456
    else
        RunPureCN="false"
    fi

#---------------------------------------------------------------------------------------------------

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | PURECN WITH MUTECT2 RESULT FINISHED." >> ${LogFile}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | VARIANT-CALLING LOFREQ START." >> ${LogFile}
#==================================================================================================# 

#---| LoFreq Tumor-Only Call |---------------------------------------------------------------------#
    
    ## indel qual bam
    ${RunLofreq} indelqual --dindel --ref ${GenomeFasta} \
        --out ${OutputPrefix_Lofreq_VCF}.indelqual.bam \
        ${TumorBamDir}/${TumorSeqID}.analysisReady.bam 

    samtools index -@ ${Threads} ${OutputPrefix_Lofreq_VCF}.indelqual.bam

    ## run lofreq (in parallel mode)
    ${RunLofreq} call-parallel \
        -f ${GenomeFasta} -l ${TargetRegionBed} \
        --call-indels --pp-threads ${Threads} \
        -o ${OutputPrefix_Lofreq_VCF}.vcf \
        ${OutputPrefix_Lofreq_VCF}.indelqual.bam

    bgzip -@ ${Threads} -c ${OutputPrefix_Lofreq_VCF}.vcf > ${OutputPrefix_Lofreq_VCF}.vcf.gz
    tabix -p vcf ${OutputPrefix_Lofreq_VCF}.vcf.gz

#---------------------------------------------------------------------------------------------------

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | VARIANT-CALLING LOFREQ FINISHED." >> ${LogFile}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | VARIANT-CALLING SOMATICSEQ START." >> ${LogFile}
#==================================================================================================# 

#---| SomaticSeq |----------------------------------------------------------------------------------

    ## run somaticseq
    ${RunSomaticseq} \
        --output-directory ${SomaticseqDir} \
        --genome-reference ${GenomeFasta} \
        --threads ${Threads} \
        --inclusion-region ${TargetRegionBed} \
        single \
        --bam-file ${TumorBamDir}/${TumorSeqID}.analysisReady.bam \
        --sample-name ${TumorSeqID} \
        --mutect2-vcf ${OutputPrefix_Mutect_VCF}.filtered.gnomad-af-db.vcf.gz \
        --lofreq-vcf ${OutputPrefix_Lofreq_VCF}.vcf.gz
        
    ${RunGatk4} UpdateVCFSequenceDictionary \
        -V ${SomaticseqDir}/Consensus.sSNV.vcf \
        --source-dictionary ${GenomeFasta}.dict \
        --output ${SomaticseqDir}/Consensus.sSNV.contig.fixed.vcf

    ${RunGatk4} UpdateVCFSequenceDictionary \
        -V ${SomaticseqDir}/Consensus.sINDEL.vcf \
        --source-dictionary ${GenomeFasta}.dict \
        --output ${SomaticseqDir}/Consensus.sINDEL.contig.fixed.vcf

    bcftools sort -Oz -o ${SomaticseqDir}/Consensus.sSNV.contig.fixed.vcf.gz ${SomaticseqDir}/Consensus.sSNV.contig.fixed.vcf
    tabix -p vcf ${SomaticseqDir}/Consensus.sSNV.contig.fixed.vcf.gz

    bcftools sort -Oz -o ${SomaticseqDir}/Consensus.sINDEL.contig.fixed.vcf.gz ${SomaticseqDir}/Consensus.sINDEL.contig.fixed.vcf
    tabix -p vcf ${SomaticseqDir}/Consensus.sINDEL.contig.fixed.vcf.gz
 
    bcftools concat -a -Oz -o ${SomaticseqDir}/${TumorSeqID}.Tonly.somaticseq.consensus.vcf.gz \
        ${SomaticseqDir}/Consensus.sSNV.contig.fixed.vcf.gz \
        ${SomaticseqDir}/Consensus.sINDEL.contig.fixed.vcf.gz
    tabix -p vcf ${SomaticseqDir}/${TumorSeqID}.Tonly.somaticseq.consensus.vcf.gz

    ## add purecn vcf
    if [ "${RunPureCN}" = "true" ]; then

        bcftools annotate -a ${PurecnDir}/${TumorSeqID}.Tonly.mutect2.purecn.vcf.gz \
            -c INFO/PureCN.PS,INFO/PureCN.ML.SOMATIC,INFO/PureCN.ML.C,INFO/PureCN.CS,INFO/PureCN.ML.LOH,INFO/PureCN.GHOMOZYGOUS,INFO/PureCN.LR,INFO/PureCN.CF,INFO/PureCN.ML.M \
            -Oz -o ${SomaticseqDir}/${TumorSeqID}.Tonly.somaticseq.consensus.purecn.vcf.gz \
            ${SomaticseqDir}/${TumorSeqID}.Tonly.somaticseq.consensus.vcf.gz
        tabix -p vcf ${SomaticseqDir}/${TumorSeqID}.Tonly.somaticseq.consensus.purecn.vcf.gz

        bgzip -dc ${SomaticseqDir}/${TumorSeqID}.Tonly.somaticseq.consensus.purecn.vcf.gz > ${TumorVcfDir}/${TumorSeqID}.Tonly.somaticseq.purecn.vcf
    else
        bgzip -dc ${SomaticseqDir}/${TumorSeqID}.Tonly.somaticseq.consensus.vcf.gz \
            > ${TumorVcfDir}/${TumorSeqID}.Tonly.somaticseq.consensus.no.purecn.vcf
    fi
    
#---------------------------------------------------------------------------------------------------

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | VARIANT-CALLING SOMATICSEQ FINISHED." >> ${LogFile}
#==================================================================================================# 

#==================================================================================================# 
    echo " " >> ${LogFile}
    echo "WES VARIANT CALLING Tumor-Only DONE." >> ${LogFile}
    echo "----------------------------------------------------------------------" >> ${LogFile}
    echo " " >> ${LogFile}
#==================================================================================================# 

