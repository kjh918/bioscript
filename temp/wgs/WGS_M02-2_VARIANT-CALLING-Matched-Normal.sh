
#!/bin/bash

# SCRIPT VERSION : v0.1
# GUIDELINE ID   : GMC-NGS-B02-A09
# DATE           : 2024-03-11
# AUTHOR         : kangsm@gencurix.com

#---| Input Options |------------------------------------------------------------------------------#
    usage() {
        echo "
        Usage: WGS_M02-2_MUTECT2-Matched-Normal.sh [ options ]
        
        --BASE_DIR          analysis base folder by NGS data type
                            default = /data/wes
        
        --BATCH_ID          [ REQUIRED ] dataset batch-id
        
        --TUMOR_SEQ_ID      [ REQURIED ] tumor sample-id
        
        --NORMAL_SEQ_ID     [ REQURIED ] matched-normal sample-id
       
        --GENOME_ASSEMBLY   reference genome assembly version
                            default = hg19

        --THREADS           parallel processing cpu cores 
                            default = 16

        --NORMAL_TYPE       type of matched-normal smaple. either 'TS' or 'ORG' 
                            TS = tissue or cells, ORG = organoid
                            default = TS

        -h | --help ]       Print this usage 
        "
    }
#--------------------------------------------------------------------------------------------------#

#---| Parse Option Arguments |---------------------------------------------------------------------#
    ARGS=$(getopt -a -o h: --long BATCH_ID:,TUMOR_SEQ_ID:,NORMAL_SEQ_ID:,BASE_DIR:,THREADS:,GENOME_ASSEMBLY:,NORMAL_TYPE:,help -- "$@" )
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
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  NormalSeqID=$2 ; shift 2 ;;
            esac  ;;
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
        --NORMAL_TYPE )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = TS" ; shift ;;
                *)  NormalType=$2 ; shift 2 ;;
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
    if [ -z "${NormalType}" ]; then NormalType="TS"; fi
#---------------------------------------------------------------------------------------------------

#---| Manual Params |------------------------------------------------------------------------------#
    BaseDir="/data/wgs/RUO_Service"
    BatchID="WGS-26-01"
    TumorSeqID="WGS-26-01_02"
    NormalSeqID="WGS-26-01_03"
    Threads=24
    GenomeAssembly="hg19"
    NormalType="TS"
#--------------------------------------------------------------------------------------------------#

#---| FOLDERS |------------------------------------------------------------------------------------#
    DataDir="${BaseDir}/${BatchID}"
    # tumor sample folders
    TumorBamDir="${DataDir}/${TumorSeqID}/bam"
    TumorVcfDir="${DataDir}/${TumorSeqID}/vcf"
    TumorQcDir="${DataDir}/${TumorSeqID}/qcfiles"
    # mutect2
    MutectDir_NT="${TumorVcfDir}/mutect2_matched_normal"
    mkdir -p "${MutectDir_NT}"
    # lofreq
    LofreqDir_NT="${TumorVcfDir}/lofreq_matched_normal"
    mkdir -p "${LofreqDir_NT}"
    # strelka2
    StrelkaDir="${TumorVcfDir}/strelka2_matched_normal"
    mkdir -p "${StrelkaDir}"
    # muse2
    Muse2Dir="${TumorVcfDir}/muse2_matched_normal"
    mkdir -p "${Muse2Dir}"
    # somaticseq
    SomaticseqDir="${TumorVcfDir}/somaticseq_matched_normal"
    mkdir -p "${SomaticseqDir}"
    # normal sample folders
    NormalBamDir="${DataDir}/${NormalSeqID}/bam"
    NormalQcDir="${DataDir}/${NormalSeqID}/qcfiles"   
    # additional folders
    TmpDir="${DataDir}/${TumorSeqID}/tmp"
    LogDir="${DataDir}/${TumorSeqID}/log"
    if [ ! -d "${LogDir}" ]; then mkdir -p "${LogDir}"; fi
#--------------------------------------------------------------------------------------------------#

#---| SWTOOLS RUN |--------------------------------------------------------------------------------#   
    GATK_VERSION="4.6.2.0"
    RunGatk4="singularity exec -B /storage,/data /storage/images/gatk-${GATK_VERSION}.sif gatk"
    RunLofreq="singularity exec -B /storage,/data /storage/images/lofreq-v2.1.sif lofreq"
    RunStrelkaSif="singularity exec -B /storage/,/data /storage/images/strelka2.sif"
    RunSomaticseq="singularity exec -B /storage,/data /storage/images/somaticseq-v3.11.1.sif somaticseq"
    RunMuse2="singularity exec -B /storage,/data /storage/images/muse2-v2.03.sif /MuSE/bin/MuSE"
#--------------------------------------------------------------------------------------------------#

#---| LogFile |------------------------------------------------------------------------------------#
    LogFile="${LogDir}/${TumorSeqID}.variant.calling.matched.normal.mode.$(date '+%Y%m%d').log"
    if [ ! -f "${LogFile}" ]; then touch "${LogFile}"; fi
#--------------------------------------------------------------------------------------------------#

#---| RESOURCES |----------------------------------------------------------------------------------#

    # reference fasta
    GenomeFasta="/storage/references_and_index/${GenomeAssembly}/fasta/${GenomeAssembly}.primary_decoy_HLA_no.ALT.fa"
    # resource vcf-files
    VCF_gnomad="/storage/references_and_index/${GenomeAssembly}/vcf/${GenomeAssembly}.gnomad.af.only.vcf.gz"
    VCF_gnomad_mutect_annot="/storage/references_and_index/${GenomeAssembly}/vcf/${GenomeAssembly}.gnomad.split.gnomad-af.vcf.gz"
    VCF_pon="/storage/references_and_index/${GenomeAssembly}/vcf/${GenomeAssembly}.panel.of.normals.vcf.gz"
    VCF_common="/storage/references_and_index/${GenomeAssembly}/vcf/${GenomeAssembly}.common.biallelic.snp.filtered.for.wgs.vcf.gz"
    VCF_dbsnp="/storage/references_and_index/${GenomeAssembly}/vcf/${GenomeAssembly}.dbsnp138.vcf.gz"
    # call regions (autosome + X , Y only & blacklist filtered)
    WgsCallRegionBed="/storage/references_and_index/${GenomeAssembly}/bed/${GenomeAssembly}.wgs.call.target.region.bed"
    WgsCallRegionIntervals="/storage/references_and_index/${GenomeAssembly}/bed/${GenomeAssembly}.wgs.call.target.region.interval_list"

#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo " " >> ${LogFile}
    echo "WES VARIANT CALLING Matched-Normal log" >> ${LogFile}
    echo "----------------------------------------------------------------------" >> ${LogFile}
    echo " RUN DATE             : $(date '+%Y-%m-%d %H:%M:%S')" >> ${LogFile}
    echo " BATCH ID             : ${BatchID}" >> ${LogFile}
    echo " TUMOR SEQ ID         : ${TumorSeqID}" >> ${LogFile}
    echo " NORMAL SEQ ID        : ${NormalSeqID}" >> ${LogFile}
    echo " GENOME ASSEMBLY      : ${GenomeAssembly}" >> ${LogFile}
    echo " VARIANT CALLING MODE : Matched Normal Mode" >> ${LogFile}
    echo "----------------------------------------------------------------------" >> ${LogFile}
    echo " " >> ${LogFile}
#==================================================================================================# 

#---| Output Prefix |------------------------------------------------------------------------------#
    
    # normal type tag
    if [ "${NormalType}" = "ORG" ]; then
        PrefixTag="NT.ORG"
    else
        PrefixTag="NT"
    fi
    # mutect2
    OutputPrefix_Mutect_VCF="${MutectDir_NT}/${TumorSeqID}.${PrefixTag}.mutect2"
    OutputPrefix_Mutect_QC="${TumorQcDir}/${TumorSeqID}.${PrefixTag}.mutect2"
    # muse2 
    OutputPrefix_Muse2_TXT="${Muse2Dir}/${TumorSeqID}.${PrefixTag}"
    OutputPrefix_Muse2_VCF="${Muse2Dir}/${TumorSeqID}.${PrefixTag}.muse2"
    # lofreq
    OutputPrefix_Lofreq_VCF="${LofreqDir_NT}/${TumorSeqID}.${PrefixTag}.lofreq"
    OutputPrefix_Lofreq_QC="${TumorQcDir}/${TumorSeqID}.${PrefixTag}.lofreq"

#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | VARIANT-CALLING MUTECT2 START." >> ${LogFile}
#==================================================================================================# 

#---| MUTECT2 Matched-Normal Call |----------------------------------------------------------------#
    
    ## run mutect2
    ${RunGatk4} --java-options "-Xmx32g" Mutect2 \
        -pairHMM "AVX_LOGLESS_CACHING_OMP" \
        --native-pair-hmm-threads "${Threads}" \
        --reference "${GenomeFasta}" \
        --intervals "${WgsCallRegionIntervals}" \
        --input "${TumorBamDir}/${TumorSeqID}.analysisReady.bam" \
        --input "${NormalBamDir}/${NormalSeqID}.analysisReady.bam" \
        -tumor "${TumorSeqID}" -normal "${NormalSeqID}" \
        --germline-resource "${VCF_gnomad}" \
        --af-of-alleles-not-in-resource "2.5e-06" \
        --panel-of-normals "${VCF_pon}" \
        --f1r2-tar-gz "${OutputPrefix_Mutect_QC}.f1r2.tar.gz" \
        --f1r2-max-depth 300 \
        --min-base-quality-score 20 \
        --base-quality-score-threshold 20 \
        --max-mnp-distance 0 \
        --output "${OutputPrefix_Mutect_VCF}.unfiltered.vcf"
    
    ## orientation-bias model
    ${RunGatk4} --java-options "-Xmx8g" LearnReadOrientationModel \
        -I "${OutputPrefix_Mutect_QC}.f1r2.tar.gz" \
        -O "${OutputPrefix_Mutect_QC}.read-orientation-model.tar.gz"

    ## pileup summary - tumor
    ${RunGatk4} --java-options "-Xmx4g" GetPileupSummaries \
        -R "${GenomeFasta}" \
        -I "${TumorBamDir}/${TumorSeqID}.analysisReady.bam" \
        -V "${VCF_common}" \
        -L "${VCF_common}" \
        -O "${OutputPrefix_Mutect_QC}.pileups.table"


    ## pileup summary - normal
    ${RunGatk4} --java-options "-Xmx4g" GetPileupSummaries \
        -R "${GenomeFasta}" \
        -I "${NormalBamDir}/${NormalSeqID}.analysisReady.bam" \
        -V "${VCF_common}" \
        -L "${VCF_common}" \
        -O "${NormalQcDir}.${NormalSeqID}.${PrefixTag}.mutect2"


    ## calculate contamination
    ${RunGatk4} CalculateContamination \
        -I "${OutputPrefix_Mutect_QC}.pileups.table" \
        -matched "${NormalQcDir}.${NormalSeqID}.${PrefixTag}.mutect2" \
        -O "${OutputPrefix_Mutect_QC}.contamination.table" \
        --tumor-segmentation "${OutputPrefix_Mutect_QC}.segments.table"

    ## filter mutect2 calls
    ${RunGatk4} FilterMutectCalls \
        -R "${GenomeFasta}" \
        -V "${OutputPrefix_Mutect_VCF}.unfiltered.vcf" \
        --contamination-table "${OutputPrefix_Mutect_QC}.contamination.table" \
        --tumor-segmentation "${OutputPrefix_Mutect_QC}.segments.table" \
        --ob-priors "${OutputPrefix_Mutect_QC}.read-orientation-model.tar.gz" \
        -O "${OutputPrefix_Mutect_VCF}.filtered.vcf"
    
    # add gnomad info to vcf   
    bcftools norm -m -any -f "${GenomeFasta}" -Oz -o "${OutputPrefix_Mutect_VCF}.filtered.norm.vcf.gz" \
        "${OutputPrefix_Mutect_VCF}.filtered.vcf"
    tabix -p vcf "${OutputPrefix_Mutect_VCF}.filtered.norm.vcf.gz"

    SomaticSeqInput_MUTECT2="${OutputPrefix_Mutect_VCF}.filtered.norm.vcf.gz"

#---------------------------------------------------------------------------------------------------

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | VARIANT-CALLING MUTECT2 FINISHED." >> ${LogFile}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | VARIANT-CALLING LOFREQ START." >> ${LogFile}
#==================================================================================================# 

#---| LoFreq Matched-Normal Call |-----------------------------------------------------------------#
    
    ## run lofreq
    lofreq_log=${LofreqDir_NT}/run.log
    touch ${lofreq_log}
    
    if [[ "${GenomeAssembly}" = "hg19" ]]; then
        CHR_LIST=($(printf "%s " {1..22} X Y))
    else
        CHR_LIST=($(printf "chr%s " {1..22} X Y))
    fi

    LofreqTmpDir="${LofreqDir_NT}/chr_tmp"
    mkdir -p "${LofreqTmpDir}"
    LofreqJobThreads=$(( Threads / 2 ))

    # create indelqual bam
    parallel -j ${LofreqJobThreads} -k "
        samtools view -bh -@ 2 ${TumorBamDir}/${TumorSeqID}.analysisReady.bam {} > ${LofreqTmpDir}/tumor.{}.bam && 
        samtools view -bh -@ 2 ${NormalBamDir}/${NormalSeqID}.analysisReady.bam {} > ${LofreqTmpDir}/normal.{}.bam &&
        samtools index -@ 2 ${LofreqTmpDir}/tumor.{}.bam &&
        samtools index -@ 2 ${LofreqTmpDir}/normal.{}.bam &&
        ${RunLofreq} indelqual --dindel -f ${GenomeFasta} ${LofreqTmpDir}/tumor.{}.bam -o ${LofreqTmpDir}/tumor.{}.indelqual.bam &&
        ${RunLofreq} indelqual --dindel -f ${GenomeFasta} ${LofreqTmpDir}/normal.{}.bam -o ${LofreqTmpDir}/normal.{}.indelqual.bam &&
        samtools index -@ 2 ${LofreqTmpDir}/tumor.{}.indelqual.bam &&
        samtools index -@ 2 ${LofreqTmpDir}/normal.{}.indelqual.bam &&
        rm ${LofreqTmpDir}/tumor.{}.bam  ${LofreqTmpDir}/tumor.{}.bam.bai ${LofreqTmpDir}/normal.{}.bam ${LofreqTmpDir}/normal.{}.bam.bai
    " ::: ${CHR_LIST[@]}
    
    # call somatic variants each chromosome
    for CHROM in ${CHR_LIST[@]}
    do
        awk -v chr="${CHROM}" '$1==chr' "${WgsCallRegionBed}" > "${LofreqTmpDir}/target.chrom.${CHROM}.bed"

        ${RunLofreq} somatic \
            -n ${LofreqTmpDir}/normal.${CHROM}.indelqual.bam \
            -t ${LofreqTmpDir}/tumor.${CHROM}.indelqual.bam \
            -f ${GenomeFasta} -l ${LofreqTmpDir}/target.chrom.${CHROM}.bed \
            --call-indels --threads ${Threads} \
            -o ${LofreqTmpDir}/${TumorSeqID}.target.chrom.${CHROM}.

        cp ${LofreqTmpDir}/${TumorSeqID}.target.chrom.${CHROM}.somatic_final.indels.vcf.gz  ${LofreqTmpDir}/final_vcfs/
        cp ${LofreqTmpDir}/${TumorSeqID}.target.chrom.${CHROM}.somatic_final.indels.vcf.gz.tbi ${LofreqTmpDir}/final_vcfs/
        cp ${LofreqTmpDir}/${TumorSeqID}.target.chrom.${CHROM}.somatic_final.snvs.vcf.gz  ${LofreqTmpDir}/final_vcfs/
        cp ${LofreqTmpDir}/${TumorSeqID}.target.chrom.${CHROM}.somatic_final.snvs.vcf.gz.tbi  ${LofreqTmpDir}/final_vcfs/

        #rm ${LofreqTmpDir}/${TumorSeqID}.target.chrom.${CHROM}.*
        #rm ${LofreqTmpDir}/tumor.${CHROM}.*
        #rm ${LofreqTmpDir}/normal.${CHROM}.*
        #rm ${LofreqTmpDir}/target.chrom.${CHROM}.bed
    done


    # merge final vcfs
    bcftools concat -a $(readlink -f ${LofreqTmpDir}/final_vcfs/* | grep snvs.vcf.gz$ | sort -V) \
        -Oz -o "${OutputPrefix_Lofreq_VCF}.somatic_final.snvs.vcf.gz" 
    tabix -p vcf "${OutputPrefix_Lofreq_VCF}.somatic_final.snvs.vcf.gz"

    bcftools concat -a $(readlink -f ${LofreqTmpDir}/final_vcfs/* | grep indels.vcf.gz$ | sort -V) \
        -Oz -o "${OutputPrefix_Lofreq_VCF}.somatic_final.indels.vcf.gz" 
    tabix -p vcf "${OutputPrefix_Lofreq_VCF}.somatic_final.indels.vcf.gz"

    bcftools norm -m -any -f "${GenomeFasta}" -Oz -o "${OutputPrefix_Lofreq_VCF}.somatic_final.norm.snvs.vcf.gz" \
        "${OutputPrefix_Lofreq_VCF}.somatic_final.snvs.vcf.gz"
    tabix -p vcf "${OutputPrefix_Lofreq_VCF}.somatic_final.norm.snvs.vcf.gz"

    bcftools norm -m -any -f "${GenomeFasta}" -Oz -o "${OutputPrefix_Lofreq_VCF}.somatic_final.norm.indels.vcf.gz" \
        "${OutputPrefix_Lofreq_VCF}.somatic_final.indels.vcf.gz"
    tabix -p vcf "${OutputPrefix_Lofreq_VCF}.somatic_final.norm.indels.vcf.gz"

    SomaticSeqInput_LOFREQ_SNV="${OutputPrefix_Lofreq_VCF}.somatic_final.norm.snvs.vcf.gz"
    SomaticSeqInput_LOFREQ_INDEL="${OutputPrefix_Lofreq_VCF}.somatic_final.norm.indels.vcf.gz"

#---------------------------------------------------------------------------------------------------

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | VARIANT-CALLING LOFREQ FINISHED." >> ${LogFile}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | VARIANT-CALLING MUSE2 START." >> ${LogFile}
#==================================================================================================# 

#---| MuSE2 Matched-Normal Mode |-----------------------------------------------------------------

    ${RunMuse2} call \
        -f "${GenomeFasta}" -n "${Threads}" \
        -O "${OutputPrefix_Muse2_TXT}" \
        "${TumorBamDir}/${TumorSeqID}.analysisReady.bam" "${NormalBamDir}/${NormalSeqID}.analysisReady.bam"

    ${RunMuse2} sump -G \
        -I "${OutputPrefix_Muse2_TXT}.MuSE.txt" \
        -O "${OutputPrefix_Muse2_VCF}.vcf" \
        -n "${Threads}" \
        -D "${VCF_dbsnp}"

    bgzip -@ "${Threads}" -c "${OutputPrefix_Muse2_VCF}.vcf" > "${OutputPrefix_Muse2_VCF}.vcf.gz"
    tabix -p vcf "${OutputPrefix_Muse2_VCF}.vcf.gz"

    Muse2CalledVarNum=$(zcat ${OutputPrefix_Muse2_VCF}.vcf.gz | grep -v ^"#" | wc -l)
    
    if (( ${Muse2CalledVarNum} > 0 )); then
        bcftools view -i 'FILTER="PASS" || FILTER="Tier1" || FILTER="Tier2" || FILTER="Tier3" || FILTER="Tier4"' "${OutputPrefix_Muse2_VCF}.vcf.gz" \
            | bcftools norm -m -any -f "${GenomeFasta}" -Oz -o "${OutputPrefix_Muse2_VCF}.filtered.norm.vcf.gz"
        tabix -p vcf "${OutputPrefix_Muse2_VCF}.filtered.norm.vcf.gz"

        SomaticSeqInput_MUSE2="${OutputPrefix_Muse2_VCF}.filtered.norm.vcf.gz"
    fi

#---------------------------------------------------------------------------------------------------

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | VARIANT-CALLING MUSE2 FINISHED." >> ${LogFile}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | VARIANT-CALLING STRELKA2 START." >> ${LogFile}
#==================================================================================================# 

#---| Strelka2 Matched-Noormal Mode |---------------------------------------------------------------

    ## strelka2 config
    ${RunStrelkaSif} python /tools/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py \
        --tumorBam "${TumorBamDir}/${TumorSeqID}.analysisReady.bam" \
        --normalBam "${NormalBamDir}/${NormalSeqID}.analysisReady.bam" \
        --referenceFasta "${GenomeFasta}" \
        --callRegions "${WgsCallRegionBed}.gz" \
        --runDir "${StrelkaDir}"

    # run strelka2
    ${RunStrelkaSif} python ${StrelkaDir}/runWorkflow.py -m local -j ${Threads} -g 64

    bcftools norm -m -any -f "${GenomeFasta}" -Oz -o "${StrelkaDir}/results/variants/somatic.snvs.norm.vcf.gz" \
        "${StrelkaDir}/results/variants/somatic.snvs.vcf.gz"
    tabix -p vcf "${StrelkaDir}/results/variants/somatic.snvs.norm.vcf.gz"

    bcftools norm -m -any -f "${GenomeFasta}" -Oz -o "${StrelkaDir}/results/variants/somatic.indels.norm.vcf.gz" \
        "${StrelkaDir}/results/variants/somatic.indels.vcf.gz" 
    tabix -p vcf "${StrelkaDir}/results/variants/somatic.indels.norm.vcf.gz"

    SomaticSeqInput_STRELKA2_SNV="${StrelkaDir}/results/variants/somatic.snvs.norm.vcf.gz"
    SomaticSeqInput_STRELKA2_INDEL="${StrelkaDir}/results/variants/somatic.indels.norm.vcf.gz"

#---------------------------------------------------------------------------------------------------

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | VARIANT-CALLING STRELKA2 FINISHED." >> ${LogFile}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | VARIANT-CALLING SOMATICSEQ START." >> ${LogFile}
#==================================================================================================# 

#---| Run SomaticSeq |------------------------------------------------------------------------------

    ## VCF check and create args
        # mutect2 
        InputArgs_VCF=("--mutect2-vcf ${SomaticSeqInput_MUTECT2}")

        # strelka2 
        strelka_snv_vcf="${SomaticSeqInput_STRELKA2_SNV}"
        strelka_called_snv=$(zcat "${strelka_snv_vcf}" | grep -v ^# | wc -l)
        if [[ -f "${strelka_snv_vcf}" ]]; then
            if (( ${strelka_called_snv} > 0 )); then 
                InputArgs_VCF+=("--strelka-snv ${SomaticSeqInput_STRELKA2_SNV}")
            fi
        fi

        strelka_indel_vcf="${SomaticSeqInput_STRELKA2_INDEL}"
        strelka_called_indel=$(zcat "${strelka_indel_vcf}" | grep -v ^# | wc -l)
        if [[ -f "${strelka_indel_vcf}" ]]; then
            if (( ${strelka_called_snv} > 0 )); then 
                InputArgs_VCF+=("--strelka-indel ${SomaticSeqInput_STRELKA2_INDEL}")
            fi
        fi

        # muse2
        muse2_vcf="${SomaticSeqInput_MUSE2}"
        muse2_called_vars=$(zcat "${muse2_vcf}" | grep -v ^# | wc -l)
        if [[ -f "${muse2_vcf}" ]]; then
            if (( ${muse2_called_vars} > 0 )); then 
                InputArgs_VCF+=("--muse-vcf ${SomaticSeqInput_MUSE2}")
            fi
        fi

        # lofreq 
        lofreq_snv_vcf="${SomaticSeqInput_LOFREQ_SNV}"
        lofreq_called_snv=$(zcat "${lofreq_snv_vcf}" | grep -v ^# | wc -l)
        if [[ -f "${lofreq_snv_vcf}" ]]; then
            if (( ${lofreq_called_snv} > 0 )); then 
                InputArgs_VCF+=("--lofreq-snv ${SomaticSeqInput_LOFREQ_SNV}")
            fi
        fi

        lofreq_indel_vcf="${SomaticSeqInput_LOFREQ_INDEL}"
        lofreq_called_indel=$(zcat "${lofreq_indel_vcf}" | grep -v ^# | wc -l)
        if [[ -f "${lofreq_indel_vcf}" ]]; then
            if (( ${lofreq_called_snv} > 0 )); then 
                InputArgs_VCF+=("--lofreq-indel ${SomaticSeqInput_LOFREQ_INDEL}")
            fi
        fi
    #
    
    ### RUN SOMATICSEQ

    ${RunSomaticseq} \
        --output-directory "${SomaticseqDir}" \
        --genome-reference "${GenomeFasta}" \
        --threads "${Threads}" \
        --inclusion-region "${WgsCallRegionBed}" \
        paired \
        --tumor-bam-file "${TumorBamDir}/${TumorSeqID}.analysisReady.bam" \
        --normal-bam-file "${NormalBamDir}/${NormalSeqID}.analysisReady.bam" \
        --tumor-sample "${TumorSeqID}" \
        --normal-sample "${NormalSeqID}" \
        ${InputArgs_VCF[@]}
        
        
    ## add full-contigs into vcf
    ${RunGatk4} UpdateVCFSequenceDictionary \
        -V "${SomaticseqDir}/Consensus.sSNV.vcf" \
        --source-dictionary "${GenomeFasta}.dict" \
        --output "${SomaticseqDir}/Consensus.sSNV.contig.fix.vcf"

    ${RunGatk4} UpdateVCFSequenceDictionary \
        -V "${SomaticseqDir}/Consensus.sINDEL.vcf" \
        --source-dictionary "${GenomeFasta}.dict" \
        --output "${SomaticseqDir}/Consensus.sINDEL.contig.fix.vcf"

    ## sort vcf 
    bcftools sort -Oz -o "${SomaticseqDir}/Consensus.sSNV.contig.fix.vcf.gz" "${SomaticseqDir}/Consensus.sSNV.contig.fix.vcf"
    tabix -p vcf "${SomaticseqDir}/Consensus.sSNV.contig.fix.vcf.gz"

    bcftools sort -Oz -o "${SomaticseqDir}/Consensus.sINDEL.contig.fix.vcf.gz" "${SomaticseqDir}/Consensus.sINDEL.contig.fix.vcf"
    tabix -p vcf "${SomaticseqDir}/Consensus.sINDEL.contig.fix.vcf.gz"

    ## merge vcf
    bcftools concat -a -Oz -o "${SomaticseqDir}/${TumorSeqID}.somaticseq.consensus.vcf.gz" \
        "${SomaticseqDir}/Consensus.sSNV.contig.fix.vcf.gz" \
        "${SomaticseqDir}/Consensus.sINDEL.contig.fix.vcf.gz"
    tabix -p vcf "${SomaticseqDir}/${TumorSeqID}.somaticseq.consensus.vcf.gz"

    ## FINAL SOMATICSEQ VCF
    bcftools view -f PASS "${SomaticseqDir}/${TumorSeqID}.somaticseq.consensus.vcf.gz" \
        | bcftools norm -m -any -f "${GenomeFasta}" -Oz -o "${SomaticseqDir}/${TumorSeqID}.somaticseq.consensus.pass.vcf.gz"

    cp "${SomaticseqDir}/${TumorSeqID}.somaticseq.consensus.pass.vcf.gz" "${TumorVcfDir}/${TumorSeqID}.${PrefixTag}.somaticseq.consensus.pass.vcf.gz"
    tabix -p vcf "${TumorVcfDir}/${TumorSeqID}.${PrefixTag}.somaticseq.consensus.pass.vcf.gz"

    ln -Tsf "${TumorVcfDir}/${TumorSeqID}.${PrefixTag}.somaticseq.consensus.pass.vcf.gz" "${TumorVcfDir}/${TumorSeqID}.${PrefixTag}.VARIANTS.vcf.gz"
    ln -Tsf "${TumorVcfDir}/${TumorSeqID}.${PrefixTag}.somaticseq.consensus.pass.vcf.gz.tbi" "${TumorVcfDir}/${TumorSeqID}.${PrefixTag}.VARIANTS.vcf.gz.tbi"
    
###-------------------------------------------------------------------------------------------------


#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | VARIANT-CALLING SOMATICSEQ FINISHED." >> ${LogFile}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | VEP-ANNATATION-READY VCF CREATED." >> ${LogFile}
#==================================================================================================# 

#==================================================================================================# 
    echo " " >> ${LogFile}
    echo "WES VARIANT CALLING Matched-Normal DONE." >> ${LogFile}
    echo "----------------------------------------------------------------------" >> ${LogFile}
    echo " " >> ${LogFile}
#==================================================================================================# 

