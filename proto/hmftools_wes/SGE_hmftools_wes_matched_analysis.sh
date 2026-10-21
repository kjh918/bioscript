

#! /bin/bash

#$ -N WES_VEP
#$ -S /bin/bash
#$ -pe mpi 12
#$ -t 1-18
#$ -tc 6
#$ -e /storage/home/kangsm/analysis/gcx.pipeline.use.hmftools.in.wes/01.preprocess/runlog/wes/fix.$TASK_ID.err
#$ -o /storage/home/kangsm/analysis/gcx.pipeline.use.hmftools.in.wes/01.preprocess/runlog/wes/fix.$TASK_ID.out

#---/ Sample List File
    SampleListFile="/storage/home/kangsm/analysis/gcx.pipeline.use.hmftools.in.wes/01.preprocess/analysis.list"
#-------------------------------------------------------------------------------


#---/ Sample-Info from Sample-List-File |---------------------------------------
    # SeqID RGID BatchID WesLibrary
    RunInfo=$(sed -n -e "${SGE_TASK_ID} p" ${SampleListFile}) 
    BatchID=$(echo ${RunInfo} | awk '{print $3}')
    TumorSeqID=$(echo ${RunInfo} | awk '{print $1}')
    NormalSeqID=$(echo ${RunInfo} | awk '{print $2}')
    WesLibrary=$(echo ${RunInfo} | awk '{print $4}')

#-------------------------------------------------------------------------------

#---/ Manual Params |-----------------------------------------------------------
    BaseDir="/data/wes"
    #BatchID="WES-panel.of.normals-hg19-sureselect.v6"
    #SampleListFile="/data/pon/WES-panel.of.normals-hg19-sureselect.v6/meta/sample.list"
    Threads=12
    GenomeAssembly="hg19"
    #WesLibrary="sureselect.v6"
#-------------------------------------------------------------------------------

#---| SWTOOLS RUN |--------------------------------------------------------------------------------#   
    GATK_VERSION="4.6.2.0"
    RunGatk4="singularity exec -B /storage,/data /storage/images/gatk-${GATK_VERSION}.sif gatk"
    RunLofreq="singularity exec -B /storage,/data /storage/images/lofreq-v2.1.sif lofreq"
    RunStrelkaSif="singularity exec -B /storage/,/data /storage/images/strelka2.sif"
    RunMuse2="singularity exec -B /storage,/data /storage/images/muse2-v2.03.sif /MuSE/bin/MuSE"
    RunVardict="java -Xmx16g -jar /storage/apps/vardict-v1.8.3/lib/VarDict-1.8.3.jar"
    RunVardictSomaticTest="/storage/apps/vardict-v1.8.3/bin/testsomatic.R"
    RunVardictCreateVcf="/storage/apps/vardict-v1.8.3/bin/var2vcf_paired.pl"

    RunSomaticseq="singularity exec -B /storage,/data /storage/images/somaticseq-v3.11.1.sif somaticseq"
    RunPurecn="singularity exec -B /storage,/data /storage/images/purecn-v2.10.0.sif"
    COVERAGE_R=/opt/PureCN/Coverage.R
    PURECN_R=/opt/PureCN/PureCN.R
#--------------------------------------------------------------------------------------------------#

#---| RESOURCES |----------------------------------------------------------------------------------#

    # reference fasta
    GenomeFasta="/storage/references_and_index/${GenomeAssembly}/fasta/${GenomeAssembly}.primary_decoy_HLA_no.ALT.fa"
    # resource vcf-files
    VCF_gnomad="/storage/references_and_index/${GenomeAssembly}/vcf/${GenomeAssembly}.gnomad.af.only.vcf.gz"
    VCF_gnomad_mutect_annot="/storage/references_and_index/${GenomeAssembly}/vcf/${GenomeAssembly}.gnomad.split.gnomad-af.vcf.gz"
    VCF_pon="/storage/references_and_index/${GenomeAssembly}/vcf/${GenomeAssembly}.panel.of.normals.vcf.gz"
    VCF_common="/storage/references_and_index/${GenomeAssembly}/vcf/${GenomeAssembly}.common.biallelic.snp.vcf.gz"
    VCF_dbsnp="/storage/references_and_index/${GenomeAssembly}/vcf/${GenomeAssembly}.dbsnp138.vcf.gz" 
    # target regions
    TargetRegionIntervals="/storage/references_and_index/${GenomeAssembly}/bed/${WesLibrary}/${GenomeAssembly}_${WesLibrary}.target.interval_list"    
    TargetRegionBed="/storage/references_and_index/${GenomeAssembly}/bed/${WesLibrary}/${GenomeAssembly}_${WesLibrary}.target.bed" 
    # purecn resources
    PurecnResourceDir="/storage/references_and_index/${GenomeAssembly}/purecn/resources_${WesLibrary}"
    PurecnNormalDb="${PurecnResourceDir}/NORMALDB.RDS"
    PurecnMappingBias="${PurecnResourceDir}/MAPPING-BIAS.RDS"
    PurecnIntervals="${PurecnResourceDir}/INTERVALS.RDS"
    PurecnSnpBlacklist="${PurecnResourceDir}/SNP-BLACKLIST.bed"

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
    StrelkaDir_NT="${TumorVcfDir}/strelka2_matched_normal"
    mkdir -p "${StrelkaDir_NT}"
    # muse2
    Muse2Dir_NT="${TumorVcfDir}/muse2_matched_normal"
    mkdir -p "${Muse2Dir_NT}"
    # vardict
    VardictDir_NT="${TumorVcfDir}/vardict_matched_normal"
    mkdir -p "${VardictDir_NT}"
    # somaticseq
    SomaticseqDir="${TumorVcfDir}/somaticseq_matched_normal"
    mkdir -p "${SomaticseqDir}"
    # purecn input vcf
    PurecnVcfPrepDir="${TumorVcfDir}/purecn_input_matched_normal"
    mkdir -p "${PurecnVcfPrepDir}"
    # purecn result
    PurecnDir="${DataDir}/${TumorSeqID}/cnv/purecn/matched_normal"
    mkdir -p "${PurecnDir}"
    # normal sample folders
    NormalBamDir="${DataDir}/${NormalSeqID}/bam"
    NormalQcDir="${DataDir}/${NormalSeqID}/qcfiles" 
    # additional folders
    TmpDir="${DataDir}/${TumorSeqID}/tmp"
#--------------------------------------------------------------------------------------------------#

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
    # lofreq
    OutputPrefix_Lofreq_VCF="${LofreqDir_NT}/${TumorSeqID}.${PrefixTag}.lofreq"
    OutputPrefix_Lofreq_QC="${TumorQcDir}/${TumorSeqID}.${PrefixTag}.lofreq"
    # muse2 
    OutputPrefix_Muse2_TXT="${Muse2Dir_NT}/${TumorSeqID}.${PrefixTag}"
    OutputPrefix_Muse2_VCF="${Muse2Dir_NT}/${TumorSeqID}.${PrefixTag}.muse2"
    # vardict
    OutputPrefix_Vardict="${VardictDir_NT}/${TumorSeqID}.${PrefixTag}.vardict"

#--------------------------------------------------------------------------------------------------#

    SomaticSeqInput_MUTECT2="${OutputPrefix_Mutect_VCF}.filtered.norm.vcf.gz"
    SomaticSeqInput_LOFREQ_SNV="${OutputPrefix_Lofreq_VCF}.somatic_final.norm.snvs.vcf.gz"
    SomaticSeqInput_LOFREQ_INDEL="${OutputPrefix_Lofreq_VCF}.somatic_final.norm.indels.vcf.gz"

    bcftools norm -m -any -f "${GenomeFasta}" -Oz -o "${StrelkaDir_NT}/results/variants/somatic.snvs.norm.vcf.gz" \
        "${StrelkaDir_NT}/results/variants/somatic.snvs.vcf.gz"
    tabix -p vcf "${StrelkaDir_NT}/results/variants/somatic.snvs.norm.vcf.gz"

    bcftools norm -m -any -f "${GenomeFasta}" -Oz -o "${StrelkaDir_NT}/results/variants/somatic.indels.norm.vcf.gz" \
        "${StrelkaDir_NT}/results/variants/somatic.indels.vcf.gz"
    tabix -p vcf "${StrelkaDir_NT}/results/variants/somatic.indels.norm.vcf.gz"

    SomaticSeqInput_STRELKA2_SNV="${StrelkaDir_NT}/results/variants/somatic.snvs.norm.vcf.gz"
    SomaticSeqInput_STRELKA2_INDEL="${StrelkaDir_NT}/results/variants/somatic.indels.norm.vcf.gz"

    SomaticSeqInput_MUSE2="${OutputPrefix_Muse2_VCF}.filtered.norm.vcf.gz"
    SomaticSeqInput_VARDICT="${OutputPrefix_Vardict}.somatic.norm.vcf.gz"

#---| Run SomaticSeq |------------------------------------------------------------------------------

    # INPUT VCF ARGS
        # mutect2 
        InputArgs_VCF=("--mutect2-vcf ${SomaticSeqInput_MUTECT2}")

        # lofreq 
        lofreq_snv_vcf=${SomaticSeqInput_LOFREQ_SNV}
        lofreq_called_snv=$(zcat ${lofreq_snv_vcf} | grep -v ^# | wc -l)
        if [[ -f "${lofreq_snv_vcf}" ]]; then
            if (( ${lofreq_called_snv} > 0 )); then 
                InputArgs_VCF+=("--lofreq-snv ${SomaticSeqInput_LOFREQ_SNV}")
            fi
        fi

        lofreq_indel_vcf=${SomaticSeqInput_LOFREQ_INDEL}
        lofreq_called_indel=$(zcat ${lofreq_indel_vcf} | grep -v ^# | wc -l)
        if [[ -f "${lofreq_indel_vcf}" ]]; then
            if (( ${lofreq_called_snv} > 0 )); then 
                InputArgs_VCF+=("--lofreq-indel ${SomaticSeqInput_LOFREQ_INDEL}")
            fi
        fi

        # strelka2 
        strelka_snv_vcf=${SomaticSeqInput_STRELKA2_SNV}
        strelka_called_snv=$(zcat ${strelka_snv_vcf} | grep -v ^# | wc -l)
        if [[ -f "${strelka_snv_vcf}" ]]; then
            if (( ${strelka_called_snv} > 0 )); then 
                InputArgs_VCF+=("--strelka-snv ${SomaticSeqInput_STRELKA2_SNV}")
            fi
        fi

        strelka_indel_vcf=${SomaticSeqInput_STRELKA2_INDEL}
        strelka_called_indel=$(zcat ${strelka_indel_vcf} | grep -v ^# | wc -l)
        if [[ -f "${strelka_indel_vcf}" ]]; then
            if (( ${strelka_called_snv} > 0 )); then 
                InputArgs_VCF+=("--strelka-indel ${SomaticSeqInput_STRELKA2_INDEL}")
            fi
        fi

        # muse2
        muse2_vcf=${SomaticSeqInput_MUSE2}
        muse2_called_vars=$(zcat ${muse2_vcf} | grep -v ^# | wc -l)
        if [[ -f "${muse2_vcf}" ]]; then
            if (( ${muse2_called_vars} > 0 )); then 
                InputArgs_VCF+=("--muse-vcf ${SomaticSeqInput_MUSE2}")
            fi
        fi

        # vardict
        vardict_vcf=${SomaticSeqInput_VARDICT}
        vardict_called_vars=$(zcat ${vardict_vcf} | grep -v ^# | wc -l)
        if [[ -f "${vardict_vcf}" ]]; then
            if (( ${vardict_called_vars} > 0 )); then 
                InputArgs_VCF+=("--vardict-vcf ${SomaticSeqInput_VARDICT}")
            fi
        fi
    #
    ## run somaticseq
    ${RunSomaticseq} \
        --output-directory ${SomaticseqDir} \
        --genome-reference ${GenomeFasta} \
        --threads ${Threads} \
        --inclusion-region ${TargetRegionBed} \
        paired \
        --tumor-bam-file ${TumorBamDir}/${TumorSeqID}.analysisReady.bam \
        --normal-bam-file ${NormalBamDir}/${NormalSeqID}.analysisReady.bam \
        --tumor-sample ${TumorSeqID} \
        --normal-sample ${NormalSeqID} \
        ${InputArgs_VCF[@]}


    ## add full-contigs into vcf
    ${RunGatk4} UpdateVCFSequenceDictionary \
        -V ${SomaticseqDir}/Consensus.sSNV.vcf \
        --source-dictionary ${GenomeFasta}.dict \
        --output ${SomaticseqDir}/Consensus.sSNV.contig.fix.vcf

    ${RunGatk4} UpdateVCFSequenceDictionary \
        -V ${SomaticseqDir}/Consensus.sINDEL.vcf \
        --source-dictionary ${GenomeFasta}.dict \
        --output ${SomaticseqDir}/Consensus.sINDEL.contig.fix.vcf

    ## sort vcf 
    bcftools sort -Oz -o ${SomaticseqDir}/Consensus.sSNV.contig.fix.vcf.gz ${SomaticseqDir}/Consensus.sSNV.contig.fix.vcf
    tabix -p vcf ${SomaticseqDir}/Consensus.sSNV.contig.fix.vcf.gz

    bcftools sort -Oz -o ${SomaticseqDir}/Consensus.sINDEL.contig.fix.vcf.gz ${SomaticseqDir}/Consensus.sINDEL.contig.fix.vcf
    tabix -p vcf ${SomaticseqDir}/Consensus.sINDEL.contig.fix.vcf.gz

    ## merge vcf
    bcftools concat -a -Oz -o ${SomaticseqDir}/${TumorSeqID}.somaticseq.consensus.vcf.gz \
        ${SomaticseqDir}/Consensus.sSNV.contig.fix.vcf.gz \
        ${SomaticseqDir}/Consensus.sINDEL.contig.fix.vcf.gz
    tabix -p vcf ${SomaticseqDir}/${TumorSeqID}.somaticseq.consensus.vcf.gz


    ## FINAL SOMATICSEQ VCF
    bcftools view -f PASS ${SomaticseqDir}/${TumorSeqID}.somaticseq.consensus.vcf.gz \
        | bcftools norm -m -any -f ${GenomeFasta} -Oz -o ${SomaticseqDir}/${TumorSeqID}.somaticseq.consensus.pass.vcf.gz
    tabix -p vcf "${SomaticseqDir}/${TumorSeqID}.somaticseq.consensus.pass.vcf.gz"

    # gnomad af annotation
    echo '##INFO=<ID=GNOMAD_AF,Number=A,Type=Float,Description="Allele frequency from gnomAD">' > ${SomaticseqDir}/gnomad_db.header
    echo '##INFO=<ID=GNOMAD_DB,Number=0,Type=Flag,Description="Present in gnomAD population database">' >> ${SomaticseqDir}/gnomad_db.header
    bcftools annotate \
        -a "${VCF_gnomad_mutect_annot}" \
        -h "${SomaticseqDir}/gnomad_db.header" \
        -c CHROM,POS,REF,ALT,INFO/GNOMAD_AF,INFO/GNOMAD_DB:=INFO/GNOMAD_AF \
        -Oz -o "${SomaticseqDir}/${TumorSeqID}.somaticseq.consensus.pass.gnomad.vcf.gz" \
        "${SomaticseqDir}/${TumorSeqID}.somaticseq.consensus.pass.vcf.gz"
    tabix -p vcf "${SomaticseqDir}/${TumorSeqID}.somaticseq.consensus.pass.gnomad.vcf.gz"
    
    ln -Tsf "${SomaticseqDir}/${TumorSeqID}.somaticseq.consensus.pass.gnomad.vcf.gz" "${TumorVcfDir}/${TumorSeqID}.${PrefixTag}.VARIANTS.vcf.gz"
    ln -Tsf "${SomaticseqDir}/${TumorSeqID}.somaticseq.consensus.pass.gnomad.vcf.gz.tbi" "${TumorVcfDir}/${TumorSeqID}.${PrefixTag}.VARIANTS.vcf.gz.tbi"

#---------------------------------------------------------------------------------------------------

    /storage/home/kangsm/runScripts/WES_v2/WES_M04-1_VARIANT-ANNOTATION.sh \
        --BASE_DIR "${BaseDir}" \
        --BATCH_ID "${BatchID}" \
        --TUMOR_SEQ_ID "${TumorSeqID}" \
        --NORMAL_SEQ_ID "${NormalSeqID}" \
        --ANNOTATE_TARGET_VCF_TYPE "NT_TS" \
        --GENOME_ASSEMBLY "${GenomeAssembly}" \
        --THREADS "${Threads}"


#-------------------------------------------------------------------------------

