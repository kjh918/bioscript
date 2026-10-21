

#! /bin/bash

#$ -N WES_HMFTOOLS
#$ -S /bin/bash
#$ -pe mpi 16
#$ -t 1-18
#$ -tc 8
#$ -e /storage/home/kangsm/analysis/gcx.pipeline.use.hmftools.in.wes/01.preprocess/runlog/wes/hmftools.run.$TASK_ID.err
#$ -o /storage/home/kangsm/analysis/gcx.pipeline.use.hmftools.in.wes/01.preprocess/runlog/wes/hmftools.run.$TASK_ID.out

#---/ Sample List File
    SampleListFile="/storage/home/kangsm/analysis/gcx.pipeline.use.hmftools.in.wes/01.preprocess/analysis.list"
#-------------------------------------------------------------------------------

for SGE_TASK_ID in $(seq 1 18)
do
#---/ Sample-Info from Sample-List-File |---------------------------------------
    # SeqID RGID BatchID WesLibrary
    RunInfo=$(sed -n -e "${SGE_TASK_ID} p" ${SampleListFile})
  
    TumorSeqID=$(echo ${RunInfo} | awk '{print $1}')
    NormalSeqID=$(echo ${RunInfo} | awk '{print $2}')
    BatchID=$(echo ${RunInfo} | awk '{print $3}')
    WesLibrary=$(echo ${RunInfo} | awk '{print $4}')

#-------------------------------------------------------------------------------

#---/ Manual Params |-----------------------------------------------------------
    BaseDir="/data/wes"
    #BatchID="WES-panel.of.normals-hg19-sureselect.v6"
    #SampleListFile="/data/pon/WES-panel.of.normals-hg19-sureselect.v6/meta/sample.list"
    Threads=16
    GenomeAssembly="hg19"
    #WesLibrary="sureselect.v6"
#-------------------------------------------------------------------------------

    HMF_GENOME_FASTA="/storage/references_and_index/hg19/hmftools/Homo_sapiens.GRCh37.GATK.illumina.fasta"
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
    VepGenome="GRCh37"
    RetainInfoItems="Refseq_Feature,Refseq_Gene,Refseq_EXON,Refseq_INTRON,Refseq_HGVSc,Refseq_HGVSp,GT,DP,BIALLELIC,BLACKLIST_BED,BLACKLIST_VCF,GND_FREQ,HOTSPOT,PON_COUNT,PON_MAX,PURPLE_AF,PURPLE_CN,PURPLE_MACN,PURPLE_VCN,TIER"

    HMFTOOLS_RES_BASE="${BaseDir}/${BatchID}/${TumorSeqID}/hmftools"

    RunMode="NT"
    OUTPUT_PURPLE="${HMFTOOLS_RES_BASE}/purple"
    InputPurpleSomaticVCF="${OUTPUT_PURPLE}/${TumorSeqID}.purple.somatic.vcf.gz"

    # output files prefix
    OutputVCF_Refseq="${InputPurpleSomaticVCF/vcf.gz/vep.refseq.vcf}"
    OutputVCF_Ensembl="${InputPurpleSomaticVCF/vcf.gz/vep.ensembl.vcf}"
    OutputVCF_Merged="${InputPurpleSomaticVCF/vcf.gz/vep.annot.vcf}"
    OutputVCF_Normalized="${InputPurpleSomaticVCF/vcf.gz/vep.annot.norm.vcf}"
    OutputMAF_Merged="${InputPurpleSomaticVCF/vcf.gz/vep.annot.norm.maf}"

    bcftools norm -f "${HMF_GENOME_FASTA}" -m -any "${OutputVCF_Merged}.gz" -Ov -o "${OutputVCF_Normalized}"
    ${RunVcf2Maf} --inhibit-vep \
        --input-vcf "${OutputVCF_Normalized}" \
        --output-maf "${OutputMAF_Merged}" \
        --tumor-id "${TumorSeqID}" \
        --normal-id "${NormalSeqID}" \
        --ref-fasta "${HMF_GENOME_FASTA}" \
        --retain-info "${RetainInfoItems}" \
        --ncbi-build "${VepGenome}" \
        --retain-ann "HGVSg,am_class,am_pathogenicity,ClinVar_CLNSIG,ClinVar_CLNDN,ClinVar_ORIGIN,gnomADe_AF,gnomADe_EAS_AF,MAX_AF,MAX_AF_POPS,TRANSCRIPTION_FACTORS,COSMIC,COSMIC_HGVSC,COSMIC_LEGACY_ID,COSMIC_TIER"

    RunMode="Tonly"
    OUTPUT_PURPLE="${HMFTOOLS_RES_BASE}/purple_Tonly"
    InputPurpleSomaticVCF="${OUTPUT_PURPLE}/${TumorSeqID}.purple.somatic.vcf.gz"

    # output files prefix
    OutputVCF_Refseq="${InputPurpleSomaticVCF/vcf.gz/vep.refseq.vcf}"
    OutputVCF_Ensembl="${InputPurpleSomaticVCF/vcf.gz/vep.ensembl.vcf}"
    OutputVCF_Merged="${InputPurpleSomaticVCF/vcf.gz/vep.annot.vcf}"
    OutputVCF_Normalized="${InputPurpleSomaticVCF/vcf.gz/vep.annot.norm.vcf}"
    OutputMAF_Merged="${InputPurpleSomaticVCF/vcf.gz/vep.annot.norm.maf}"
    bcftools norm -f "${HMF_GENOME_FASTA}" -m -any "${OutputVCF_Merged}.gz" -Ov -o "${OutputVCF_Normalized}"
    ${RunVcf2Maf} --inhibit-vep \
        --input-vcf ${OutputVCF_Merged} \
        --output-maf ${OutputMAF_Merged} \
        --tumor-id ${TumorSeqID} \
        --ref-fasta ${HMF_GENOME_FASTA} \
        --retain-info ${RetainInfoItems} \
        --retain-ann HGVSg,am_class,am_pathogenicity,ClinVar_CLNSIG,ClinVar_CLNDN,ClinVar_ORIGIN,gnomADe_AF,gnomADe_EAS_AF,MAX_AF,MAX_AF_POPS,TRANSCRIPTION_FACTORS,COSMIC,COSMIC_HGVSC,COSMIC_LEGACY_ID,COSMIC_TIER

#-------------------------------------------------------------------------------
done


