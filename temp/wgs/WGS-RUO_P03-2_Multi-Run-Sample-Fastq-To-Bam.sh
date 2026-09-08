
#!/bin/bash

#==============================================================================#
# ScriptID : GMC-NGS-B02-A04
# Author : kangsm
# version : v1.0 
# lastest update date : 2025.12.26
#==============================================================================#

### CMD Usage
    usage() {
        echo "
        Usage: WGS_P03-1_FASTQ-TO-BAM.sh [ ARGUMENTS ]...
        
        --BASE_DIR          [ REQUIRED ] analysis base folder by NGS data type
        
        --BATCH_ID          [ REQUIRED ] dataset batch-id
        
        --SEQ_ID            [ REQURIED ] sample-id
        
        --THREADS           parallel processing cpu cores 
                            default = 15
        
        --GENOME_ASSEMBLY   reference genome GENOME_ASSEMBLY version
                            default = hg38
        
        --MULTI_RUN_SAMPLE  multi-ngs-run sample
                            default = False 

        --RGID              [ REQURIED ] read-group-id list 

        --RGSM              read-group-sample-id
                            default = 'RGID'
                        
        --RGLB              read-group-library
                            default = WGS
        
        --RGPL              read-group-platform 
                            default = ILLUMINA
        
        --RGCN              read-group-center
                            default = GCX

        --SKIP_FASTQ_QC     skip fastq-qc step if already done
                            default = False

        -h | --help         Print this usage 
        "
    }
####

#### Option Args
    ARGS=$(getopt -a -o h: --long BATCH_ID:,SEQ_ID:,BASE_DIR:,THREADS:,GENOME_ASSEMBLY:,MULTI_RUN_SAMPLE:,RGID:,RGSM:,RGLB:,RGPL:,RGCN:,SKIP_FASTQ_QC:,help -- "$@" )
    VALID_ARGS=$?
    if [ "$VALID_ARGS" != "0" ]; then 
        usage >&2 
        exit 2
    fi

    eval set -- "$ARGS"
    while :
    do
        case "$1" in
            --BASE_DIR ) case "$2" in -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;; *) BaseDir="$2" ; shift 2 ;; esac ;;
            --BATCH_ID ) case "$2" in -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;; *) BatchID="$2" ; shift 2 ;; esac ;;
            --SEQ_ID ) case "$2" in -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;; *) SeqID="$2" ; shift 2 ;; esac ;;
            --THREADS ) case "$2" in -* | --* | "") echo "no value at $1, default = 15" ; shift ;; *) Threads="$2" ; shift 2 ;; esac ;;
            --GENOME_ASSEMBLY ) case "$2" in -* | --* | "") echo "no value at $1, default = hg19" ; shift ;; *) GenomeAssembly="$2" ; shift 2 ;; esac ;;
            --MULTI_RUN_SAMPLE ) case "$2" in -* | --* | "") echo "no value at $1, default = False" ; shift ;; *) MultiRunSample="$2" ; shift 2 ;; esac ;;
            --RGID ) case "$2" in -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;; *) RGID="$2" ; shift 2 ;; esac ;;
            --RGSM ) case "$2" in -* | --* | "") echo "no value at $1, default = 'RGID'" ; shift ;; *) RGSM=$2 ; shift 2 ;; esac ;;
            --RGLB ) case "$2" in -* | --* | "") echo "no value at $1, default = WGS" ; shift ;; *) RGLB=$2 ; shift 2 ;; esac ;;
            --RGPL ) case "$2" in -* | --* | "") echo "no value at $1, default = ILLUMINA" ; shift ;; *) RGPL=$2 ; shift 2 ;; esac ;;
            --RGCN ) case "$2" in -* | --* | "") echo "no value at $1, default = GCX" ; shift ;; *) RGCN=$2 ; shift 2 ;; esac ;;
            --SKIP_FASTQ_QC ) case "$2" in -* | --* | "") echo "no value at $1, default = False" ; shift ;; *) SkipFastqQC=$2 ; shift 2 ;; esac ;;
            -h | --help ) usage >&2 ; exit 2 ;;   
            --) shift ; break ;;
            *)  usage >&2 ; exit 2 ;;
        esac
    done
####

#### Check Params and Set Default 
    if [ -z "${BaseDir}" ]; then exit 2; fi
    if [ -z "${BatchID}" ]; then exit 2; fi
    if [ -z "${SeqID}" ]; then exit 2; fi
    if [ -z "${RGID}" ]; then exit 2; fi
    if [ -z "${Threads}" ]; then Threads=15; fi
    if [ -z "${GenomeAssembly}" ]; then GenomeAssembly="hg38"; fi
    if [ -z "${MultiRunSample}" ]; then MultiRunSample="False"; fi
    if [ -z "${RGSM}" ]; then RGSM=${RGID}; fi
    if [ -z "${RGLB}" ]; then RGLB="WGS"; fi
    if [ -z "${RGPL}" ]; then RGPL="ILLUMINA"; fi
    if [ -z "${RGCN}" ]; then RGCN="GCX"; fi
####

#### Manual Params 
    # BaseDir="/data/wgs/RUO_Service"
    # BatchID="WGS-26-02"
    # SeqID="WGS-26-02_01"
    # Threads=16
    # GenomeAssembly="hg19"
    # MultiRunSample="True"
    # RGID="23N72FLT4,23NCVJLT4"
    # RGSM="WGS-26-02_02"
    # RGPL="ILLUMINA"
    # RGLB="WGS"
    # RGCN="GCX"
    # SkipFastqQC="False"
####

#### Folders 
    DataDir=${BaseDir}/${BatchID}
    RawFastqDir=${DataDir}/fastq
    FastqDir=${DataDir}/${SeqID}/fastq
    BamDir=${DataDir}/${SeqID}/bam
    QcDir=${DataDir}/${SeqID}/qcfiles
    TmpDir=${DataDir}/${SeqID}/tmp
    LogDir=${DataDir}/${SeqID}/log
    mkdir -p ${FastqDir}
    mkdir -p ${BamDir}
    mkdir -p ${TmpDir}
    mkdir -p ${LogDir}
    mkdir -p ${QcDir}/qc_res
####

#### swtools
    export LC_ALL=C
    export LANG=C
    GATK_VERSION="4.6.2.0"
    RunSingularity="singularity exec -B /storage,/data /storage/images"
    RUN_FASTQ="${RunSingularity}/fastqc-0.12.1.sif fastqc"
    RUN_FASTP="${RunSingularity}/fastp-0.23.4.sif fastp"
    RUN_GATK4="${RunSingularity}/gatk-${GATK_VERSION}.sif gatk"
    RUN_BWA="${RunSingularity}/bwa-0.7.17.sif bwa mem"
    RUN_BWA_MEM2="${RunSingularity}/bwa-mem2-v2.2.1.sif bwa-mem2"
    HMF_RUN_REDUX="${RunSingularity}/hmftools-redux-v1.2.3.sif redux"
####

#### resources 
    ResourceDir="/storage/references_and_index/${GenomeAssembly}"
    GenomeFasta="${ResourceDir}/fasta/${GenomeAssembly}.primary_decoy_HLA_no.ALT.fa"
    BwaIndex="${ResourceDir}/fasta/${GenomeAssembly}.primary_decoy_HLA_no.ALT.fa"
    VCF_dbsnp="${ResourceDir}/vcf/${GenomeAssembly}.dbsnp138.vcf.gz"
    VCF_hcsnp="${ResourceDir}/vcf/${GenomeAssembly}.1000g.phase1.snps.high.confidence.vcf.gz"
    VCF_knwonindels="${ResourceDir}/vcf/${GenomeAssembly}.known.indels.vcf.gz"
    VCF_gsindels="${ResourceDir}/vcf/${GenomeAssembly}.gold.standard.indels.vcf.gz"
    VCF_gnomad="${ResourceDir}/vcf/${GenomeAssembly}.gnomad.af.only.vcf.gz"
    VCF_omni="${ResourceDir}/vcf/${GenomeAssembly}.omni.vcf.gz"
    VCF_hapmap="${ResourceDir}/vcf/${GenomeAssembly}.hapmap.vcf.gz"
    WgsCallTargetIntervals="${ResourceDir}/bed/${GenomeAssembly}.wgs.call.target.region.interval_list"
    
    if [ "${GenomeAssembly}" = "hg19" ]; then
        HMF_GenomeFasta="/storage/references_and_index/hg19/hmftools/Homo_sapiens.GRCh37.GATK.illumina.fasta"
        HMF_BwaMem2Index="/storage/references_and_index/hg19/hmftools/bwa-mem2_index-2.2.1/Homo_sapiens.GRCh37.GATK.illumina.fasta"
        HMF_GenomeVersion="37"
    else
        HMF_GenomeFasta="/storage/references_and_index/hg38/hmftools/GRCh38_masked_exclusions_alts_hlas.fasta"
        HMF_BwaMem2Index="/storage/references_and_index/hg38/hmftools/bwa-mem2_index-2.2.1/GRCh38_masked_exclusions_alts_hlas.fasta"
        HMF_GenomeVersion="38"
    fi
    HMF_ResourceDir="/storage/references_and_index/${GenomeAssembly}/hmftools/resources_links"
    HMF_UnmapRegion="${HMF_ResourceDir}/UNMAP-REGIONS.tsv"
    HMF_MsiJitterSites="${HMF_ResourceDir}/MSI-JITTER-SITES.tsv.gz"

####

#### RunLog 
    RunLog=${LogDir}/${SeqID}.wgs.fastq-to-bam.$(date '+%Y%m%d').log

    echo " " > ${RunLog}
    echo " WGS FASTQ-TO-BAM" >> ${RunLog}
    echo "----------------------------------------------------------------------" >> ${RunLog}
    echo " RUN DATE             $(date '+%Y-%m-%d %H:%M:%S')" >> ${RunLog}
    echo " BATCH ID             ${BatchID}" >> ${RunLog}
    echo " SEQ ID               ${SeqID}" >> ${RunLog}
    echo " READ GROUP ID        ${RGID}" >> ${RunLog}
    echo " SAMPLE NAME          ${RGSM}" >> ${RunLog}
    echo " REFERENCE VERSION    ${GenomeAssembly}" >> ${RunLog}
    echo " IS MULTI-RUN SAMPLES ${MultiRunSample}" >> ${RunLog}
    echo "----------------------------------------------------------------------" >> ${RunLog}
    echo " " >> ${RunLog}   
####

#### FASTQ-QC 
    if [[ "${SkipFastqQC}" = "False" || "${SkipFastqQC}" = "false" ]]; then

        if [[ "${MultiRunSample}" = "True" || "${MultiRunSample}" = "true"  ]]; then 
        
            echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Run Fastq-QC (Multi-Run-Samples)" >> ${RunLog}
        
            # multi-run-list (flowcell list)
            IFS=',' read -ra ReadGroupID <<< "$RGID"

            # fastq-qc each multi-run-sample
            for readgroupid in "${ReadGroupID[@]}"
            do
                # file-prefix
                RawFastqPrefix="${RawFastqDir}/${SeqID}-${readgroupid}"
                TrimFastqPrefix="${FastqDir}/${SeqID}-${readgroupid}.trimmed"
                QcPrefix="${QcDir}/${SeqID}-${readgroupid}"           
                # fastqc args
                FASTQC_ARGS=()
                FASTQC_ARGS+=("--extract")
                FASTQC_ARGS+=("-t ${Threads}")
                FASTQC_ARGS+=("-o ${QcDir}")
                FASTQC_ARGS+=("${RawFastqPrefix}_R1.fastq.gz")
                FASTQC_ARGS+=("${RawFastqPrefix}_R2.fastq.gz")
                # fastp args
                FASTP_ARGS=()
                FASTP_ARGS+=("--in1 ${RawFastqPrefix}_R1.fastq.gz")
                FASTP_ARGS+=("--in2 ${RawFastqPrefix}_R2.fastq.gz")
                FASTP_ARGS+=("--json ${QcPrefix}.fastp.json")
                FASTP_ARGS+=("--html ${QcPrefix}.fastp.html")
                FASTP_ARGS+=("--out1 ${TrimFastqPrefix}_R1.fastq.gz")
                FASTP_ARGS+=("--out2 ${TrimFastqPrefix}_R2.fastq.gz")
                FASTP_ARGS+=("--thread ${Threads}")
                FASTP_ARGS+=("--length_required 50")
                FASTP_ARGS+=("--cut_tail")
                FASTP_ARGS+=("--cut_tail_window_size 4")
                FASTP_ARGS+=("--cut_tail_mean_quality 20")
                FASTP_ARGS+=("--trim_poly_g")
                FASTP_ARGS+=("--detect_adapter_for_pe")
                FASTP_ARGS+=("--n_base_limit 5")
                # run fastqc 
                echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Run fastqc : ${SeqID}-${readgroupid}" >> ${RunLog}
                ${RUN_FASTQ} ${FASTQC_ARGS[@]} 
                # run fastp
                echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Run fastp : ${SeqID}-${readgroupid}" >> ${RunLog}
                ${RUN_FASTP} ${FASTP_ARGS[@]}
                #            
            done
        
        else

            echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Run Fastq-QC (Single-Run-Sample)" >> ${RunLog}
            # file-prefix
            RawFastqPrefix="${RawFastqDir}/${SeqID}"
            TrimFastqPrefix="${FastqDir}/${SeqID}.trimmed"
            QcPrefix="${QcDir}/${SeqID}"
            # fastqc args
            FASTQC_ARGS=()
            FASTQC_ARGS+=("--extract")
            FASTQC_ARGS+=("-t ${Threads}")
            FASTQC_ARGS+=("-o ${QcDir}")
            FASTQC_ARGS+=("${RawFastqPrefix}_R1.fastq.gz")
            FASTQC_ARGS+=("${RawFastqPrefix}_R2.fastq.gz")
            # fastp args
            FASTP_ARGS=()
            FASTP_ARGS+=("--in1 ${RawFastqPrefix}_R1.fastq.gz")
            FASTP_ARGS+=("--in2 ${RawFastqPrefix}_R2.fastq.gz")
            FASTP_ARGS+=("--json ${QcPrefix}.fastp.json")
            FASTP_ARGS+=("--html ${QcPrefix}.fastp.html")
            FASTP_ARGS+=("--out1 ${TrimFastqPrefix}_R1.fastq.gz")
            FASTP_ARGS+=("--out2 ${TrimFastqPrefix}_R2.fastq.gz")
            FASTP_ARGS+=("--thread ${Threads}")
            FASTP_ARGS+=("--length_required 50")
            FASTP_ARGS+=("--cut_tail")
            FASTP_ARGS+=("--cut_tail_window_size 4")
            FASTP_ARGS+=("--cut_tail_mean_quality 20")
            FASTP_ARGS+=("--trim_poly_g")
            FASTP_ARGS+=("--detect_adapter_for_pe")
            FASTP_ARGS+=("--n_base_limit 5")
            # run fastqc 
            echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Run fastqc : ${SeqID}" >> ${RunLog}
            ${RUN_FASTQ} ${FASTQC_ARGS[@]} 
            # run fastp
            echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Run fastqc : ${SeqID}" >> ${RunLog}
            ${RUN_FASTP} ${FASTP_ARGS[@]}

        fi
    
    fi
####

#### HMFTOOLS pipeline
    if [[ "${MultiRunSample}" = "True" || "${MultiRunSample}" = "true"  ]];

        echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] HMFTOOLS pipeline -- Multi-Run-Samples" >> ${RunLog}
        #
        for readgroupid in "${ReadGroupID[@]}"
        do
            # file prefix
            TrimFastqPrefix="${FastqDir}/${SeqID}-${readgroupid}.trimmed"
            BamPrefix="${BamDir}/${SeqID}-${readgroupid}"
            BwaThreads=12
            SortThreads=4
            # bwa-mem2 args
            BWA_MEM2_ARGS=()
            BWA_MEM2_ARGS+=("-M")
            BWA_MEM2_ARGS+=("-Y")
            BWA_MEM2_ARGS+=("-q")
            BWA_MEM2_ARGS+=("-K 100000000")
            BWA_MEM2_ARGS+=("-t ${BwaThreads}")
            BWA_MEM2_ARGS+=("-t ${BwaThreads}")
            BWA_MEM2_ARGS+=("-R @RG\tID:${readgroupid}\tSM:${SeqID}\tLB:WGS\tPL:ILLUMINA\tCN:GCX")
            BWA_MEM2_ARGS+=("${HMF_BwaMem2Index}")
            BWA_MEM2_ARGS+=("${TrimFastqPrefix}_R1.fastq.gz")
            BWA_MEM2_ARGS+=("${TrimFastqPrefix}_R2.fastq.gz")
            # samtools args
            SAMTOOLS_ARGS=()
            SAMTOOLS_ARGS+=("-@ ${SortThreads}")
            SAMTOOLS_ARGS+=("-o ${BamPrefix}.bwa.mem2.bam")
            # run bwa-mem2
            echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Run BWA-mem2 : ${SeqID}-${readgroupid}" >> ${RunLog}
            ${RUN_BWA_MEM2} ${BWA_MEM2_ARGS[@]} | samtools sort ${SAMTOOLS_ARGS[@]} -
            #
        done
        #
        # merge bwa-mem2 multi-run-bams args
        BWA_MEM2_BAM_MERGE_ARGS=()
        BWA_MEM2_BAM_MERGE_ARGS+=("--REFERENCE_SEQUENCE ${HMF_GenomeFasta}")
        BWA_MEM2_BAM_MERGE_ARGS+=("--CREATE_INDEX true")
        BWA_MEM2_BAM_MERGE_ARGS+=("--SORT_ORDER coordinate")
        BWA_MEM2_BAM_MERGE_ARGS+=("--USE_THREADING true")
        BWA_MEM2_BAM_MERGE_ARGS+=("--VALIDATION_STRINGENCY LENIENT")
        BWA_MEM2_BAM_MERGE_ARGS+=("--TMP_DIR ${TmpDir}")
        BWA_MEM2_BAM_MERGE_ARGS+=("--OUTPUT ${BamDir}/${SeqID}.bwa.mem2.bam")
        for readgroupid in ${ReadGroupID[@]}
        do
            BWA_MEM2_BAM_MERGE_ARGS+=("-I ${BamDir}/${SeqID}-${readgroupid}.bwa.mem2.bam")
        done
        # run bam merge
        echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Run Multi-Run-BWA-mem2 Bams Merge : ${SeqID}" >> ${RunLog}
        ${RUN_GATK4} --java-options "-Xmx24g" MergeSamFiles ${BWA_MEM2_BAM_MERGE_ARGS[@]}
        #
    else        
        echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] HMFTOOLS pipeline -- Single-Run-Samples" >> ${RunLog}
        #
        # file prefix
        TrimFastqPrefix="${FastqDir}/${SeqID}.trimmed"
        BamPrefix="${BamDir}/${SeqID}"
        BwaThreads=12
        SortThreads=4
        # bwa-mem2 args
        BWA_MEM2_ARGS=()
        BWA_MEM2_ARGS+=("-M")
        BWA_MEM2_ARGS+=("-Y")
        BWA_MEM2_ARGS+=("-q")
        BWA_MEM2_ARGS+=("-K 100000000")
        BWA_MEM2_ARGS+=("-t ${BwaThreads}")
        BWA_MEM2_ARGS+=("-t ${BwaThreads}")
        BWA_MEM2_ARGS+=("-R @RG\tID:${readgroupid}\tSM:${SeqID}\tLB:WGS\tPL:ILLUMINA\tCN:GCX")
        BWA_MEM2_ARGS+=("${HMF_BwaMem2Index}")
        BWA_MEM2_ARGS+=("${TrimFastqPrefix}_R1.fastq.gz")
        BWA_MEM2_ARGS+=("${TrimFastqPrefix}_R2.fastq.gz")
        # samtools args
        SAMTOOLS_ARGS=()
        SAMTOOLS_ARGS+=("-@ ${SortThreads}")
        SAMTOOLS_ARGS+=("-o ${BamPrefix}.bwa.mem2.bam")
        # run bwa-mem2
        echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Run BWA-mem2 : ${SeqID}" >> ${RunLog}
        ${RUN_BWA_MEM2} ${BWA_MEM2_ARGS[@]} | samtools sort ${SAMTOOLS_ARGS[@]} -
        #
    fi
    #
    # hmftools redux 
        #
        HMF_ReduxOutDir=${DataDir}/${SeqID}/hmftools/01.redux
        mkdir -p ${HMF_ReduxOutDir}
        #
        HMF_REDUX_ARGS=()
        HMF_REDUX_ARGS+=("-sample ${SeqID}")
        HMF_REDUX_ARGS+=("-input_bam ${BamDir}/${SeqID}.bwa.mem2.bam")
        HMF_REDUX_ARGS+=("-output_bam ${HMF_ReduxOutDir}/${SeqID}.redux.bam")
        HMF_REDUX_ARGS+=("-output_bam ${HMF_ReduxOutDir}")
        HMF_REDUX_ARGS+=("-ref_genome_version ${HMF_GenomeVersion}")
        HMF_REDUX_ARGS+=("-ref_genome ${HMF_GenomeFasta}")
        HMF_REDUX_ARGS+=("-unmap_regions ${HMF_UnmapRegion}")
        HMF_REDUX_ARGS+=("-ref_genome_msi_file ${HMF_MsiJitterSites}")
        HMF_REDUX_ARGS+=("-form_consensus")
        HMF_REDUX_ARGS+=("-write_stats")
        HMF_REDUX_ARGS+=("-bamtool /opt/conda/bin/samtools")
        HMF_REDUX_ARGS+=("-threads ${Threads}")
        #
        ${HMF_RUN_REDUX} ${HMF_REDUX_ARGS[@]}
    # hmftools bam-qc (gatk wgsmetric)
        #
        WGS_METRIC_ARGS=()
        WGS_METRIC_ARGS+=("--INPUT ${HMF_ReduxOutDir}/${SeqID}.redux.bam")
        WGS_METRIC_ARGS+=("--OUTPUT ${HMF_ReduxOutDir}/${SeqID}.redux.wgs.metrics.txt")
        WGS_METRIC_ARGS+=("--REFERENCE_SEQUENCE ${HMF_GenomeFasta}")
        WGS_METRIC_ARGS+=("--VALIDATION_STRINGENCY LENIENT")
        #
        ${RUN_GATK4} --java-options "-Xmx16g" CollectWgsMetrics ${WGS_METRIC_ARGS[@]}

####

#### GATK pipeline
    if [[ "${MultiRunSample}" = "True" || "${MultiRunSample}" = "true"  ]]; then
        echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] GATK4 pipeline -- Multi-Run-Samples" >> ${RunLog}
        #
        for readgroupid in "${ReadGroupID[@]}"
        do
            # file prefix
            TrimFastqPrefix="${FastqDir}/${SeqID}-${readgroupid}.trimmed"
            BamPrefix="${BamDir}/${SeqID}-${readgroupid}"
            QcPrefix="${QcDir}/${SeqID}-${readgroupid}"
            # u-bam args
            FASTQ_TO_SAM_ARGS=()
            FASTQ_TO_SAM_ARGS+=("--FASTQ ${TrimFastqPrefix}_R1.fastq.gz")
            FASTQ_TO_SAM_ARGS+=("--FASTQ2 ${TrimFastqPrefix}_R2.fastq.gz")
            FASTQ_TO_SAM_ARGS+=("--OUTPUT ${BamPrefix}.fastqtosam.bam")
            FASTQ_TO_SAM_ARGS+=("--READ_GROUP_NAME ${readgroupid}")
            FASTQ_TO_SAM_ARGS+=("--SAMPLE_NAME ${RGSM}")
            FASTQ_TO_SAM_ARGS+=("--LIBRARY_NAME ${RGLB}")
            FASTQ_TO_SAM_ARGS+=("--PLATFORM ${RGPL}")
            FASTQ_TO_SAM_ARGS+=("--SEQUENCING_CENTER ${RGCN}")
            FASTQ_TO_SAM_ARGS+=("--QUALITY_FORMAT Standard")
            FASTQ_TO_SAM_ARGS+=("--TMP_DIR ${TmpDir}")
            # bwa-bam args
            BWA_ARGS=()
            BWA_ARGS+=("-M")
            BWA_ARGS+=("-Y")
            BWA_ARGS+=("-q")
            BWA_ARGS+=("-K 100000000")
            BWA_ARGS+=("-t ${Threads}")
            BWA_ARGS+=("${BwaIndex}")
            BWA_ARGS+=("${TrimFastqPrefix}_R1.fastq.gz")
            BWA_ARGS+=("${TrimFastqPrefix}_R2.fastq.gz")
            # merged-bam args
            MERGE_BAM_ALIGNS_ARGS=()
            MERGE_BAM_ALIGNS_ARGS+=("--REFERENCE_SEQUENCE ${GenomeFasta}")
            MERGE_BAM_ALIGNS_ARGS+=("--UNMAPPED_BAM ${BamPrefix}.fastqtosam.bam")
            MERGE_BAM_ALIGNS_ARGS+=("--ALIGNED_BAM ${BamPrefix}.bwa.mem.sam")
            MERGE_BAM_ALIGNS_ARGS+=("--OUTPUT ${BamPrefix}.merged.bam")
            MERGE_BAM_ALIGNS_ARGS+=("--EXPECTED_ORIENTATIONS FR")
            MERGE_BAM_ALIGNS_ARGS+=("--PRIMARY_ALIGNMENT_STRATEGY BestMapq")
            MERGE_BAM_ALIGNS_ARGS+=("--SORT_ORDER coordinate")
            MERGE_BAM_ALIGNS_ARGS+=("--CREATE_INDEX true")
            MERGE_BAM_ALIGNS_ARGS+=("--INCLUDE_SECONDARY_ALIGNMENTS true")
            MERGE_BAM_ALIGNS_ARGS+=("--ALIGNER_PROPER_PAIR_FLAGS false")
            MERGE_BAM_ALIGNS_ARGS+=("--CLIP_ADAPTERS true")
            MERGE_BAM_ALIGNS_ARGS+=("--ADD_MATE_CIGAR true")
            MERGE_BAM_ALIGNS_ARGS+=("--MAX_INSERTIONS_OR_DELETIONS 1")
            MERGE_BAM_ALIGNS_ARGS+=("--VALIDATION_STRINGENCY LENIENT")
            MERGE_BAM_ALIGNS_ARGS+=("--TMP_DIR ${TmpDir}")
            # run fastq-to-sam
            echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Run GATK FastqToSam : ${SeqID}-${readgroupid}" >> ${RunLog}
            ${RUN_GATK4} --java-options "-Xmx8g" FastqToSam ${FASTQ_TO_SAM_ARGS[@]}
            # run bwa
            echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Run BWA-mem : ${SeqID}-${readgroupid}" >> ${RunLog}
            ${RUN_BWA} ${BWA_ARGS[@]} > ${BamPrefix}.bwa.mem.sam 
            # run merge-bam-alignment
            echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Run GATK MergeBamAlignment : ${SeqID}-${readgroupid}" >> ${RunLog}
            ${RUN_GATK4} --java-options "-Xmx24g" MergeBamAlignment ${MERGE_BAM_ALIGNS_ARGS[@]}
            #
        done
        #
        # merge multi-run-bams args
        GATK_BAM_MERGE_ARGS=("--REFERENCE_SEQUENCE ${GenomeFasta}")
        GATK_BAM_MERGE_ARGS+=("--CREATE_INDEX true")
        GATK_BAM_MERGE_ARGS+=("--SORT_ORDER coordinate")
        GATK_BAM_MERGE_ARGS+=("--USE_THREADING true")
        GATK_BAM_MERGE_ARGS+=("--VALIDATION_STRINGENCY LENIENT")
        GATK_BAM_MERGE_ARGS+=("--TMP_DIR ${TmpDir}")
        GATK_BAM_MERGE_ARGS+=("--OUTPUT ${BamDir}/${SeqID}.merged.bam")
        for readgroupid in "${ReadGroupID[@]}"
        do
            GATK_BAM_MERGE_ARGS+=("-I ${BamDir}/${SeqID}-${readgroupid}.merged.bam")
        done
        # run bam merge
        echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Run GATK MergeSamFiles : ${SeqID}" >> ${RunLog}
        ${RUN_GATK4} --java-options "-Xmx24g" MergeSamFiles ${GATK_BAM_MERGE_ARGS[@]}
        # 
        # remove multi-run tmp-files
        echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Remove multi-run-sample temp-files" >> ${RunLog}
        for readgroupid in "${ReadGroupID[@]}"
        do
            rm ${BamDir}/${SeqID}-${readgroupid}.bwa.mem.sam
            rm ${BamDir}/${SeqID}-${readgroupid}.fastqtosam.bam
            rm ${BamDir}/${SeqID}-${readgroupid}.merged.bam
            rm ${BamDir}/${SeqID}-${readgroupid}.merged.bai
        done 
        #
    else
        # file prefix
        TrimFastqPrefix="${FastqDir}/${SeqID}.trimmed"
        BamPrefix="${BamDir}/${SeqID}"
        # u-bam args
        FASTQ_TO_SAM_ARGS=()
        FASTQ_TO_SAM_ARGS+=("--FASTQ ${TrimFastqPrefix}_R1.fastq.gz")
        FASTQ_TO_SAM_ARGS+=("--FASTQ2 ${TrimFastqPrefix}_R2.fastq.gz")
        FASTQ_TO_SAM_ARGS+=("--OUTPUT ${BamPrefix}.fastqtosam.bam")
        FASTQ_TO_SAM_ARGS+=("--READ_GROUP_NAME ${readgroupid}")
        FASTQ_TO_SAM_ARGS+=("--SAMPLE_NAME ${RGSM}")
        FASTQ_TO_SAM_ARGS+=("--LIBRARY_NAME ${RGLB}")
        FASTQ_TO_SAM_ARGS+=("--PLATFORM ${RGPL}")
        FASTQ_TO_SAM_ARGS+=("--SEQUENCING_CENTER ${RGCN}")
        FASTQ_TO_SAM_ARGS+=("--QUALITY_FORMAT Standard")
        FASTQ_TO_SAM_ARGS+=("--TMP_DIR ${TmpDir}")
        # bwa-bam args
        BWA_ARGS=()
        BWA_ARGS+=("-M")
        BWA_ARGS+=("-Y")
        BWA_ARGS+=("-q")
        BWA_ARGS+=("-K 100000000")
        BWA_ARGS+=("-t ${Threads}")
        BWA_ARGS+=("${BwaIndex}")
        BWA_ARGS+=("${TrimFastqPrefix}_R1.fastq.gz")
        BWA_ARGS+=("${TrimFastqPrefix}_R2.fastq.gz")
        # merged-bam args
        MERGE_BAM_ALIGNS_ARGS=()
        MERGE_BAM_ALIGNS_ARGS+=("--REFERENCE_SEQUENCE ${GenomeFasta}")
        MERGE_BAM_ALIGNS_ARGS+=("--UNMAPPED_BAM ${BamPrefix}.fastqtosam.bam")
        MERGE_BAM_ALIGNS_ARGS+=("--ALIGNED_BAM ${BamPrefix}.bwa.mem.sam")
        MERGE_BAM_ALIGNS_ARGS+=("--OUTPUT ${BamPrefix}.merged.bam")
        MERGE_BAM_ALIGNS_ARGS+=("--EXPECTED_ORIENTATIONS FR")
        MERGE_BAM_ALIGNS_ARGS+=("--PRIMARY_ALIGNMENT_STRATEGY BestMapq")
        MERGE_BAM_ALIGNS_ARGS+=("--SORT_ORDER coordinate")
        MERGE_BAM_ALIGNS_ARGS+=("--CREATE_INDEX true")
        MERGE_BAM_ALIGNS_ARGS+=("--INCLUDE_SECONDARY_ALIGNMENTS true")
        MERGE_BAM_ALIGNS_ARGS+=("--ALIGNER_PROPER_PAIR_FLAGS false")
        MERGE_BAM_ALIGNS_ARGS+=("--CLIP_ADAPTERS true")
        MERGE_BAM_ALIGNS_ARGS+=("--ADD_MATE_CIGAR true")
        MERGE_BAM_ALIGNS_ARGS+=("--MAX_INSERTIONS_OR_DELETIONS 1")
        MERGE_BAM_ALIGNS_ARGS+=("--VALIDATION_STRINGENCY LENIENT")
        MERGE_BAM_ALIGNS_ARGS+=("--TMP_DIR ${TmpDir}")
        # run uBAM
        echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Run GATK FastqToSam : ${SeqID}" >> ${RunLog}
        ${RUN_GATK4} --java-options "-Xmx8g" FastqToSam ${FASTQ_TO_SAM_ARGS[@]}
        # run aBAM (bwa)
        echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Run BWA-mem : ${SeqID}" >> ${RunLog}
        ${RUN_BWA} ${BWA_ARGS[@]} > ${BamPrefix}.bwa.mem.sam 
        # run merge-bam-alignment
        echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Run GATK MergeBamAlignment : ${SeqID}" >> ${RunLog}
        ${RUN_GATK4} --java-options "-Xmx24g" MergeBamAlignment ${MERGE_BAM_ALIGNS_ARGS[@]}
        #
    fi
    #
    # dedup args
    MARK_DUPS_ARGS=()
    MARK_DUPS_ARGS+=("--INPUT ${BamDir}/${SeqID}.merged.bam")
    MARK_DUPS_ARGS+=("--METRICS_FILE ${QcDir}/${SeqID}.mark.duplicates.metrics.txt")
    MARK_DUPS_ARGS+=("--OUTPUT ${BamDir}/${SeqID}.merged.dup.marked.bam")
    MARK_DUPS_ARGS+=("--CREATE_INDEX true")
    MARK_DUPS_ARGS+=("--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500")
    MARK_DUPS_ARGS+=("--TAGGING_POLICY All")
    MARK_DUPS_ARGS+=("--MAX_RECORDS_IN_RAM 2000000")
    # base-recal args
    BASE_RECAL_ARGS=()
    BASE_RECAL_ARGS+=("--reference ${GenomeFasta}")
    BASE_RECAL_ARGS+=("--input ${BamDir}/${SeqID}.merged.dup.marked.bam")
    BASE_RECAL_ARGS+=("--known-sites ${VCF_dbsnp}")
    BASE_RECAL_ARGS+=("--known-sites ${VCF_hcsnp}")
    BASE_RECAL_ARGS+=("--known-sites ${VCF_knwonindels}")
    BASE_RECAL_ARGS+=("--known-sites ${VCF_gsindels}")
    BASE_RECAL_ARGS+=("--known-sites ${VCF_gnomad}")
    BASE_RECAL_ARGS+=("--known-sites ${VCF_omni}")
    BASE_RECAL_ARGS+=("--known-sites ${VCF_hapmap}")
    BASE_RECAL_ARGS+=("--intervals ${WgsCallTargetIntervals}")
    BASE_RECAL_ARGS+=("--output ${QcDir}/${SeqID}.recal.table")
    # apply-bqrs args
    APPLY_BQSR_ARGS=()
    APPLY_BQSR_ARGS+=("--reference ${GenomeFasta}")
    APPLY_BQSR_ARGS+=("--input ${BamDir}/${SeqID}.merged.dup.marked.bam")
    APPLY_BQSR_ARGS+=("--bqsr-recal-file ${QcDir}/${SeqID}.recal.table")
    APPLY_BQSR_ARGS+=("--intervals ${WgsCallTargetIntervals}")
    APPLY_BQSR_ARGS+=("-O ${BamDir}/${SeqID}.merged.dup.marked.recal.bam")
    # run mark-dups
    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Run GATK MarkDuplicates : ${SeqID}" >> ${RunLog}
    ${RUN_GATK4} --java-options "-Xmx32g" MarkDuplicates ${MARK_DUPS_ARGS[@]}
    # run base-recalibration
    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Run GATK BaseRecalibrator : ${SeqID}" >> ${RunLog} 
    ${RUN_GATK4} --java-options "-Xmx16g" BaseRecalibrator ${BASE_RECAL_ARGS[@]}
    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Run GATK ApplyBQSR : ${SeqID}" >> ${RunLog}
    ${RUN_GATK4} --java-options "-Xmx16g" ApplyBQSR ${APPLY_BQSR_ARGS[@]}
    #
    # analysisReady BAM
    ANALYSIS_READY_BAM=${BamDir}/${SeqID}.analysisReady.bam
    cp ${BamDir}/${SeqID}.merged.dup.marked.recal.bam ${ANALYSIS_READY_BAM}
    samtools index -@ ${Threads} ${BamDir}/${SeqID}.analysisReady.bam
    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] AnalysisReady BAM prepared : ${SeqID}" >> ${RunLog}
    # remove tmp
    if [[ -s ${ANALYSIS_READY_BAM} ]]; then
        rm ${BamDir}/${SeqID}.merged.bam   
        rm ${BamDir}/${SeqID}.merged.bai 
        rm ${BamDir}/${SeqID}.merged.dup.marked.bam
        rm ${BamDir}/${SeqID}.merged.dup.marked.bai
        rm ${BamDir}/${SeqID}.merged.dup.marked.recal.bam
        rm ${BamDir}/${SeqID}.bwa.mem.sam
        rm ${BamDir}/${SeqID}.fastqtosam.bam
        rm -r ${BamDir}/bqsr_tmp
    fi
####




####
    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] WGS FASTQ-TO-BAM PROCESS FINISHED." >> ${RunLog}
    echo " " >> ${RunLog}
####

