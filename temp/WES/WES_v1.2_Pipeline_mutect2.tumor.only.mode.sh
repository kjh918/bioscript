
#!/bin/bash

#---| Input Options |
    usage() {
        echo "Usage: GCX-WES-TumorOnlyVariantCalling.sh [ options ]
        [ -d | --seqFolder ]        : data processing base directory
        [ -s | --seqID ]       : sample ID (Used as sample's processing product folder name)
        [ -l | --logFile ]        : logFile name
        [ -h | --help ]           : Print this usage 
        [ --RGID ]                : ReadGroup ID (default = GCX)
        [ --RGLB ]                : ReadGroup Library (default = twistexome2)
        [ --RGPU ]                : ReadGroup Platform Unit (default = NULL)
        [ --RGPL ]                : ReadGroup Platform (default = illumina)
        [ --RGCN ]                : ReadGroup Platform (default = GCX)
        "
    }
  #---| Set Defualt Values |
    RGID="GCX" && RGLB="twistexome2" && RGPL="ILLUMINA" && RGCN="GCX" && RGPU=NULL
  #---/
#---/

#---| Parse Option Arguments |
    ARGS=$(getopt -a -o d:s:l:h: --long seqFolder:,seqID:,RGID:,RGLB:,RGPL:,RGPU:,RGCN:,logFile:,help -- "$@" )
    VALID_ARGS=$?
    if [ "$VALID_ARGS" != "0" ]; then 
        usage >&2 
        exit 2
    fi

    eval set -- "$ARGS"
    while :
    do
    case "$1" in
        -d | --seqFolder )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  SEQ_FOLDER=$2 ; shift 2 ;;
            esac ;;
        -s | --seqID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  SEQ_ID=$2 ; shift 2 ;;
            esac  ;;
        --RGID )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = GCX" ; shift ;;
                *)  RGID=$2 ; shift 2 ;;
            esac ;;
        --RGLB )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = twistexome2" ; shift ;;
                *)  RGLB=$2 ; shift 2 ;;
            esac ;;
        --RGPL )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = illumina" ; shift ;;
                *)  RGPL=$2 ; shift 2 ;;
            esac ;;
        --RGPU )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = NULL" ; shift ;;
                *)  RGPU=$2 ; shift 2 ;;
            esac ;;
        --RGCN )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = GCX" ; shift ;;
                *)  RGCN=$2 ; shift 2 ;;
            esac ;;
        -l | --logFile )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  LOG_NAME=$2 ; shift 2 ;;
            esac ;;
        -h | --help )
        usage >&2 ; exit 2 ;;   
        --) shift ; break ;;
        *)  usage >&2 ; exit 2 ;;
    esac
    done
#---/

#---| FOLDERS |
    DATA_DIR=/data/wes/${SEQ_FOLDER}
    mkdir -p ${DATA_DIR}/${SEQ_ID}/bam
    mkdir -p ${DATA_DIR}/${SEQ_ID}/vcf
    mkdir -p ${DATA_DIR}/${SEQ_ID}/qcfiles
    mkdir -p ${DATA_DIR}/${SEQ_ID}/qcfiles/qc_res
    mkdir -p ${DATA_DIR}/${SEQ_ID}/tmp
  ##
    FASTQ_DIR=${DATA_DIR}/${SEQ_ID}/fastq
    BAM_DIR=${DATA_DIR}/${SEQ_ID}/bam
    VCF_DIR=${DATA_DIR}/${SEQ_ID}/vcf
    QC_DIR=${DATA_DIR}/${SEQ_ID}/qcfiles
    TMP_DIR=${DATA_DIR}/${SEQ_ID}/tmp
#---/
#---| REFERENCES
    REF_DIR=/storage/references_and_index/hg19
    BWA_INDEX=${REF_DIR}/fasta/human_g1k_v37_decoy.fasta
    GENOME_FASTA=${REF_DIR}/fasta/human_g1k_v37_decoy.fasta
    BAIT_BED=${REF_DIR}/bed/Twist_WES_2.0/probe.bed
    TARGET_BED=${REF_DIR}/bed/Twist_WES_2.0/target.bed
    BAIT_INTERVAL=${REF_DIR}/bed/Twist_WES_2.0/probe.interval
    TARGET_INTERVAL=${REF_DIR}/bed/Twist_WES_2.0/target.interval
  ##
    VCF_DBSNP=${REF_DIR}/vcf/dbsnp_138.b37.vcf.gz
    VCF_1KG_SNP=${REF_DIR}/vcf/1000G_phase1.snps.high_confidence.b37.vcf.gz
    VCF_GNOMAD=${REF_DIR}/vcf/af-only-gnomad.raw.sites.grch37.vcf.gz
    VCF_INDEL=${REF_DIR}/vcf/Homo_sapiens_assembly19.known_indels.vcf.gz
    VCF_1KG_INDEL=${REF_DIR}/vcf/Mills_and_1000G_gold_standard.indels.b37.vcf.gz
    VCF_PON=${REF_DIR}/vcf/gatk-practice-mutect2-panel-of-normal.vcf.gz
#---/
#---| READGROUP |
    RGSM=${SEQ_ID}
#---/
#---| READGROUP |
    LOGFILE=/data/wes/${SEQ_FOLDER}/${SEQ_ID}/${LOG_NAME}
    touch ${LOGFILE}
#---/

#---\ 
    echo " " >> ${LOGFILE}
    echo "GCX WES processing log" >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " SEQ FOLDER           : ${DATA_DIR}" >> ${LOGFILE}
    echo " SEQ ID               : ${SEQ_ID}" >> ${LOGFILE}
    echo " READ GROUP ID        : ${RGID}" >> ${LOGFILE}
    echo " LIBRARY              : ${RGLB}" >> ${LOGFILE}
    echo " REFERENCE            : hg19 (GRCh37)" >> ${LOGFILE}
    echo " VARIANT CALLING MODE : Mutect2 Tumor-Only Mode" >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo "$(date) | Start Processing." >> ${LOGFILE}
#---/
#---| FASTQC |
#---\
    echo "$(date) | FastQC Start." >> ${LOGFILE}
#---/
    /storage/apps/bin/fastqc -t 2 -o ${QC_DIR} ${FASTQ_DIR}/*.fastq.gz
#---\
    echo "$(date) | FastQC Finished." >> ${LOGFILE}
#---/
#---| FASTP FILTERING |
#---\
    echo "$(date) | Fastq Trimming and Filtering Start." >> ${LOGFILE}
    echo "                                  - minimun length cut-off = 105" >> ${LOGFILE}
    echo "                                  - base q-score cut-off   = 20" >> ${LOGFILE}
#---/
    /storage/apps/bin/fastp \
        --in1 ${FASTQ_DIR}/${SEQ_ID}_R1.fastq.gz --in2 ${FASTQ_DIR}/${SEQ_ID}_R2.fastq.gz \
        --json ${DATA_DIR}/${SEQ_ID}/qcfiles/${SEQ_ID}.fastp.json \
        --length_required 105 \
        --qualified_quality_phred 20 \
        --out1 ${FASTQ_DIR}/${SEQ_ID}.trimmed_R1.fastq.gz --out2 ${FASTQ_DIR}/${SEQ_ID}.trimmed_R2.fastq.gz
#---\
    echo "$(date) | Fastq Trimming and Filtering Finished." >> ${LOGFILE}
#---/
#---| uBMA (Unligned BAM) GENERATION |
#---\
    echo "$(date) | uBAM Generation Start." >> ${LOGFILE}
#---/ 
    java -XX:ParallelGCThreads=20 -Xmx16G -jar /storage/apps/bin/picard.jar FastqToSam \
        --FASTQ ${FASTQ_DIR}/${SEQ_ID}.trimmed_R1.fastq.gz \
        --FASTQ2 ${FASTQ_DIR}/${SEQ_ID}.trimmed_R2.fastq.gz \
        --OUTPUT ${BAM_DIR}/${SEQ_ID}.fastqtosam.bam \
        --READ_GROUP_NAME ${RGID} \
        --SAMPLE_NAME ${RGSM} \
        --LIBRARY_NAME ${RGLB} \
        --PLATFORM ${RGPL} \
        --SEQUENCING_CENTER ${RGCN}
#---\
    echo "$(date) | uBAM Generation Finished." >> ${LOGFILE}
#---/ 
#---| BAM (Aligned BAM) GENERATION |
#---\
    echo "$(date) | Aligned BAM Generation using BWA Start." >> ${LOGFILE}
#---/ 
    /storage/apps/bin/bwa mem \
        -M -t 25 \
        -R "@RG\tID:${RGID}\tPL:${RGPL}\tLB:${RGLB}\tSM:${RGSM}\tCN:${RGCN}" \
        ${BWA_INDEX} \
        ${FASTQ_DIR}/${SEQ_ID}.trimmed_R1.fastq.gz \
        ${FASTQ_DIR}/${SEQ_ID}.trimmed_R2.fastq.gz > ${BAM_DIR}/${SEQ_ID}.bwa.mem.sam
#---\
    echo "$(date) | Aligned BAM Generation using BWA Finished." >> ${LOGFILE}
#---/ 
#---| BAM MERGE ALIGNMENT | 
#---\
    echo "$(date) | Merge Aligned-BAM and uBAM Start." >> ${LOGFILE}
#---/ 
    java -Xmx16G -jar /storage/apps/bin/picard.jar MergeBamAlignment \
        --REFERENCE_SEQUENCE ${GENOME_FASTA} \
        --UNMAPPED_BAM ${BAM_DIR}/${SEQ_ID}.fastqtosam.bam \
        --ALIGNED_BAM ${BAM_DIR}/${SEQ_ID}.bwa.mem.sam \
        --OUTPUT ${BAM_DIR}/${SEQ_ID}.merged.bam \
        --CREATE_INDEX true \
        --ADD_MATE_CIGAR true \
        --CLIP_ADAPTERS false \
        --CLIP_OVERLAPPING_READS true \
        --INCLUDE_SECONDARY_ALIGNMENTS true \
        --MAX_INSERTIONS_OR_DELETIONS -1 \
        --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
        --ATTRIBUTES_TO_RETAIN XS \
        --TMP_DIR ${TMP_DIR}
#---\
    echo "$(date) | Merge Aligned-BAM and uBAM Finished." >> ${LOGFILE}
#---/ 
#---| MARK-DUPLICATES |
#---\
    echo "$(date) | Mark Duplicates Start." >> ${LOGFILE}
#---/ 
    java -XX:+UseParallelGC -XX:ParallelGCThreads=10 -Xmx16g -jar /storage/apps/bin/gatk-local.jar MarkDuplicates \
        --INPUT ${BAM_DIR}/${SEQ_ID}.merged.bam \
        --METRICS_FILE ${QC_DIR}/${SEQ_ID}.merged.dup.marked.metrics.txt \
        --OUTPUT ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.bam
#---\
    echo "$(date) | Mark Duplicates Finished." >> ${LOGFILE}
#---/         
#---| BASE-RECALIBRATION |
#---\
    echo "$(date) | Base-Recalibration Start." >> ${LOGFILE}
#---/ 
    java -XX:+UseParallelGC -XX:ParallelGCThreads=10 -Xmx16g -jar /storage/apps/bin/gatk-local.jar BaseRecalibrator \
        --input ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.bam \
        --reference ${GENOME_FASTA} \
        --known-sites ${VCF_DBSNP} \
        --known-sites ${VCF_1KG_SNP} \
        --known-sites ${VCF_GNOMAD} \
        --known-sites ${VCF_INDEL} \
        --known-sites ${VCF_1KG_INDEL} \
        --output ${QC_DIR}/${SEQ_ID}.recal.table.txt
#---\
    echo "$(date) | base-recalibration-table generated. " >> ${LOGFILE}
#---/ 
    java -XX:+UseParallelGC -XX:ParallelGCThreads=10 -Xmx16g -jar /storage/apps/bin/gatk-local.jar ApplyBQSR \
        --input ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.bam \
        --output ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.recal.bam \
        --bqsr-recal-file ${QC_DIR}/${SEQ_ID}.recal.table.txt
#---\
    echo "$(date) | Base-Recalibration Finished." >> ${LOGFILE}
#---/ 
#---| LOCAL-REALIGNMENT |
#---\
    echo "$(date) | Local Realignment Start." >> ${LOGFILE}
    echo "          - run using java-8" >> ${LOGFILE}
#---/ 
    /usr/lib/jvm/java-1.8.0-openjdk-amd64/bin/java -Xmx16G -jar /storage/apps/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
        -T RealignerTargetCreator \
        -R ${GENOME_FASTA} \
        -L ${TARGET_BED} \
        -known ${VCF_INDEL} \
        -known ${VCF_1KG_INDEL} \
        -I ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.recal.bam \
        -o ${QC_DIR}/${SEQ_ID}.merged.dup.marked.recal.bam.realignertargetcreator.intervals
#---\
    echo "$(date) | realign-target generated." >> ${LOGFILE}
#---/
    /usr/lib/jvm/java-1.8.0-openjdk-amd64/bin/java -Xmx16G -jar /storage/apps/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
        -T IndelRealigner \
        -R ${GENOME_FASTA} \
        -targetIntervals ${QC_DIR}/${SEQ_ID}.merged.dup.marked.recal.bam.realignertargetcreator.intervals \
        -known ${VCF_INDEL} \
        -known ${VCF_1KG_INDEL} \
        -I ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.recal.bam \
        -o ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.recal.realign.bam  
#---\
    echo "$(date) | Local Realignment Finished." >> ${LOGFILE}
#---/ 
#---| BAM QC |
#---\
    echo "$(date) | QC-Process Start." >> ${LOGFILE}
    echo "$(date) |   QC-1. [ Sequecning Artifact ] Start." >> ${LOGFILE}
#---/ 
  #-| CollectSequencingArtifactMetrics : only final bam
    java -XX:+UseParallelGC -XX:ParallelGCThreads=10 -Xmx16g -jar /storage/apps/bin/gatk-local.jar CollectSequencingArtifactMetrics \
        --INPUT ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.recal.realign.bam  \
        --OUTPUT ${QC_DIR}/${SEQ_ID} --FILE_EXTENSION .txt \
        --REFERENCE_SEQUENCE ${GENOME_FASTA}
#---\
    echo "$(date) |   QC-1. [ Sequecning Artifact ] Finished." >> ${LOGFILE}
    echo "$(date) |   QC-2. [ Alignment Summary ] Start." >> ${LOGFILE}
#---/ 
  #-| CollectAlignmentSummaryMetrics : all bam
    BAM_SUFFIX="merged.bam merged.dup.marked.bam merged.dup.marked.recal.bam merged.dup.marked.recal.realign.bam"
    parallel -j 3 -k "java -XX:+UseParallelGC -XX:ParallelGCThreads=10 -Xmx16g -jar /storage/apps/bin/gatk-local.jar CollectAlignmentSummaryMetrics \
        --INPUT ${BAM_DIR}/${SEQ_ID}.{} \
        --OUTPUT ${QC_DIR}/${SEQ_ID}.{}.AlignSummary.metrics.txt" ::: ${BAM_SUFFIX}
#---\
    echo "$(date) |   QC-2. [ Alignment Summary ] Finished." >> ${LOGFILE}
    echo "$(date) |   QC-3. [ HS-Metric ] Start." >> ${LOGFILE}
#---/ 
  #-| CollectHsMetrics : all bam
    parallel -j 3 -k "java -XX:+UseParallelGC -XX:ParallelGCThreads=10 -Xmx16g -jar /storage/apps/bin/gatk-local.jar CollectHsMetrics \
        --INPUT ${BAM_DIR}/${SEQ_ID}.{} \
        --OUTPUT ${QC_DIR}/${SEQ_ID}.{}.HS.metrics.txt \
        --BAIT_INTERVALS ${BAIT_BED} \
        --TARGET_INTERVALS ${TARGET_BED}" ::: ${BAM_SUFFIX}  
#---\
    echo "$(date) |   QC-3. [ HS-Metric ] Finished." >> ${LOGFILE}
    echo "$(date) |   QC-4. [ Mosdepth & Alfred ] Start." >> ${LOGFILE}
#---/ 
  #-| mosdepth & alfred (optional) : only final bam
    source /storage/apps/miniconda3/bin/activate wes
    #
    /storage/apps/mosdepth-0.3.5/mosdepth --threads 15 \
        --by ${TARGET_BED} \
        --mapq 20 \
        --thresholds 1,10,30,50,100,150,200,250,300,400,500 \
        ${QC_DIR}/${SEQ_ID} \
        ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.recal.realign.bam 
    ##
    /storage/apps/alfred-0.2.7/bin/alfred qc \
        --reference ${GENOME_FASTA} \
        --bed ${TARGET_BED} \
        --outfile ${QC_DIR}/${SEQ_ID}.alfred.qc.tsv.gz \
        ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.recal.realign.bam 
    ##
    zgrep ^CM ${QC_DIR}/${SEQ_ID}.alfred.qc.tsv.gz | cut -f 2- > ${QC_DIR}/${SEQ_ID}.alfred.chr.map.stats.txt
    zgrep ^TC ${QC_DIR}/${SEQ_ID}.alfred.qc.tsv.gz | cut -f 2- > ${QC_DIR}/${SEQ_ID}.alfred.target.coverage.txt
    #
    conda deactivate
#---\
    echo "$(date) |   QC-4. [ Mosdepth & Alfred ] Finished." >> ${LOGFILE}
    echo "$(date) | QC-Process Finished." >> ${LOGFILE}
#---/ 

#---| TUMOR-ONLY VARIANT CALLING WITH MUTECT2
#---\
    echo "$(date) | Variant Calling Start." >> ${LOGFILE}
    echo "$(date) |   VC-1. [ MUTECT2 ] Start." >> ${LOGFILE}
#---/ 
  #-| Run Mutect2
    java -Xmx32G -jar /storage/apps/bin/gatk-local.jar Mutect2 \
        --reference ${GENOME_FASTA} \
        --intervals ${TARGET_BED} \
        --input ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.recal.realign.bam  \
        --germline-resource ${VCF_GNOMAD} \
        --panel-of-normals ${VCF_PON} \
        --f1r2-tar-gz ${QC_DIR}/${SEQ_ID}.f1r2.tar.gz \
        --output ${VCF_DIR}/${SEQ_ID}.mutect2.vcf
#---\
    echo "$(date) |   VC-1. [ MUTECT2 ] Finished." >> ${LOGFILE}
    echo "$(date) |   VC-2. [ Read-Orientation-Model ] Start." >> ${LOGFILE}
#---/ 
  #-| Orientation-Bias Filter 
    java -XX:+UseParallelGC -XX:ParallelGCThreads=10 -Xmx8G -jar /storage/apps/bin/gatk-local.jar LearnReadOrientationModel \
        -I ${QC_DIR}/${SEQ_ID}.f1r2.tar.gz \
        -O ${QC_DIR}/${SEQ_ID}.read-orientation-model.tar.gz
#---\
    echo "$(date) |   VC-2. [ Read-Orientation-Model ] Finished." >> ${LOGFILE}
    echo "$(date) |   VC-3. [ Pileup-Ggeneration ] Start." >> ${LOGFILE}
#---/ 
  #-| Generate Pileup Summaries on Tumor Sample
    java -XX:+UseParallelGC -XX:ParallelGCThreads=10 -Xmx8G -jar /storage/apps/bin/gatk-local.jar GetPileupSummaries \
        -I ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.recal.realign.bam \
        -O ${QC_DIR}/${SEQ_ID}.targeted_sequencing.table \
        -V ${VCF_GNOMAD} \
        -L ${TARGET_BED} \
        -R ${GENOME_FASTA}
#---\
    echo "$(date) |   VC-3. [ Pileup-Generation ] Finished." >> ${LOGFILE}
    echo "$(date) |   VC-4. [ Contamination-Calculation ] Start." >> ${LOGFILE}
#---/ 
  #-| Calculate Contamination on Tumor Sample
    java -XX:+UseParallelGC -XX:ParallelGCThreads=10 -Xmx16G -jar /storage/apps/bin/gatk-local.jar CalculateContamination \
        -I ${QC_DIR}/${SEQ_ID}.targeted_sequencing.table \
        -O ${QC_DIR}/${SEQ_ID}.targeted_sequencing.contamination.table
#---\
    echo "$(date) |   VC-4. [ Contamination-Calculation ] Finished." >> ${LOGFILE}
    echo "$(date) |   VC-5. [ FILTER-MUTECT2 ] Start." >> ${LOGFILE}
#---/ 
  #-| Run FilterMutectCalls
    java -XX:+UseParallelGC -XX:ParallelGCThreads=10 -Xmx16G -jar /storage/apps/bin/gatk-local.jar FilterMutectCalls \
        -V ${VCF_DIR}/${SEQ_ID}.mutect2.vcf \
        -L ${TARGET_BED} \
        --reference ${GENOME_FASTA} \
        --contamination-table ${QC_DIR}/${SEQ_ID}.targeted_sequencing.contamination.table \
        --ob-priors ${QC_DIR}/${SEQ_ID}.read-orientation-model.tar.gz \
        -O ${VCF_DIR}/${SEQ_ID}.mutect2.filtered.vcf
#---\
    echo "$(date) |   VC-5. [ FILTER-MUTECT2 ] Finished." >> ${LOGFILE}
    echo "$(date) | Variant Calling Finished." >> ${LOGFILE}
#---/ 

#---| VARIANT ANNOTATION : VEP
#---\
    echo "$(date) | Variant Annotation Start." >> ${LOGFILE}
    echo "$(date) |   VA-1. [ Ensembl-VEP ] Start." >> ${LOGFILE}
#---/ 
    CACHE_DIR=/storage/references_and_index/ensembl-vep
  ##
    singularity exec -B $CACHE_DIR:/opt/vep/.vep,$DATA_DIR:/vep-input /storage/images/vep-110.sif /opt/vep/src/ensembl-vep/vep \
        --cache --dir_cache /opt/vep/.vep --force_overwrite --stats_text --species homo_sapiens --assembly GRCh37 \
        --fork 30 --show_ref_allele --uploaded_allele --numbers \
        --everything \
        --plugin AlphaMissense,file=/opt/vep/.vep/alpha-missense-prediction/AlphaMissense_b37.tsv.gz \
        --input_file /vep-input/${SEQ_ID}/vcf/${SEQ_ID}.mutect2.filtered.vcf \
        --vcf --output_file /vep-input/${SEQ_ID}/vcf/${SEQ_ID}.mutect2.filtered.vep.annotate.vcf
#---\
    echo "$(date) |   VA-1. [ Ensembl-VEP ] Finished." >> ${LOGFILE}
    echo "$(date) |   VA-2. [ VCF-to-MAF ] Start." >> ${LOGFILE}
#---/  
    perl /storage/apps/vcf2maf-1.6.21/vcf2maf.pl \
        --inhibit-vep \
        --input-vcf ${VCF_DIR}/${SEQ_ID}.mutect2.filtered.vep.annotate.vcf \
        --output-maf ${VCF_DIR}/${SEQ_ID}.mutect2.filtered.vep.annotate.maf \
        --tumor-id ${SEQ_ID} \
        --retain-info DP,FILTER \
        --retain-ann am_class,am_pathogenicity
#---\
    echo "$(date) |   VA-2. [ VCF-to-MAF ] Finished." >> ${LOGFILE}
    echo "$(date) | Variant Annotation Finished." >> ${LOGFILE}
#---/  

#---| QC RESULT SUMMARY AND IMPORT INTO DB
#---\
    echo "$(date) | QC-Rsults Summarization Start." >> ${LOGFILE}
#---/  
    SCRIPT_DIR=/storage/home/gcx/run_scripts
    source /storage/apps/miniconda3/bin/activate wes
    # singularity exec -B ${QC_DIR}:/multiqc,${SCRIPT_DIR}:/multiqc-config /storage/images/multiqc --force \
    #     --filename ${SEQ_ID}_qc \
    #     --outdir /multiqc/qc_res \
    #     --config /multiqc-config/gcx_wes_v1.2_multiqc_custom_config.yaml \
    #     --data-dir /multiqc
    multiqc --force \
        --filename ${SEQ_ID}_qc \
        --outdir ${QC_DIR}/qc_res \
        --config ${SCRIPT_DIR}/gcx_wes_v1.2_multiqc_custom_config.yaml \
        --data-dir ${QC_DIR}
    conda deactivate  
#---\
    echo "$(date) | QC-Rsults Summarization Finished." >> ${LOGFILE}
    echo "$(date) | Summarized Result Import into DB Start." >> ${LOGFILE}
#---/  
  ## Import Into DB  
    Rscript ${SCRIPT_DIR}/gcx_wes_v1.2_multiqc_post_processing.R --seqFolder ${SEQ_FOLDER} --seqID ${SEQ_ID}
#---\
    echo "$(date) | Summarized Result Import into DB Finished." >> ${LOGFILE}
    echo "$(date) | raw MAF data Import into DB Start." >> ${LOGFILE}
#---/  
#---| MAF IMPORT TO DB
    Rscript ${SCRIPT_DIR}/gcx_wes_v1.2_raw_maf_db_import.R --seqFolder ${SEQ_FOLDER} --seqID ${SEQ_ID}
#---\
    echo "$(date) | raw MAF data Import into DB Finished." >> ${LOGFILE}
#---/
#---\
    echo "$(date) | ALL PROCESSING FINISHED." >> ${LOGFILE}
#---/

