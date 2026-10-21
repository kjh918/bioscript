
#!/bin/bash

#---/ Input Options /
    usage() {
        echo "Usage: GCX-WES-Processing.sh [ options ]
        [ -d | --dataDir ]        : data processing base directory
        [ -s | --sampleInfoFile ] : sample info file (sample-id file)
        [ -r | --referenceVer ]   : Reference Genome Version (default = hg19)
        [ -l | --logFile ]        : logFile name (file will be created or overwritten)
        [ -j | --jobs ]           : jobs for parallel processing (default = the numbers of samples)  
        [ -h | --help ]           : Print this usage 
        [ --RGID ]                : ReadGroup ID (default = GCX)
        [ --RGLB ]                : ReadGroup Library (default = TwistExome2)
        [ --RGPL ]                : ReadGroup Platform (default = ILLUMINA)
        [ --runFastQC ]           : perform factQC (default = true)
        [ --runFastp ]            : perform fastp (default = true)
        [ --runBwa ]              : perform BAW alignment (default = true)
        [ --runBamProcessing ]    : perform BAM Post-processing (default = true)
        [ --runBamQC ]            : perform BAM-QC (default = true)        
        "
    }

  #---/ Set Defualt Values /
    RGID="GCX" && RGLB="TwistExome2" && RGPL="ILLUMINA"
    REF_VER="hg19"
    runFASTQC="yes" && runFASTP="yes" && runBWA="yes" && runBAMPROC="yes" && runBAMQC="yes"
    JOBS=1
  #---/
#---/

#---/ Parse Option Arguments /
    ARGS=$(getopt -a -o d:s:r:l:j:h: --long dataDir:,sampleInfoFile:,RGID:,RGLB:,RGPL:,referenceVer:,logFile:,jobs:,runFastQC:,runFastp:,runBwa:,runBamProcessing:,runBamQC:,help -- "$@" )
    VALID_ARGS=$?
    if [ "$VALID_ARGS" != "0" ]; then 
        usage >&2 
        exit 2
    fi

    eval set -- "$ARGS"
    while :
    do
    case "$1" in
        -d | --dataDir )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  DATA_DIR=$2 ; shift 2 ;;
            esac ;;
        -s | --sampleInfoFile )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  SAMPLE_ID_FILE=$2 ; shift 2 ;;
            esac  ;;
        --RGID )
            case "$2" in
                -* | --* | "") echo "no value at $1, Default ID = GCX" ; shift ;;
                *)  RGID=$2 ; shift 2 ;;
            esac ;;
        --RGLB )
            case "$2" in
                -* | --* | "") echo "no value at $1, Default Library = TwistExome2" ; shift ;;
                *)  RGLB=$2 ; shift 2 ;;
            esac ;;
        --RGPL )
            case "$2" in
                -* | --* | "") echo "no value at $1, Default Platform = ILLUMINA" ; shift ;;
                *)  RGPL=$2 ; shift 2 ;;
            esac ;;
        -r | --referenceVer )
            case "$2" in
                -* | --* | "") echo "no value at $1, Default Reference Version = hg19" ; shift ;;
                *)  REF_VER=$2 ; shift 2 ;;
            esac ;;
        -l | --logFile )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  LOG_FILE=$2 ; shift 2 ;;
            esac ;;
        -j | --jobs )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  JOBS=$2 ; shift 2 ;;
            esac ;;
        --runFastQC )
            case "$2" in
                -* | --* | "") echo "no value at $1, Default = yes" ; shift ;;
                *)  runFASTQC=$2 ; shift 2 ;;
            esac ;;
        --runFastp )
            case "$2" in
                -* | --* | "") echo "no value at $1, Default = yes" ; shift ;;
                *)  runFASTP=$2 ; shift 2 ;;
            esac ;;
        --runBwa )
            case "$2" in
                -* | --* | "") echo "no value at $1, Default = no" ; shift ;;
                *)  runBWA=$2 ; shift 2 ;;
            esac ;;
        --runBamProcessing )
            case "$2" in
                -* | --* | "") echo "no value at $1, Default = no" ; shift ;;
                *)  runBAMPROC=$2 ; shift 2 ;;
            esac ;;
        --runBamQC )
            case "$2" in
                -* | --* | "") echo "no value at $1, Default = no" ; shift ;;
                *)  runBAMQC=$2 ; shift 2 ;;
            esac ;;
        -h | --help )
        usage >&2 ; exit 2 ;;   
        --) shift ; break ;;
        *)  usage >&2 ; exit 2 ;;
    esac
    done
#---/

#---/ LOG FILE CREATE
    LOG_FILE=${DATA_DIR}/WES-processig.log
    if [ -f ${LOG_FILE} ]; then 
        mv ${LOG_FILE} ${LOG_FILE}.prev
        touch ${LOG_FILE}
    else
        touch ${LOG_FILE}
    fi 
    echo "---------------------------------------------------------------------" >> ${LOG_FILE}
    echo "### WES Processing Log ###                                           " >> ${LOG_FILE}
    echo "                                                                     " >> ${LOG_FILE}
    echo "$(date) : Processing START                                        " >> ${LOG_FILE}                                                
#---/

#---/ DIRECTORIES 
  #---/ FASTQ folders
    MERGED_FASTQ_DIR=${DATA_DIR}/fastq_merged
    TRIM_FASTQ_DIR=${DATA_DIR}/fastq_trimmed
    if [ ! -d ${TRIM_FASTQ_DIR} ]; then 
        mkdir -p ${TRIM_FASTQ_DIR}
        echo "Trimmed-Fastq-Folder not existed, newly created now."
    fi
  #---/  
  #---/ BAM folders
    BAM_DIR=${DATA_DIR}/bam_files
    if [ ! -d ${BAM_DIR} ]; then 
        mkdir -p ${BAM_DIR}
        echo "BAM-File-Folder not existed, newly created now."
    fi
  #---/ 
  #---/ VCF Result folders
    VCF_RES_DIR=${DATA_DIR}/vcf_files
    if [ ! -d ${VCF_RES_DIR} ]; then 
        mkdir -p ${VCF_RES_DIR}
        echo "VCF-Result-Folder not existed, newly created now."
    fi
  #---/ 
  #---/ QC-Result folders
    QC_RES_DIR=${DATA_DIR}/qc_files
    if [ ! -d ${QC_RES_DIR} ]; then
        mkdir -p ${QC_RES_DIR}
        mkdir -p ${QC_RES_DIR}/reports
        mkdir -p ${QC_RES_DIR}/bam_stats
        echo "QC-Result-File-Folder not existed, newly created now."
    fi
  #---/ 
  #---/ FIANL-RESULT folder
    FINAL_RES_DIR=${DATA_DIR}/Results
    if [ ! -d ${FINAL_RES_DIR} ]; then 
        mkdir -p ${FINAL_RES_DIR}
        echo "Final-Result-Folder not existed, newly created now."
    fi
  #---/ 
    echo "-- Data Folders --                                                   " >> ${LOG_FILE}                                                
    echo " Initial Fastq Files : ${MERGED_FASTQ_DIR}                           " >> ${LOG_FILE}
    echo " Trimmed Fastq Files : ${TRIM_FASTQ_DIR}                             " >> ${LOG_FILE}
    echo "           BAM Files : ${BAM_DIR}                                    " >> ${LOG_FILE}
    echo "           VCF Files : ${VCF_RES_DIR}                                " >> ${LOG_FILE}
    echo "    QC-Results Files : ${QC_RES_DIR}                                 " >> ${LOG_FILE}
    echo "  Final Result Files : ${FINAL_RES_DIR}                              " >> ${LOG_FILE}
    echo "                                                                     " >> ${LOG_FILE}
#---/

#---/ REFERENCES
    if [ "${REF_VER}" -eq "hg38" ]; then
        REF_DIR="/BiO/Reference_and_Index/gcx.hg38"
    else 
        REF_DIR="/BiO/Reference_and_Index/gcx.hg19"
    fi
  #---/ Files
    REF_DIR=/BiO/Reference_and_Index/gcx.hg19
    BWA_INDEX=${REF_DIR}/gcx.hg19
    REF_GENOME_FASTA=${REF_DIR}/gcx.hg19.fasta
    BAIT_BED=${REF_DIR}/TwistExome2.0/gcx.hg19.twist.exome2.probe.bed
    BAIT_INTERVAL=${REF_DIR}/TwistExome2.0/gcx.hg19.twist.exome2.probe.interval_list
    TARGET_BED=${REF_DIR}/TwistExome2.0/gcx.hg19.twist.exome2.target.bed
    TARGET_INTERVAL=${REF_DIR}/TwistExome2.0/gcx.hg19.twist.exome2.target.interval_list
    VCF_DBSNP=${REF_DIR}/VCF/gcx.hg19.dbsnp138.vcf.gz 
    VCF_1000GSNP=${REF_DIR}/VCF/gcx.hg19.1000G_phase1.snps.high_confidence.vcf.gz
    VCF_1000GOMNI=${REF_DIR}/VCF/gcx.hg19.1000G_omni2.5.vcf.gz 
    VCF_HAPMAP=${REF_DIR}/VCF/gcx.hg19.hapmap_3.3.vcf.gz 
    VCF_INDEL=${REF_DIR}/VCF/gcx.hg19.known_indels.vcf.gz
    VCF_1000GINDEL=${REF_DIR}/VCF/gcx.hg19.Mills_and_1000G_gold_standard.indels.sites.vcf.gz 
  #---\
    echo "-- Reference --                                                      " >> ${LOG_FILE}                                                
    echo " Reference Version : ${REF_VER}                                      " >> ${LOG_FILE}
    echo "                                                                     " >> ${LOG_FILE}
#---/

#---/ READ SAMPLE ID LIST 
    SAMPLE_LIST=$(cat ${SAMPLE_ID_FILE})
    SAMPLE_NUMBERS=$(echo ${SAMPLE_LIST} | wc -w)
    if [[ "${JOBS}" -eq "1"  &&  "${SAMPLE_NUMBERS}" -eq "1" ]]; then
        JOBS=${JOBS}
    else
        JOBS=${SAMPLE_NUMBERS}
    fi
  #---\
    echo "-- Samples --                                                        " >> ${LOG_FILE}                                                
    echo "    Total Sample Numbers : ${SAMPLE_NUMBERS}                         " >> ${LOG_FILE} 
    echo "Parallel Processing Jobs : ${JOBS}                                   " >> ${LOG_FILE}
    echo "                                                                     " >> ${LOG_FILE}
    echo "---------------------------------------------------------------------" >> ${LOG_FILE}
    echo "                                                                     " >> ${LOG_FILE}
  #---/ 
#---/

#---/ FASTQC
    if [ "${runFASTQC}" == "yes" ]; then
      #---\
        echo "$(date) : FastQC START.                                          " >> ${LOG_FILE}
      #---/  
        FASTQC_THREADS=`expr ${SAMPLE_NUMBERS} * 2`
        /BiO/SWTools/fastqc --noextract -t ${FASTQC_THREADS} -o ${QC_RES_DIR} ${MERGED_FASTQ_DIR}/*.fastq.gz
      #---\
        echo "$(date) : FastQC FINISHED.                                       " >> ${LOG_FILE}
      #---/  
    else
        echo "$(date) : FastQC is Skipped by User.                             " >> ${LOG_FILE}
    fi
#---/

#---\ FASTQ-TRIMMING & QUALITY CUT-OFF
    if [ "${runFASTP}" == "yes" ]; then
      #---\
        echo "$(date) : fastq Trimming and Filtering START.                    " >> ${LOG_FILE}
        echo "        Quanilty Cut-off = 20                                    " >> ${LOG_FILE}
        echo "          Length Cut-off = 105                                   " >> ${LOG_FILE}
      #---/  
        parallel -j ${JOBS} -k "/BiO/SWTools/fastp \
            --in1 ${MERGED_FASTQ_DIR}/{}_R1.fastq.gz \
            --in2 ${MERGED_FASTQ_DIR}/{}_R2.fastq.gz \
            --html ${QC_RES_DIR}/{}.fastp.html \
            --json ${QC_RES_DIR}/{}.fastp.json \
            --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
            --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            --length_required 105 \
            --qualified_quality_phred 20 \
            --out1 ${TRIM_FASTQ_DIR}/{}.trimmed_R1.fastq.gz \
            --out2 ${TRIM_FASTQ_DIR}/{}.trimmed_R2.fastq.gz" ::: ${SAMPLE_LIST}
      #---\
        echo "$(date) : fastq Trimming and Filtering FINISHED.                 " >> ${LOG_FILE}
      #---/
    else
        echo "$(date) : fastq Trimming and Filtering is Skipped by User.       " >> ${LOG_FILE}
    fi
#---/

#---\ BWA ALIGNMENT
    if [ "${runBWA}" == "yes" ]; then
      #---\
        echo "$(date) : BWA Alignment START.                                   " >> ${LOG_FILE}
      #---/
        parallel -j ${JOBS} -k "/BiO/SWTools/bwa mem -aM \
            -R '@RG\tID:${RGID}\tLB:${RGLB}\tSM:{}\tPL:${RGPL}' \
            -t 20 \
            -k 70 \
            ${BWA_INDEX} \
            ${TRIM_FASTQ_DIR}/{}.trimmed_R1.fastq.gz \
            ${TRIM_FASTQ_DIR}/{}.trimmed_R2.fastq.gz |
            /BiO/SWTools/samtools view -Shb -o ${BAM_DIR}/{}.bam - " ::: ${SAMPLE_LIST}
      #---\
        echo "$(date) : BWA Alignment FINISHED.                                " >> ${LOG_FILE}     
      #---/
    else
      #---\
        echo "$(date) : BWA Alignment is Skipped by User.                      " >> ${LOG_FILE}
      #---/  
#---/        

#---/ BAM POST PROCESSING
    if [ "${runBAMPROC}" == "yes" ]; then
        #---\   
        echo "$(date) : BAM Post-Processing START.                             " >> ${LOG_FILE}
        #---/  
      #---| SORT & INDEX |
        #---\   
        echo "$(date) : BAM  Sort....started                                   " >> ${LOG_FILE}
        #---/  
        parallel -j ${JOBS} -k "/BiO/SWTools/samtools sort \
            --threads 20 \
            -o ${BAM_DIR}/{}.sorted.bam \
            ${BAM_DIR}/{}.bam" ::: ${SAMPLE_LIST}
        #---\      
        echo "$(date) : BAM  Sort....done.                                     " >> ${LOG_FILE}  
        echo "$(date) : BAM Index....started                                   " >> ${LOG_FILE}
        #---/  
        parallel -j ${SAMPLE_NUMBERS} -k "/BiO/SWTools/samtools index \
            --threads 20 \
            ${BAM_DIR}/{}.sorted.bam" ::: ${SAMPLE_LIST}
        #---\      
        echo "$(date) : BAM Index....done                                      " >> ${LOG_FILE}
        #---/
      #---|    
      #---| MARK-DUPLICATES |
        #---\ 
        echo "$(date) : BAM MarkDuplicates....started                          " >> ${LOG_FILE}
        #---/  
        parallel -j ${JOBS} -k "java -XX:+UseParallelGC -XX:ParallelGCThreads=20 -Xmx60g -jar /BiO/SWTools/gatk-local.jar MarkDuplicates \
            --INPUT ${BAM_DIR}/{}.sorted.bam \
            --METRICS_FILE ${QC_RES_DIR}/{}.duplicates.marked.metrics.txt \
            --OUTPUT ${BAM_DIR}/{}.sorted.dup.marked.bam" ::: ${SAMPLE_LIST} 
        #---\      
        echo "$(date) : BAM MarkDuplicates....Done.                            " >> ${LOG_FILE}
        #---/
      #---|  
      #---| BASE RECALIBRATION |
        #---\  
        echo "$(date) : BAM BaseRecalibration....started                       " >> ${LOG_FILE}
        #---/  
        parallel -j ${JOBS} -k  "java -XX:+UseParallelGC -XX:ParallelGCThreads=10 -Xmx60g -jar /BiO/SWTools/gatk-local.jar BaseRecalibrator \
            --input ${BAM_DIR}/{}.sorted.dup.marked.bam \
            --reference ${REF_GENOME_FASTA} \
            --known-sites ${VCF_DBSNP} \ 
            --known-sites ${VCF_1000GSNP} \
            --known-sites ${VCF_1000GOMNI} \
            --known-sites ${VCF_HAPMAP} \
            --known-sites ${VCF_INDEL} \
            --known-sites ${VCF_1000GINDEL} \
            --output ${BAM_DIR}/{}.recal.table.txt " ::: ${SAMPLE_LIST} 
        #
        parallel -j ${JOBS} -k "java -XX:+UseParallelGC -XX:ParallelGCThreads=10 -Xmx60g -jar /BiO/SWTools/gatk-local.jar ApplyBQSR \
            --input ${BAM_DIR}/{}.sorted.dup.marked.bam \
            --output ${BAM_DIR}/{}.sorted.dup.marked.recal.bam \
            --bqsr-recal-file ${BAM_DIR}/{}.recal.table.txt" ::: ${SAMPLE_LIST}
        #---\
        echo "$(date) : BAM BaseRecalibration....Done.                         " >> ${LOG_FILE}      
        #---/
      #---|  
        #---\
        echo "$(date) : BAM Post-Processing FINISHED.                          " >> ${LOG_FILE}
        #---/
    else
        echo "$(date) : BAM Post-Processing is Skipped by User.                " >> ${LOG_FILE}   
    fi 
#---/       

#---/ BAM PROCESSING QC & STATS 
    if [ "${runBAMQC}" == "yes" ]; then
        #---\   
        echo "$(date) : BAM-QC Analysis START.                                 " >> ${LOG_FILE}
        #---/  
        BAM_SUFFIX="sorted.bam sorted.dup.marked.bam sorted.dup.marked.recal.bam"
      #---| CollectAlignmentSummaryMetrics |
        #---\
        echo "$(date) : BAM-QC : alignment summary metrics....started          " >> ${LOG_FILE}
        #---/      
        for bamSuffix in ${BAM_SUFFIX}
        do
            parallel -j ${JOBS} -k "java -XX:+UseParallelGC -XX:ParallelGCThreads=20 -Xmx60g -jar /BiO/SWTools/gatk-local.jar CollectAlignmentSummaryMetrics \
                --INPUT ${BAM_DIR}/{}.${bamSuffix} \
                --ADAPTER_SEQUENCE AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
                --ADAPTER_SEQUENCE AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
                --OUTPUT ${QC_RES_DIR}/{}.${bamSuffix}.AlignSummary.metrics.txt" ::: ${SAMPLE_LIST}   
        done
        #---\
        echo "$(date) : BAM-QC : alignment summary metrics....Done.            " >> ${LOG_FILE}      
        #---/
      #---| CollectHsMetrics |
        #---\
        echo "$(date) : BAM-QC : collect hybrid-selection metrics....started" >> ${LOG_FILE} 
        #---/
        for bamSuffix in ${BAM_SUFFIX}
        do
            parallel -j ${JOBS} -k "java -XX:+UseParallelGC -XX:ParallelGCThreads=20 -Xmx60g -jar /BiO/SWTools/gatk-local.jar CollectHsMetrics \
                --INPUT ${BAM_DIR}/{}.${bamSuffix} \
                --OUTPUT ${QC_RES_DIR}/{}.${bamSuffix}.HS.metrics.txt \
                --BAIT_INTERVALS ${BAIT_INTERVAL} \
                --TARGET_INTERVALS ${TARGET_INTERVAL} " ::: ${SAMPLE_LIST} 
        done
        #---\
        echo "$(date) : BAM-QC : collect hybrid-selection metrics....Done.     " >> ${LOG_FILE}       
        #---/
      #---| CollectInsertSizeMetrics |
        #---\
        echo "$(date) : BAM-QC : collect insert-size metrics....started" >> ${LOG_FILE} 
        #---/
        for bamSuffix in ${BAM_SUFFIX}
        do
            parallel -j ${JOBS} -k "java -XX:+UseParallelGC -XX:ParallelGCThreads=20 -Xmx60g -jar /BiO/SWTools/gatk-local.jar CollectInsertSizeMetrics \
                --INPUT ${BAM_DIR}/{}.${bamSuffix} \
                --OUTPUT ${QC_RES_DIR}/{}.${bamSuffix}.insertsize.metrics.txt \
                -H ${QC_RES_DIR}/{}.${bamSuffix}.insertsize.histogram.pdf " ::: ${SAMPLE_LIST} 
        done
        #---\
        echo "$(date) : BAM-QC : collect insert-size metrics....Done.          " >> ${LOG_FILE}       
        #---/
      #---| CollectOxoGMetrics |
        #---\
        echo "$(date) : BAM-QC : collect assess oxidative artifacts metrics....started" >> ${LOG_FILE} 
        #---/
        for bamSuffix in ${BAM_SUFFIX}
        do
            parallel -j ${JOBS} -k "java -XX:+UseParallelGC -XX:ParallelGCThreads=20 -Xmx60g -jar /BiO/SWTools/gatk-local.jar CollectOxoGMetrics \
                --INPUT ${BAM_DIR}/{}.${bamSuffix} \
                --OUTPUT ${QC_RES_DIR}/{}.${bamSuffix}.oxoG.metrics.txt \
                --REFERENCE_SEQUENCE ${REF_GENOME_FASTA} " ::: ${SAMPLE_LIST} 
        done
        #---\
        echo "$(date) : BAM-QC : collect assess oxidative artifacts metrics....Done.  " >> ${LOG_FILE}       
        #---/
      #---| mosdepth
        #---\
        echo "$(date) : BAM-QC : run mosdepth....started" >> ${LOG_FILE} 
        #---/
        source /BiO/Install/miniconda3/bin/activate wes
        parallel -j ${JOBS} -k "mosdepth --threads 10 \
            --by /BiO/Reference_and_Index/gcx.hg19/TwistExome2.0/gcx.hg19.twist.exome2.target.bed \
            ${QC_RES_DIR}/{} ${BAM_DIR}/{}.sorted.dup.marked.recal.bam" ::: ${SAMPLE_LIST}
        #---\
        echo "$(date) : BAM-QC : run mosdepth....DONE." >> ${LOG_FILE} 
        #---/
        #---\   
        echo "$(date) : BAM-QC Analysis FINISHED.                              " >> ${LOG_FILE}
        #---/  
    else
        echo "$(date) : BAM-QC Analysis is Skipped by User.                    " >> ${LOG_FILE}
    fi
#---/

#---/ WES QC-REPORT
  #---\
    echo "                                                                     " >> ${LOG_FILE}                                                
    echo "---------------------------------------------------------------------" >> ${LOG_FILE}
    echo "                                                                     " >> ${LOG_FILE}                                                
    echo " QC-Report & DB Import                                               " >> ${LOG_FILE}
    echo "                                                                     " >> ${LOG_FILE}
  #---/  
    MQC_REPORT_DIR=${QC_RES_DIR}/reports
    MQC_BAM_DIR=${QC_RES_DIR}/bam_stats
    BAM_SUFFIX2="sorted.bam sorted.dup.marked.bam sorted.dup.marked.recal.bam"
    for sampleID in ${SAMPLE_LIST}
    do
        for bamSuffix2 in ${BAM_SUFFIX2}
        do 
            cp ${QC_RES_DIR}/${sampleID}.${bamSuffix2}.AlignSummary.metrics.txt ${MQC_BAM_DIR}/
            cp ${QC_RES_DIR}/${sampleID}.${bamSuffix2}.HS.metrics.txt ${MQC_BAM_DIR}/
            cp ${QC_RES_DIR}/${sampleID}.${bamSuffix2}.insertsize.metrics.txt ${MQC_BAM_DIR}/
            cp ${QC_RES_DIR}/${sampleID}.${bamSuffix2}.oxoG.metrics.txt ${MQC_BAM_DIR}/
        done 
        cp ${QC_RES_DIR}/${sampleID}.duplicates.marked.metrics.txt ${MQC_BAM_DIR}/
      ###
        cp ${QC_RES_DIR}/${sampleID}.fastp.json ${MQC_REPORT_DIR}/
        cp ${QC_RES_DIR}/${sampleID}.duplicates.marked.metrics.txt ${MQC_REPORT_DIR}/  
        cp ${QC_RES_DIR}/${sampleID}.sorted.dup.marked.recal.bam.AlignSummary.metrics.txt ${MQC_REPORT_DIR}/  
        cp ${QC_RES_DIR}/${sampleID}.sorted.dup.marked.recal.bam.HS.metrics.txt ${MQC_REPORT_DIR}/  
        cp ${QC_RES_DIR}/${sampleID}.sorted.dup.marked.recal.bam.insertsize.metrics.txt ${MQC_REPORT_DIR}/  
        cp ${QC_RES_DIR}/${sampleID}.sorted.dup.marked.recal.bam.oxoG.metrics.txt ${MQC_REPORT_DIR}/
    done
  #---/ conda environment activation
    source /BiO/Install/miniconda3/bin/activate mqc
  #---/ Run MultiQC for WES-QC-REPORT
    multiqc --force \
        --filename qc_report \
        --outdir ${QC_RES_DIR} \
        --data-dir ${MQC_REPORT_DIR}   
  
  #---/ Run MultiQC for BAM-QC 
    multiqc --force \
        --filename qc_bam \
        --outdir ${QC_RES_DIR} \
        --config /BiO/SWTools/runScripts/multiqc-gcx-custom-config.yaml \
        --data-dir ${MQC_BAM_DIR}   
  #---/   

  #---/
    Rscript WES_qc_report.R --sampleListFile ${SAMPLE_ID_FILE} --dataDir ${DATA_DIR}
  #
  








