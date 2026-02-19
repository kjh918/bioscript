
#! /bin/bash

#---| SCRIPT USAGE |-------------------------------------------------------------------------------#
    usage() {
        echo "Usage: WTS_P04-1_BAM-QC.sh [ARGUMENTS]...

            --baseDir               Base Work Dir

            --batchID               Data processing base directory

            --seqID                 Sample ID (Used as sample's processing product folder name)

            --threads               cpu cores 

            --stranded              stranded seq or not 'true' or 'false', default= 'true'

            --genomeAseembly        Reference genome assembly. default = 'hg38'

            --paired                paired-end sequencing or not. default = 'true'

            --runMultiqc            run multiqc for qc results summarization and integration. default = 'true'

            --qcResDbImport         qc results database import. default = 'true' 

            --xenograftPipeline     run as Xenograft pipeline or not. default = 'false'  

            -h | --help             Print this usage 
        "
    }
#--------------------------------------------------------------------------------------------------#

#---| OPTION PARSE |-------------------------------------------------------------------------------#
    ARGS=$(getopt -a -o h: --long batchID:,seqID:,threads:,baseDir:,stranded:,genomeAseembly:,paired:,runMultiqc:,qcResDbImport:,xenograftPipeline:,help -- "$@" )
    VALID_ARGS=$?
    if [ "$VALID_ARGS" != "0" ]; then 
        usage >&2 
        exit 2
    fi

    eval set -- "$ARGS"
    while :
    do
    case "$1" in
        --baseDir )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = /data/wts" ; exit 2 ;;
                *)  BaseDir=$2 ; shift 2 ;;
            esac ;;
        --batchID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  BatchID=$2 ; shift 2 ;;
            esac ;;
        --seqID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  SeqID=$2 ; shift 2 ;;
            esac  ;;
        --threads )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 15" ; exit 2 ;;
                *)  Threads=$2 ; shift 2 ;;
            esac ;;
        --stranded )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = false" ; exit 2 ;;
                *)  Stranded=$2 ; shift 2 ;;
            esac ;; 
        --genomeAseembly )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = hg19" ; exit 2 ;;
                *)  GenomeAseembly=$2 ; shift 2 ;;
            esac ;;
        --paired )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = true" ; exit 2 ;;
                *)  Paired=$2 ; shift 2 ;;
            esac ;; 
        --runMultiqc )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = true" ; exit 2 ;;
                *)  RunMultiqc=$2 ; shift 2 ;;
            esac ;;
        --qcResDbImport )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = true" ; exit 2 ;;
                *)  QcResDbImport=$2 ; shift 2 ;;
            esac ;;
        --xenograftPipeline )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = true" ; exit 2 ;;
                *)  XenograftPipeline=$2 ; shift 2 ;;
            esac ;;
        -h | --help )
        usage >&2 ; exit 2 ;;   
        --) shift ; break ;;
        *)  usage >&2 ; exit 2 ;;
    esac
    done
#--------------------------------------------------------------------------------------------------#

#---| Manual Params |------------------------------------------------------------------------------#
    # BaseDir="/data/wts"
    # BatchID="WTS_25_02"
    # SeqID="WTS_25_02_03"
    # Threads=15
    # Stranded="true"
    # GenomeAseembly="hg38"
    # Paired="true"
    # RunMultiqc="ture"
    # QcResDbImport="true"
    # XenograftPipeline="true"
#--------------------------------------------------------------------------------------------------#


#---| FOLDERS |------------------------------------------------------------------------------------#
    DataDir=${BaseDir}/${BatchID}
    BamDir=${DataDir}/${SeqID}/bam
    QcDir=${DataDir}/${SeqID}/qcfiles
    TmpDir=${DataDir}/${SeqID}/tmp
    LogDir=${DataDir}/${SeqID}/log
    #--------------------------------------------------------------------------- 
    if [ ! -d ${TmpDir} ]; then mkdir -p ${TmpDir}; fi
    if [ ! -d ${LogDir} ]; then mkdir -p ${LogDir}; fi
    if [ ! -d ${QcDir} ]; then mkdir -p ${QcDir}/qc_res; fi
#--------------------------------------------------------------------------------------------------#

#---| REFERENCES |---------------------------------------------------------------------------------#
    GenomeFasta=/storage/references_and_index/${GenomeAseembly}/fasta/rnaseq/${GenomeAseembly}.fa
    GeneModelGtfFile=/storage/references_and_index/${GenomeAseembly}/gtf/${GenomeAseembly}.stranded.gene.model.gtf
    RefFlatFile=/storage/references_and_index/${GenomeAseembly}/gtf/${GenomeAseembly}.refFlat.txt
    RrnaIntervalList=/storage/references_and_index/${GenomeAseembly}/gtf/${GenomeAseembly}.rRNA.interval_list
#--------------------------------------------------------------------------------------------------#

#---| SWTOOLS RUN |--------------------------------------------------------------------------------#
    RunSingularity="singularity exec --env LC_ALL=C.UTF-8 --env LANG=C.UTF-8 -B /storage,/data /storage/images"
    RunGatk4="${RunSingularity}/gatk-4.4.0.0.sif gatk"
    RunRnaseqc="${RunSingularity}/rna-seqc-2.4.2.sif rnaseqc"
    RunScript="/storage/home/kangsm/runScripts/WTS"
#--------------------------------------------------------------------------------------------------#

#---| GATK COLLECT-RNASEQ-METRICS |----------------------------------------------------------------#
    # stranded option 
    if [ "${Stranded}" = "true" ]; then GatkStrandOption="SECOND_READ_TRANSCRIPTION_STRAND"; else GatkStrandOption="NONE"; fi
    #---------------------------------------------------------------------------

    # bam file by pipeline
    if [ "${XenograftPipeline}" = "true" ]; then

        # Xenograft piepline

        HumanBam=${BamDir}/${SeqID}.human.xenofilter.bam
        MouseBam=${BamDir}/${SeqID}.mouse.xenofilter.bam

        ${RunGatk4} CollectRnaSeqMetrics \
            --java-options "-XX:ParallelGCThreads=14 -Xmx32768m" \
            --INPUT ${HumanBam} \
            --OUTPUT ${QcDir}/${SeqID}.human.gatk.rnaseq.metrics \
            --REF_FLAT /storage/references_and_index/hg38/gtf/hg38.refFlat.txt \
            --STRAND_SPECIFICITY ${GatkStrandOption} \
            --CHART_OUTPUT ${QcDir}/${SeqID}.hg38.gatk.rnaseq.metric.pdf \
            --RIBOSOMAL_INTERVALS /storage/references_and_index/hg38/gtf/hg38.rRNA.interval_list \
            --REFERENCE_SEQUENCE /storage/references_and_index/hg38/fasta/rnaseq/hg38.fa

        ${RunGatk4} CollectRnaSeqMetrics \
            --java-options "-XX:ParallelGCThreads=14 -Xmx32768m" \
            --INPUT ${MouseBam} \
            --OUTPUT ${QcDir}/${SeqID}.mouse.gatk.rnaseq.metrics \
            --REF_FLAT /storage/references_and_index/mm10/gtf/mm10.refFlat.txt \
            --STRAND_SPECIFICITY ${GatkStrandOption} \
            --CHART_OUTPUT ${QcDir}/${SeqID}.mm10.gatk.rnaseq.metric.pdf \
            --RIBOSOMAL_INTERVALS /storage/references_and_index/mm10/gtf/mm10.rRNA.interval_list \
            --REFERENCE_SEQUENCE /storage/references_and_index/mm10/fasta/rnaseq/mm10.fa

    else

        # general-pipeline

        InputBam=${BamDir}/${SeqID}.Aligned.out.Sorted.bam
    
        ${RunGatk4} CollectRnaSeqMetrics \
            --java-options "-XX:ParallelGCThreads=14 -Xmx32768m" \
            --INPUT ${InputBam} \
            --OUTPUT ${QcDir}/${SeqID}.gatk.rnaseq.metrics \
            --REF_FLAT ${RefFlatFile} \
            --STRAND_SPECIFICITY ${GatkStrandOption} \
            --CHART_OUTPUT ${QcDir}/${SeqID}.gatk.rnaseq.metric.pdf \
            --RIBOSOMAL_INTERVALS ${RrnaIntervalList} \
            --REFERENCE_SEQUENCE ${GenomeFasta}
    fi
    #---------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------#

#---| RNA-SEQC |-----------------------------------------------------------------------------------#
    
    if [ "${XenograftPipeline}" = "true" ]; then

        # Xenograft piepline

        HumanBam=${BamDir}/${SeqID}.human.xenofilter.bam
        MouseBam=${BamDir}/${SeqID}.mouse.xenofilter.bam
        HumanGeneModel=/storage/references_and_index/hg38/gtf/hg38.stranded.gene.model.gtf
        MouseGeneModel=/storage/references_and_index/mm10/gtf/mm10.stranded.gene.model.gtf

        SeqcQuantDir=${DataDir}/${SeqID}/quant/rnaseqc
        mkdir -p ${SeqcQuantDir}
        
        ${RunRnaseqc} ${HumanGeneModel} ${HumanBam} ${SeqcQuantDir} \
            --sample ${SeqID}.human --stranded RF --coverage --detection-threshold=3 \
            --fasta /storage/references_and_index/hg38/fasta/rnaseq/hg38.fa
        
        ${RunRnaseqc} ${MouseGeneModel} ${MouseBam} ${SeqcQuantDir} \
            --sample ${SeqID}.mouse --stranded RF --coverage --detection-threshold=3 \
            --fasta /storage/references_and_index/mm10/fasta/rnaseq/mm10.fa

        cp ${SeqcQuantDir}/${SeqID}.human.metrics.tsv ${QcDir}/${SeqID}.human.rnaseqc.metrics.tsv
        cp ${SeqcQuantDir}/${SeqID}.human.gc_content.tsv ${QcDir}/${SeqID}.human.rnaseqc.gc_content.tsv
        cp ${SeqcQuantDir}/${SeqID}.human.coverage.tsv ${QcDir}/${SeqID}.human.rnaseqc.coverage.tsv

        cp ${SeqcQuantDir}/${SeqID}.mouse.metrics.tsv ${QcDir}/${SeqID}.mouse.rnaseqc.metrics.tsv
        cp ${SeqcQuantDir}/${SeqID}.mouse.gc_content.tsv ${QcDir}/${SeqID}.mouse.rnaseqc.gc_content.tsv
        cp ${SeqcQuantDir}/${SeqID}.mouse.coverage.tsv ${QcDir}/${SeqID}.mouse.rnaseqc.coverage.tsv

    else

        # general-pipeline

        SeqcQuantDir=${DataDir}/${SeqID}/quant/rnaseqc
        mkdir -p ${SeqcQuantDir}

        if [ "${Paired}" = "true" ]; then 
            ${RunRnaseqc} ${GeneModelGtfFile} ${BamDir}/${SeqID}.Aligned.out.Sorted.bam ${SeqcQuantDir} \
                --sample ${SeqID} --fasta ${GenomeFasta} --stranded RF --coverage --detection-threshold=3   
        else
            ${RunRnaseqc} ${GeneModelGtfFile} ${BamDir}/${SeqID}.Aligned.out.Sorted.bam ${SeqcQuantDir} \
                --sample ${SeqID} --fasta ${GenomeFasta} --stranded RF --coverage --detection-threshold=3 --unpaired  
        fi

        cp ${SeqcQuantDir}/${SeqID}.metrics.tsv ${QcDir}/${SeqID}.rnaseqc.metrics.tsv
        cp ${SeqcQuantDir}/${SeqID}.gc_content.tsv ${QcDir}/${SeqID}.rnaseqc.gc_content.tsv
        cp ${SeqcQuantDir}/${SeqID}.coverage.tsv ${QcDir}/${SeqID}.rnaseqc.coverage.tsv

    fi

#--------------------------------------------------------------------------------------------------#

#---| RUN MULTIQC |--------------------------------------------------------------------------------#
    if [ "${RunMultiqc}" = "true" ]; then       
        MultiqcConfig=/storage/home/kangsm/runScripts/WTS/WTS_P03-2_QC-POST_PROCESSING-MULTIQC-CONFIG.yaml
        #-----------------------------------------------------------------------
        ${RunSingularity}/multiqc-1.18.sif multiqc --force \
            --filename ${SeqID}_qc \
            --outdir ${QcDir}/qc_res \
            --config ${MultiqcConfig} \
            --data-dir ${QcDir}
    fi
#--------------------------------------------------------------------------------------------------#    

#---| QC-POST-PROCESSING |-------------------------------------------------------------------------#
    Rscript ${RunScript}/WTS_P03-3_QC-POST-PROCESSING.R \
        --baseDir ${BaseDir} \
        --batchID ${BatchID} \
        --seqID ${SeqID} \
        --databaseImport ${QcResDbImport}
#--------------------------------------------------------------------------------------------------#    


