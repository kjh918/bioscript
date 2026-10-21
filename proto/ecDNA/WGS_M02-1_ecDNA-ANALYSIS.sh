
#! /bin/bash

#---| USAGE |---------------------------------------------------------------------------------------
    usage() {echo "
        Usage: WGS_M02-1_ecDNA-ANALYSIS.sh [ ARGUMENTS ]...
        
        --BASE_DIR                  [ REQUIRED ] analysis data type base folder

        --BATCH_ID                  [ REQUIRED ] analysis batch-id 
        
        --SEQ_ID                    [ REQUIRED ] sample-id 
        
        --THREADS                   parallel processing cpu cores
                                    default = 15
        
        --GENOME_ASSEMBLY           reference genome version. hg19 or hg38
                                    default = hg38
        
        --RUN_AMPLICON_ARCHITECT    run 'amplicon architect'
                                    default = true
        
        --RUN_CIRCLE_MAP            run 'circle-map'.
                                    default = false

        --RUN_CIRCLE_FINDER         run 'circle-finder'. 
                                    default = false

        -h | --help        : Print this usage 
        "
    }
#---------------------------------------------------------------------------------------------------

#---| OPTION PARAMETERS PARSING --------------------------------------------------------------------
    ARGS=$(getopt -a -o d:s:h: --long BASE_DIR:,BATCH_ID:,SEQ_ID:,THREADS:,GENOME_ASSEMBLY:,RUN_AMPLICON_ARCHITECT:,RUN_CIRCLE_MAP:,RUN_CIRCLE_FINDER:,help -- "$@" )
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
        --SEQ_ID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  SeqID=$2 ; shift 2 ;;
            esac ;;
        --THREADS )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  Threads=$2 ; shift 2 ;;
            esac ;;
        --GENOME_ASSEMBLY )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = hg38" ; shift ;;
                *)  GenomeAssembly=$2 ; shift 2 ;;
            esac ;;
        --RUN_AMPLICON_ARCHITECT )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = true" ; shift ;;
                *)  RunAmpliconArchitect=$2 ; shift 2 ;;
            esac ;;
        --RUN_CIRCLE_MAP )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = false" ; shift ;;
                *)  RunCircleMap=$2 ; shift 2 ;;
            esac ;;
        --RUN_CIRCLE_FINDER )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = false" ; shift ;;
                *)  RunCircleFinder=$2 ; shift 2 ;;
            esac ;;
        -h | --help )
        usage >&2 ; exit 2 ;;   
        --) shift ; break ;;
        *)  usage >&2 ; exit 2 ;;
    esac
    done
#--------------------------------------------------------------------------------------------------

#---| Manual Params --------------------------------------------------------------------------------
    # BaseDir="/data/wgs"
    # BatchID="seqc2_wgs"
    # SeqID="HCC1395"
    # Threads=16
    # GenomeAssembly="hg38"
    # RunAmpliconArchitect="true"
    # RunCircleMap="false"
    # RunCircleFinder="false"
#---------------------------------------------------------------------------------------------------

#---| Param Check ----------------------------------------------------------------------------------
    if [ -z "${BaseDir}" ]; then echo "No 'BASE_DIR'. REQUIRED." ; exit 2; fi
    if [ -z "${BatchID}" ]; then echo "No 'BATCH_ID'. REQUIRED." ; exit 2; fi
    if [ -z "${SeqID}" ]; then echo "No 'SEQ_ID'. REQUIRED." ; exit 2; fi
    if [ -z "${Threads}" ]; then Threads=15 fi
    if [ -z "${GenomeAssembly}" ]; GenomeAssembly="hg38" fi
    if [ -z "${RunAmpliconArchitect}" ]; then RunAmpliconArchitect="ture" fi
    if [ -z "${RunCircleMap}" ]; then RunCircleMap="false" fi
    if [ -z "${RunCircleFinder}" ]; then RunCircleFinder="false" fi
#---------------------------------------------------------------------------------------------------

#---| Swtools Run ----------------------------------------------------------------------------------
    
    # environment variables
    export MOSEKLM_LICENSE_FILE=/storage/home/kangsm/myDB/Mosek_Lisence/Trial_20260309/mosek.lic
    export AA_DATA_REPO=/storage/references_and_index/${GenomeAssembly}/amplicon_architect
    
    # app-image and run command
    AmpliconArchitect_SIF="/storage/images/ampliconsuite-pipeline_v1.5.0.sif"
    RunCircleMap="singularity exec -B /storage,/data /storage/images/circleMap.sif"
    RunCircleFinder="/storage/apps/circle_finder/circle_finder-pipeline-bwa-mem-samblaster.sh"

    # resources
    GenomeFasta=/storage/references_and_index/${GenomeAssembly}/fasta/${GenomeAssembly}.primary_decoy_HLA_no.ALT.fa
    GenomeIndex=/storage/references_and_index/${GenomeAssembly}/fasta/${GenomeAssembly}.primary_decoy_HLA_no.ALT.fa

#---------------------------------------------------------------------------------------------------

#---| Input Files and Folders ----------------------------------------------------------------------
    
    # BAM file folder
    BamDir=${BaseDir}/${BatchID}/${SeqID}/bam
    
#---------------------------------------------------------------------------------------------------
    
#---| RUN AMPLICON-CRCHITECT -----------------------------------------------------------------------
    if [ "$RunAmpliconArchitect" = "true" ]; then

        # reference genome tag
        if [ "${GenomeAssembly}" = "hg19" ]; then 
            AA_ReferenceVersion="GRCh37"
        else
            AA_ReferenceVersion="GRCh38"
        fi

        # result folder
        AmpArchResDir=${BaseDir}/${BatchID}/${SeqID}/ecDNA/ampliconarchitect
        mkdir -p ${AmpArchResDir}
        cd ${AmpArchResDir}

        # input bam file
        AmpAarchInputBam=${BaseDir}/${BatchID}/${SeqID}/bam/${SeqID}.analysisReady.bam

        # run
        python3 /storage/apps/AmpliconSuite-pipeline/singularity/run_paa_singularity.py \
            --sif ${AA_SIF} \
            --output_directory ${AmpArchResDir} \
            --sample_name ${SeqID} \
            --nthreads ${Threads} \
            --ref ${AA_ReferenceVersion} \
            --bam ${AmpAarchInputBam} \
            --run_AA \
            --run_AC 
    
    fi

#---------------------------------------------------------------------------------------------------

#---| RUN CIRCLE-MAP -------------------------------------------------------------------------------

    if [ "$RunCircleMap" = "true" ]; then

        # reulst folder
        CircleMapDir=${BaseDir}/${BatchID}/${SeqID}/ecDNA/circlemap
        mkdir -p ${CircleMapDir}

        # create queryname sorted bam 
        samtools sort -@ ${Threads} -n -O BAM ${BamDir}/${SeqID}.analysisReady.bam \
            | samtools fixmate -@ ${Threads} -m - - \
            | samtools view -@ ${Threads} -b -f1 -F12 \
            | samtools sort -@ ${Threads} -n -o ${BamDir}/${SeqID}.analysisReady.name.sorted.bam

        # Run Circle-Map (1)
        ${RunCircleMap} Circle-Map ReadExtractor \
            -i ${SeqID}.analysisReady.name.sorted.bam \
            -o ${SeqID}.circlemap.candidate.reads.bam \
            -dir ${BamDir}

        # Run Circle-Map (2)
        samtools sort --threads ${Threads} -o ${SeqID}.circlemap.candidate.reads.sorted.bam ${SeqID}.circlemap.candidate.reads.bam
        samtools index --threads ${Threads} ${SeqID}.circlemap.candidate.reads.sorted.bam 
        
        # Run Circle-Map (3)
        cd ${BamDir}
        ${RunCircleMap} Circle-Map Realign \
            -i ${SeqID}.circlemap.candidate.reads.sorted.bam  \
            -qbam ${SeqID}.analysisReady.name.sorted.bam \
            -sbam ${SeqID}.analysisReady.bam \
            -fasta ${GenomeFasta} \
            -o ${SeqID}.ecDNA.circlemap.result.bed \
            -dir ${BamDir} \
            -t ${Threads}

        # copy initial result to circle-map folder
        cp ${BamDir}/${SeqID}.ecDNA.circlemap.result.bed ${CircleMapDir}/${SeqID}.ecDNA.circlemap.result.bed

    fi

#---------------------------------------------------------------------------------------------------

#---| RUN CIRCLE-FINDER ----------------------------------------------------------------------------

    if [ "$RunCircleFinder" = "true" ]; then

        # result folder
        CircleFinderDir=${BaseDir}/${BatchID}/${SeqID}/ecDNA/circlefinder
        mkdir -p ${CircleFinderDir}
        cd ${CircleFinderDir}

        RunCircleFinderInputFastqR1=${BaseDir}/${BatchID}/${SeqID}/fastq/${SeqID}.trimmed_R1.fastq.gz
        RunCircleFinderInputFastqR2=${BaseDir}/${BatchID}/${SeqID}/fastq/${SeqID}.trimmed_R2.fastq.gz

        # run 
        ${RunCircleFinder} ${Threads} ${GenomeIndex} ${RunCircleFinderInputFastqR1} ${RunCircleFinderInputFastqR2} 10 ${SeqID} ${GenomeAssembly}

    fi

#---------------------------------------------------------------------------------------------------

#---| Post-Analysis --------------------------------------------------------------------------------
    
    # Rscript /storage/home/kangsm/runScripts/00_Pipeline_Development/ecDNA/WGS_M02-2_ecDNA-POST-ANALYSIS.R \

#---------------------------------------------------------------------------------------------------




