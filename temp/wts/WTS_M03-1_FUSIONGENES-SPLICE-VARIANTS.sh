
#! /bin/bash

# WTS_run.FastqQC.sh <-- USE WES_run.FastqQC.sh

#---| SCRIPT USAGE |-------------------------------------------------------------------------------#
    usage() {
        echo "Usage: WTS_run.FusionGenes_and_SpliceVariants.sh [ARGUMENTS]...
        [ --BASE_DIR          ] : base folder. (data type or work type)
        [ --SEQ_FOLDER        ] : dataset id folder ( = analysis batch folder )
        [ --SEQ_ID            ] : sample ID ( = seq id )
        [ --THREADS           ] : N-CPU threads
        [ --GENOME_ASSEMBLY   ] : genome assembly. default = 'hg19'
        [ --RUN_ARRIBA        ] : run arriba. default = 'true'
        [ --RUN_STARFUSION    ] : run star-fusion. default = 'true'
        [ --RUN_SPLICEVARIANT ] : run ctat_splicing. default = 'true'
        [ --DRAW_ARRIBA_PLOTS ] : draw fusion-gene plots using arriba. default = 'true'
        [ --PAIRED_END        ] : sequencing type. true = PE, false = SE. default = 'true'
        [ -h | --help       ] : Print this usage 
        "
    }
#--------------------------------------------------------------------------------------------------#

#---| DEFAULT VALUES |-----------------------------------------------------------------------------#
    BASE_DIR="/data/wts"
    THREADS=15
    GENOME_ASSEMBLY="hg19"
    RUN_ARRIBA="true"
    RUN_STARFUSION="true"
    RUN_SPLICEVARIANT="true"
    DRAW_ARRIBA_PLOTS="true"
    PAIRED_END="true"
#--------------------------------------------------------------------------------------------------#

#---| OPTION PARSE |-------------------------------------------------------------------------------#
    ARGS=$(getopt -a -o h: --long SEQ_FOLDER:,SEQ_ID:,THREADS:,BASE_DIR:,RUN_ARRIBA:,RUN_STARFUSION:,RUN_SPLICEVARIANT:,GENOME_ASSEMBLY:,DRAW_ARRIBA_PLOTS:,PAIRED_END:,help -- "$@" )
    VALID_ARGS=$?
    if [ "$VALID_ARGS" != "0" ]; then 
        usage >&2 
        exit 2
    fi

    eval set -- "$ARGS"
    while :
    do
    case "$1" in
        --SEQ_FOLDER )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  SEQ_FOLDER=$2 ; shift 2 ;;
            esac ;;
        --SEQ_ID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  SEQ_ID=$2 ; shift 2 ;;
            esac  ;;
        --THREADS )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 15" ; exit 2 ;;
                *)  THREADS=$2 ; shift 2 ;;
            esac ;;
        --BASE_DIR )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = /data/wts" ; exit 2 ;;
                *)  BASE_DIR=$2 ; shift 2 ;;
            esac ;;
        --GENOME_ASSEMBLY )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = hg19" ; exit 2 ;;
                *)  GENOME_ASSEMBLY=$2 ; shift 2 ;;
            esac ;; 
        --RUN_ARRIBA )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = true" ; exit 2 ;;
                *)  RUN_ARRIBA=$2 ; shift 2 ;;
            esac ;; 
        --RUN_STARFUSION )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = true" ; exit 2 ;;
                *)  RUN_STARFUSION=$2 ; shift 2 ;;
            esac ;; 
        --RUN_SPLICEVARIANT )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = true" ; exit 2 ;;
                *)  RUN_SPLICEVARIANT=$2 ; shift 2 ;;
            esac ;; 
        --DRAW_ARRIBA_PLOTS )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = true" ; exit 2 ;;
                *)  DRAW_ARRIBA_PLOTS=$2 ; shift 2 ;;
            esac ;; 
        --PAIRED_END )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = true" ; exit 2 ;;
                *)  PAIRED_END=$2 ; shift 2 ;;
            esac ;;
        -h | --help )
        usage >&2 ; exit 2 ;;   
        --) shift ; break ;;
        *)  usage >&2 ; exit 2 ;;
    esac
    done
#--------------------------------------------------------------------------------------------------#

#---| FOLDERS |------------------------------------------------------------------------------------#
    DATA_DIR=${BASE_DIR}/${SEQ_FOLDER}
    ##  
    FASTQ_DIR=${DATA_DIR}/${SEQ_ID}/fastq
    BAM_DIR=${DATA_DIR}/${SEQ_ID}/bam
    FUSIONGENE_RES=${DATA_DIR}/${SEQ_ID}/fusionGenes_spliceVariants
    ARRIBA_RES_DIR=${FUSIONGENE_RES}/arriba
    STARFUSION_RES_DIR=${FUSIONGENE_RES}/star_fusion
    SPLICE_VAR_RES_DIR=${FUSIONGENE_RES}/splice_variants
    LOG_DIR=${DATA_DIR}/${SEQ_ID}/log
    ##
    if [ ! -d ${FUSIONGENE_RES} ]; then mkdir -p ${FUSIONGENE_RES}; fi
   
    if [ ! -d ${LOG_DIR} ]; then mkdir -p ${LOG_DIR}; fi
#--------------------------------------------------------------------------------------------------#

#---| REFERENCES AND INDEX |-----------------------------------------------------------------------#
    GENOME_FASTA=/storage/references_and_index/${GENOME_ASSEMBLY}/fasta/human_g1k_v37_decoy.fasta
    GTF_FILE=/storage/references_and_index/${GENOME_ASSEMBLY}/gtf/${GENOME_ASSEMBLY}.gtf
    ## Arriba
    ARRIBA_REFERENCE=/storage/references_and_index/${GENOME_ASSEMBLY}/arriba
    ARRIBA_BLACKLIST=${ARRIBA_REFERENCE}/blacklist_${GENOME_ASSEMBLY}.tsv.gz
    ARRIBA_KNOWN_FUSIONS=${ARRIBA_REFERENCE}/known_fusions_${GENOME_ASSEMBLY}.tsv.gz
    ARRIBA_PROTEIN_FEATURES=${ARRIBA_REFERENCE}/protein_domains_${GENOME_ASSEMBLY}.gff3
    ARRIBA_CYTOBAND=${ARRIBA_REFERENCE}/cytobands_${GENOME_ASSEMBLY}.tsv
    ARRIBA_GTF=${ARRIBA_REFERENCE}/${GENOME_ASSEMBLY}.gtf
    ## STAR-Fusion & CTAT_SPLICING
    CTAT_RESOURCE_LIB=/storage/references_and_index/${GENOME_ASSEMBLY}/star_fusion/CTAT_resource_lib
#--------------------------------------------------------------------------------------------------#

#---| SWTOOLS RUN |--------------------------------------------------------------------------------#
    RUN_SINGULARITY="singularity exec -B /storage,/data /storage/images"
    ARRIBA_RUN="${RUN_SINGULARITY}/arriba-2.4.0.sif /arriba_v2.4.0"
    STAR_FUSION_RUN="${RUN_SINGULARITY}/star-fusion-0.13.0.sif STAR-Fusion"
    CTAT_SPLICING_RUN="${RUN_SINGULARITY}/ctat-splicing-0.0.2.sif /usr/local/src/CTAT-SPLICING"
    RUN_SAMTOOLS="${RUN_SINGULARITY}/samtools-1.18.sif samtools"
#--------------------------------------------------------------------------------------------------#

#---| LOGFILE |------------------------------------------------------------------------------------#
    LOGFILE=${LOG_DIR}/${SEQ_ID}.fusion.genes.and.splicing.variants.$(date '+%Y%m%d').log
    if [ ! -f ${LOGFILE} ]; then touch ${LOGFILE}; fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo " " >> ${LOGFILE}
    echo "WTS Fusion Genes and Splicing Variants log" >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " RUN DATE             : $(date '+%Y-%m-%d %H:%M:%S')" >> ${LOGFILE}
    echo " SEQ FOLDER           : ${DATA_DIR}" >> ${LOGFILE}
    echo " SEQ ID               : ${SEQ_ID}" >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " " >> ${LOGFILE}   
#==================================================================================================# 

#---| RUN ARRIBA |---------------------------------------------------------------------------------#
    if [ "${RUN_ARRIBA}" = "true" ]; then
        #==========================================================================================# 
            echo "$(date '+%Y-%m-%d %H:%M:%S') | FusionGenes: Arriba START." >> ${LOGFILE}
        #==========================================================================================#
        if [ ! -d ${ARRIBA_RES_DIR} ]; then mkdir -p ${ARRIBA_RES_DIR}; fi
        #----------------------------------------------------------------------#
        ${ARRIBA_RUN}/arriba \
            -x ${BAM_DIR}/${SEQ_ID}.Aligned.out.Sorted.bam \
            -o ${ARRIBA_RES_DIR}/${SEQ_ID}.arriba.fusion.genes.tsv \
            -O ${ARRIBA_RES_DIR}/${SEQ_ID}.arriba.filtered.out.fusion.genes.tsv \
            -a ${GENOME_FASTA} \
            -g ${ARRIBA_GTF} \
            -b ${ARRIBA_BLACKLIST} \
            -k ${ARRIBA_KNOWN_FUSIONS} \
            -p ${ARRIBA_PROTEIN_FEATURES} \
            -i "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y" 
        #==========================================================================================# 
            echo "$(date '+%Y-%m-%d %H:%M:%S') | FusionGenes: Arriba FINISHED." >> ${LOGFILE}
        #==========================================================================================#
    else
        #==========================================================================================# 
            echo "$(date '+%Y-%m-%d %H:%M:%S') | FusionGenes: Arriba SKIPPED." >> ${LOGFILE}
        #==========================================================================================#
    fi 
#--------------------------------------------------------------------------------------------------#

#---| RUN STAR-FUSION |----------------------------------------------------------------------------#
    if [ "${RUN_STARFUSION}" = "true" ]; then
        #==========================================================================================# 
            echo "$(date '+%Y-%m-%d %H:%M:%S') | FusionGenes: STAR-Fusion START." >> ${LOGFILE}
        #==========================================================================================# 
        if [ ! -d ${STARFUSION_RES_DIR} ]; then mkdir -p ${STARFUSION_RES_DIR}; fi
        #----------------------------------------------------------------------#
        if [ "${PAIRED_END}" = "true" ]; then
            ${STAR_FUSION_RUN} \
                --genome_lib_dir ${CTAT_RESOURCE_LIB} \
                --left_fq ${FASTQ_DIR}/${SEQ_ID}.trimmed_R1.fastq.gz \
                --right_fq ${FASTQ_DIR}/${SEQ_ID}.trimmed_R2.fastq.gz \
                --output_dir ${STARFUSION_RES_DIR} \
                --CPU ${THREADS} \
                --STAR_twopass    
        else
            ${STAR_FUSION_RUN} \
                --genome_lib_dir ${CTAT_RESOURCE_LIB} \
                --left_fq ${FASTQ_DIR}/${SEQ_ID}.trimmed.fastq.gz \
                --output_dir ${STARFUSION_RES_DIR} \
                --CPU ${THREADS} \
                --STAR_twopass 
        fi
        #==========================================================================================# 
            echo "$(date '+%Y-%m-%d %H:%M:%S') | FusionGenes: STAR-Fusion FINISHED." >> ${LOGFILE}
        #==========================================================================================#
    else
        #==========================================================================================# 
            echo "$(date '+%Y-%m-%d %H:%M:%S') | FusionGenes: STAR-Fusion SKIPPED." >> ${LOGFILE}
        #==========================================================================================#
    fi 
#--------------------------------------------------------------------------------------------------#

#---| RUN CTAT_SPLICING |--------------------------------------------------------------------------#
    if [ "${RUN_SPLICEVARIANT}" = "true" ]; then
        #==========================================================================================# 
            echo "$(date '+%Y-%m-%d %H:%M:%S') | SpliceVariants: CTAT_SPLICING START." >> ${LOGFILE}
        #==========================================================================================#
        if [ ! -d ${SPLICE_VAR_RES_DIR} ]; then mkdir -p ${SPLICE_VAR_RES_DIR}; fi
        #----------------------------------------------------------------------#
        if [ -f ${STARFUSION_RES_DIR}/Chimeric.out.junction ]; then
            #======================================================================================# 
                echo "$(date '+%Y-%m-%d %H:%M:%S') | STAR-Fusion junction information files found." >> ${LOGFILE}
                echo "$(date '+%Y-%m-%d %H:%M:%S') | Used STAR-Fusion junction information files." >> ${LOGFILE}
                echo "STAR-Fusion junction information files found."
            #======================================================================================#
        else
            #======================================================================================# 
                echo "$(date '+%Y-%m-%d %H:%M:%S') | No STAR-Fusion chimeric junction file found." >> ${LOGFILE}
                echo "$(date '+%Y-%m-%d %H:%M:%S') | Run STAR-Fusion for junction information files." >> ${LOGFILE}
                echo "STAR-Fusion junction information files NOT found."
                echo "Run STAR-Fusion for junction information files."
            #======================================================================================#
            if [ ! -d ${STARFUSION_RES_DIR} ]; then mkdir -p ${STARFUSION_RES_DIR}; fi
            #----------------------------------------------------------------------#
            if [ "${PAIRED_END}" = "true" ]; then
                ${STAR_FUSION_RUN} \
                    --genome_lib_dir ${CTAT_RESOURCE_LIB} \
                    --left_fq ${FASTQ_DIR}/${SEQ_ID}.trimmed_R1.fastq.gz \
                    --right_fq ${FASTQ_DIR}/${SEQ_ID}.trimmed_R2.fastq.gz \
                    --output_dir ${STARFUSION_RES_DIR} \
                    --CPU ${THREADS} \
                    --STAR_twopass    
            else
                ${STAR_FUSION_RUN} \
                    --genome_lib_dir ${CTAT_RESOURCE_LIB} \
                    --left_fq ${FASTQ_DIR}/${SEQ_ID}.trimmed.fastq.gz \
                    --output_dir ${STARFUSION_RES_DIR} \
                    --CPU ${THREADS} \
                    --STAR_twopass 
            fi 
            #======================================================================================# 
                echo "$(date '+%Y-%m-%d %H:%M:%S') | Run STAR-Fusion for junction information files Done." >> ${LOGFILE}
            #======================================================================================#
        fi  

        ${RUN_SAMTOOLS} sort --threads ${THREADS} --output-fmt BAM -o ${STARFUSION_RES_DIR}/Aligned.out.Sorted.bam ${STARFUSION_RES_DIR}/Aligned.out.bam  
        ${RUN_SAMTOOLS} index --threads ${THREADS} -b ${STARFUSION_RES_DIR}/Aligned.out.Sorted.bam
        
        ${CTAT_SPLICING_RUN}/STAR_to_cancer_introns.py \
            --ctat_genome_lib ${CTAT_RESOURCE_LIB} \
            --SJ_tab_file  ${STARFUSION_RES_DIR}/SJ.out.tab \
            --chimJ_file ${STARFUSION_RES_DIR}/Chimeric.out.junction \
            --vis \
            --bam_file ${STARFUSION_RES_DIR}/Aligned.out.Sorted.bam \
            --output_prefix ${SPLICE_VAR_RES_DIR}/${SEQ_ID} \
            --sample_name ${SEQ_ID}    
        #==========================================================================================# 
            echo "$(date '+%Y-%m-%d %H:%M:%S') | SpliceVariants: CTAT_SPLICING FINISHED." >> ${LOGFILE}
        #==========================================================================================#

    else
        #==========================================================================================# 
            echo "$(date '+%Y-%m-%d %H:%M:%S') | SpliceVariants: CTAT_SPLICING SKIPPED." >> ${LOGFILE}
        #==========================================================================================#
    fi 
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | Result Integration and Summarization START." >> ${LOGFILE}
#==================================================================================================# 

#---| RESULT INTEGRATION & SUMMARY + CIRCOS-PLOT |-------------------------------------------------#
    Rscript /storage/home/kangsm/runScripts/WTS_run.FusionGenes_SpliceVariants_Summary.R \
        --BASE_DIR ${BASE_DIR} \
        --SEQ_FOLDER ${SEQ_FOLDER} \
        --SEQ_ID ${SEQ_ID} \
        --DRAW_STAR_FUSION_CIRCOS TRUE \
        --GENOME_ASSEMBLY ${GENOME_ASSEMBLY}
#--------------------------------------------------------------------------------------------------#

#---| ARRIBA-FORMAT FUSION-GENES PLOTS |-----------------------------------------------------------#
    if [ "${DRAW_ARRIBA_PLOTS}" = "true" ]; then
        TARGET_FUSION_ARRIBA_FILE=${FUSIONGENE_RES}/${SEQ_ID}.FusionGenes.Common.Arriba.Format.tsv
        if [ ! -f ${TARGET_FUSION_ARRIBA_FILE} ]; then
            TARGET_FUSION_ARRIBA_FILE=${FUSIONGENE_RES}/arriba/${SEQ_ID}.arriba.fusion.genes.tsv
            if [ ! -f ${TARGET_FUSION_ARRIBA_FILE} ]; then 
                echo "No Arriba-Format FusionGene Result File. Plotting Skipped."
            else
                FUSION_GENE_PLOT_RES_DIR=${FUSIONGENE_RES}/FusionGenesPlots
                if [ ! -d ${FUSION_GENE_PLOT_RES_DIR} ]; then mkdir -p ${FUSION_GENE_PLOT_RES_DIR}; fi
            
                ${ARRIBA_RUN}/draw_fusions.R \
                    --fusions=${TARGET_FUSION_ARRIBA_FILE} \
                    --alignments=${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/bam/${SEQ_ID}.Aligned.out.Sorted.bam \
                    --output=${FUSION_GENE_PLOT_RES_DIR}/${SEQ_ID}.FusionGene.ALL.Plots.pdf \
                    --annotation=${ARRIBA_GTF} \
                    --cytobands=${ARRIBA_CYTOBAND} \
                    --proteinDomains=${ARRIBA_PROTEIN_FEATURES} \
                    --fontSize=1 2>&1 | tee ${FUSION_GENE_PLOT_RES_DIR}/plot.log
                #------------------------------------------------------------------#
                grep "^Drawing fusion" ${FUSION_GENE_PLOT_RES_DIR}/plot.log > ${FUSION_GENE_PLOT_RES_DIR}/plot.arriba.only.fusionGenes.list
                #------------------------------------------------------------------#
                for PAGENUM in $(seq 1 $(cat ${FUSION_GENE_PLOT_RES_DIR}/plot.arriba.only.fusionGenes.list | wc -l))
                do
                    FUSION_NAME=$(sed -n -e "${PAGENUM} p" ${FUSION_GENE_PLOT_RES_DIR}/plot.arriba.only.fusionGenes.list | awk '{print $4}' | sed 's/:/-/' -)
                    pdftoppm -f ${PAGENUM} -l ${PAGENUM} -png ${FUSION_GENE_PLOT_RES_DIR}/${SEQ_ID}.FusionGene.ALL.Plots.pdf > ${FUSION_GENE_PLOT_RES_DIR}/${SEQ_ID}.arriba.only.FusionGene.${FUSION_NAME}.plot.png
                    echo $FUSION_NAME
                done
            fi
        else
            FUSION_GENE_PLOT_RES_DIR=${FUSIONGENE_RES}/FusionGenesPlots
            if [ ! -d ${FUSION_GENE_PLOT_RES_DIR} ]; then mkdir -p ${FUSION_GENE_PLOT_RES_DIR}; fi
           
            ${ARRIBA_RUN}/draw_fusions.R \
                --fusions=${TARGET_FUSION_ARRIBA_FILE} \
                --alignments=${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/bam/${SEQ_ID}.Aligned.out.Sorted.bam \
                --output=${FUSION_GENE_PLOT_RES_DIR}/${SEQ_ID}.FusionGene.ALL.Plots.pdf \
                --annotation=${ARRIBA_GTF} \
                --cytobands=${ARRIBA_CYTOBAND} \
                --proteinDomains=${ARRIBA_PROTEIN_FEATURES} \
                --fontSize=1 2>&1 | tee ${FUSION_GENE_PLOT_RES_DIR}/plot.log
            #------------------------------------------------------------------#
            grep "^Drawing fusion" ${FUSION_GENE_PLOT_RES_DIR}/plot.log > ${FUSION_GENE_PLOT_RES_DIR}/plot.fusionGenes.list
            #------------------------------------------------------------------#
            for PAGENUM in $(seq 1 $(cat ${FUSION_GENE_PLOT_RES_DIR}/plot.fusionGenes.list | wc -l))
            do
                FUSION_NAME=$(sed -n -e "${PAGENUM} p" ${FUSION_GENE_PLOT_RES_DIR}/plot.fusionGenes.list | awk '{print $4}' | sed 's/:/-/' -)
                pdftoppm -f ${PAGENUM} -l ${PAGENUM} -png ${FUSION_GENE_PLOT_RES_DIR}/${SEQ_ID}.FusionGene.ALL.Plots.pdf > ${FUSION_GENE_PLOT_RES_DIR}/${SEQ_ID}.FusionGene.${FUSION_NAME}.plot.png
                echo $FUSION_NAME
            done
            #------------------------------------------------------------------#
        fi
    else
        #==========================================================================================#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | Creation of Arriba-Format FusionGene Plot Skipped." >> ${LOGFILE}
        #==========================================================================================#
    fi
# #--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo " " >> ${LOGFILE}
    echo "WTS Fusion Genes and Splicing Variants DONE." >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " " >> ${LOGFILE}
#==================================================================================================# 




