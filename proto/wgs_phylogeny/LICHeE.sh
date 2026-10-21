

#! /bin/bash

#---| Input Options |------------------------------------------------------------------------------#
    usage() {
        echo "
        Usage: hmftools-purple-somatic-vcf-to-var-reads-table.sh [ ARGUMENTS ]...

        --BASE_DIR                  [ REQUIRED ] analysis base folder by NGS data type
        
        --BATCH_ID                  [ REQUIRED ] dataset batch-id
        
        --ANALYSIS_LIST             [ REQURIED ] analysis samples list

        --OUTPUT_DIR                [ REQURIED ] output folder
        
        -h | --help                 Print this usage 
        "
    }
#--------------------------------------------------------------------------------------------------#

#---| OPTION PARAMETERS PARSING --------------------------------------------------------------------
    ARGS=$(getopt -a -o h: --long BASE_DIR:,BATCH_ID:,ANALYSIS_LIST:,OUTPUT_DIR:,help -- "$@" )
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
        --ANALYSIS_LIST )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  AnalysisList=$2 ; shift 2 ;;
            esac ;;
        --OUTPUT_DIR )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  OutputDir=$2 ; shift 2 ;;
            esac ;;
        -h | --help )
        usage >&2 ; exit 2 ;;   
        --) shift ; break ;;
        *)  usage >&2 ; exit 2 ;;
    esac
    done
#---------------------------------------------------------------------------------------------------

### check params
    if [[ -z "${BaseDir}" ]]; then exit 2; fi
    if [[ -z "${BatchID}" ]]; then exit 2; fi
    if [[ -z "${AnalysisList}" ]]; then exit 2; fi
    if [[ -z "${OutputDir}" ]]; then exit 2; fi
###
### manual params
    BaseDir="/data/wgs"
    BatchID="WGS-26-01"
    AnalysisList="/data/wgs/WGS-26-01/meta/analysis.list"
    OutputDir="/data/wgs/WGS-26-01/analysis/pairtree"
###
### 
    lichee_dir=${BaseDir}/${BatchID}/analysis/lichee/wgs

    GenomeFasta="/storage/references_and_index/hg19/hmftools/Homo_sapiens.GRCh37.GATK.illumina.fasta"

    
    TumorList=($(cat ${AnalysisList} | grep -v "normal" | awk '{print $1}'))
    for TID in ${TumorList[@]}
    do
        hmftools_purple_vcf=${BaseDir}/${BatchID}/${TID}/hmftools/purple/${TID}.purple.somatic.vcf.gz

        bcftools norm -m -any -f "${GenomeFasta}" "${hmftools_purple_vcf}" \
        | bcftools view -f PASS -s ${TID} \
            -Oz -o ${lichee_dir}/${TID}.purple.somatic.filtered.vcf.gz
        tabix -p vcf ${lichee_dir}/${TID}.purple.somatic.filtered.vcf.gz
    done

    VCFS=()
    for TID in ${TumorList[@]}
    do
        VCFS+=("${lichee_dir}/${TID}.purple.somatic.filtered.vcf.gz")
    done
    
    bcftools merge -m none -Oz -o "${lichee_dir}/${BatchID}.tumors.variants.vcf.gz" "${VCFS[@]}"
    {
        printf "CHROM\tPOS\tREF\tALT\tFILTER"
        for ID in ${TumorList[@]}; do
            printf "\t%s" "$ID"
        done
        printf "\n"

        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%AD]\n' \
            "${lichee_dir}/${BatchID}.tumors.variants.vcf.gz"
    } > "${lichee_dir}/${BatchID}.tumors.variants.vaf.tsv" 




    Rscript /storage/home/kangsm/runScripts/00_Pipeline_Development/wgs_phylogeny/create-lichee-input.R \
        --WORK_DIR "${lichee_dir}" \
        --ANALYSIS_LIST "analysis.list" \
        --VAF_TABLE "${BatchID}.tumors.variants.vaf.tsv" 

    RunLichee="/storage/apps/lichee/LICHeE/release/lichee"
    ln -s /storage/apps/lichee/LICHeE/release/lichee.jar lichee.jar 

    ${RunLichee} -build \
        -i ${lichee_dir}/lichee_analysis_input.txt \
        -o ${lichee_dir} \
        -maxVAFAbsent 0.0005 \
        -minVAFPresent 0.0005 \
        -n 0 \
        -showTree 1


    cd lichee/LICHe

    mkdir -p build/temp
    cd build/temp
    find ../../lib -name "*.jar" -exec jar xf {} \;
    cd ../..

    javac -cp "lib/*" -d build/temp src/lineage/*.java src/util/*.java

    echo "Main-Class: lineage.LineageEngine" > build/MANIFEST.MF

    jar cfm release/lichee.jar build/MANIFEST.MF -C build/temp .
