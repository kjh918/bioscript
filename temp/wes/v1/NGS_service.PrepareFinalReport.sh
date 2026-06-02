

#! /bin/bash 

# Script Usage -------------------------------------------------------------------------------------
    usage() {
        echo "Usage: NGS_service.PrepareFinalReport.sh [ ARGUMENTS ]...
            --InputPdf = input initial report PDF file. (only filename. not include absolute path) 
            --LocalRun = run script in local PC. default = 'true'
            --ReportType = analysis type. 'standard' or 'advanced'.
            --NgsApplication = NGS application type. 'wes', 'wts' or 'tso'
            --OrderId = NGS service order id
            --ReportResourceDir = report resource folder
            --OutputDir = final report output folder
            --RemoveTmp = remove temp folder
            -h | --help = print this usage 
        "
    }
#---------------------------------------------------------------------------------------------------
# Default Parameters -------------------------------------------------------------------------------
    REPORT_RESOURCE_DIR=""
    LOCAL_RUN="true"
    OUTPUT_DIR="./"
    REMOVE_TMP="false"
#---------------------------------------------------------------------------------------------------
# Arguments ----------------------------------------------------------------------------------------
    ARGS=$(getopt -a -o h: --long InputPdf:,LocalRun:,ReportType:,NgsApplication:,OrderId:,ReportResourceDir:,OutputDir:,RemoveTmp:,help -- "$@" )
    VALID_ARGS=$?
    if [ "$VALID_ARGS" != "0" ]; then 
        usage >&2 
        exit 2
    fi

    eval set -- "$ARGS"
    while :
    do
    case "$1" in
        --InputPdf )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  INPUT_PDF=$2 ; shift 2 ;;
            esac ;;
        --LocalRun )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = true " ; shift ;;
                *)  LOCAL_RUN=$2 ; shift 2 ;;
            esac ;;
        --ReportType )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  REPORT_TYPE=$2 ; shift 2 ;;
            esac ;;
        --NgsApplication )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  NGS_APPLICATION=$2 ; shift 2 ;;
            esac ;;
        --OrderId )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  ORDER_ID=$2 ; shift 2 ;;
            esac ;;
        --ReportResourceDir )
            case "$2" in
                -* | --* | "") echo "no value at $1, default folder will be used" ; shift ;;
                *)  REPORT_RESOURCE_DIR=$2 ; shift 2 ;;
            esac ;;
        --OutputDir )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = './''" ; shift ;;
                *)  OUTPUT_DIR=$2 ; shift 2 ;;
            esac ;;
        --RemoveTmp )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = './''" ; shift ;;
                *)  REMOVE_TMP=$2 ; shift 2 ;;
            esac ;;
        -h | --help )
        usage >&2 ; exit 2 ;;   
        --) shift ; break ;;
        *)  usage >&2 ; exit 2 ;;
    esac
    done
#---------------------------------------------------------------------------------------------------

# Manual Params ------------------------------------------------------------------------------------
    INPUT_PDF="advanced_report.pdf"
    LOCAL_RUN="true"
    REPORT_TYPE="advanced"
    NGS_APPLICATION="wts"
    ORDER_ID="GCX-C11-251215"
    OUTPUT_DIR="./"
    #REPORT_RESOURCE_DIR=""
    REMOVE_TMP="false"
#---------------------------------------------------------------------------------------------------

# Resources Folder ---------------------------------------------------------------------------------
    if [[ ! -n ${REPORT_RESOURCE_DIR} ]]; then
        # default resource dir if not specified in args
        if [ "${LOCAL_RUN}" = "true" ]; then
            RESOURCE_DIR="/mnt/d/01.NGS_RUO_SERVICE/ReportResources"
        else
            RESOURCE_DIR="/storage/home/kangsm/shinyWeb/resources"
        fi
    else
        # use resource dir from args
        RESOURCE_DIR=${REPORT_RESOURCE_DIR}
    fi
#---------------------------------------------------------------------------------------------------

# Verify Resource Folder
    # if [ ! -f ${RESOURCE_DIR}/default_stamp.tex ]; then 
    #     echo "No required files found in Resource folder. Check again. Stop." >&2 
    #     exit 2
    # fi
#

# Cover PDF ----------------------------------------------------------------------------------------
    if [ "${NGS_APPLICATION}" = "wes" ]; then COVER_PDF=${RESOURCE_DIR}/wes_${REPORT_TYPE}_report_cover.pdf; fi
    if [ "${NGS_APPLICATION}" = "wts" ]; then COVER_PDF=${RESOURCE_DIR}/wts_report_cover.pdf; fi
    if [ "${NGS_APPLICATION}" = "tso" ]; then COVER_PDF=${RESOURCE_DIR}/tso_report_cover.pdf; fi
    if [ "${NGS_APPLICATION}" = "qc"  ]; then COVER_PDF=${RESOURCE_DIR}/ngs_qc_report_cover.pdf; fi
#---------------------------------------------------------------------------------------------------

# PROCESSING ---------------------------------------------------------------------------------------
    # Output Folder and Temp Folder
    cd ${OUTPUT_DIR}
    mkdir -p ./tmp
    # Copy Input PDF into Temp Folder
    cp ./${INPUT_PDF} tmp/
    # Change Folder
    cd ./tmp
    # Remove Front Page from Initial Report PDF
    pdftk ${INPUT_PDF} cat 2-end output rm_cover_report.pdf
    # Count Report Pages
    REPORT_PAGES=$(pdftk rm_cover_report.pdf dump_data | grep NumberOfPages | awk '{print $2}')
    # Prepapre Stamp PDF -------------------------------------------------------------------------------
    cp ${RESOURCE_DIR}/default_stamp.tex ./report_stamp.tex
    sed -i -e "s/ORDER_ID_HERE/${ORDER_ID}/" report_stamp.tex
    sed -i -e "s/PDF_PAGE_HERE/${REPORT_PAGES}/" report_stamp.tex
    sed -i 's/_/-/g' report_stamp.tex
    pdflatex -output-format "pdf_document" -jobname "report_stamp" report_stamp.tex
    # Add Stamp PDF to Report PDF
    pdftk rm_cover_report.pdf multistamp report_stamp.pdf output rm_cover_paged_report.pdf
    # Add Cover to Report PDF
    pdftk ${COVER_PDF} rm_cover_paged_report.pdf cat output rm_cover_paged_add_cover_report.pdf
    # Compression Final Report PDF
    if [ "${REPORT_TYPE}" = "tso" ]; then
        REPORT_TITLE="TSO500"
    elif [ "${REPORT_TYPE}" = "standard" ]; then 
        REPORT_TITLE=$(echo ${REPORT_TYPE} | sed 's/^s/S/')
    elif [ "${REPORT_TYPE}" = "advanced" ]; then 
        REPORT_TITLE=$(echo ${REPORT_TYPE} | sed 's/^a/A/')
    fi
    ps2pdf rm_cover_paged_add_cover_report.pdf ${ORDER_ID}.${REPORT_TITLE}_Analysis_Report.pdf
    # Copy Final Report to Output Folder 
    cp ${ORDER_ID}.${REPORT_TITLE}_Analysis_Report.pdf ../
    # Remove Temp Folder
    cd ../
    if [ "${REMOVE_TMP}" = "ture" ]; then rm -r tmp; fi
#---------------------------------------------------------------------------------------------------




# report_dir="/mnt/d/NGS_Pipeline_Setup/WTS_MAQC"
# order_id="GCX-C01-250101"

# # remove first page (cover page)
# report_pdf=${report_dir}/WTS_analysis_sample_report.pdf
# cover_rm_report_pdf=${report_dir}/cover_rm_sample_report.pdf

# pdftk ${report_pdf} cat 2-end output ${cover_rm_report_pdf}

# # count pdf pages
# report_page_n=$(pdftk ${cover_rm_report_pdf} dump_data | grep NumberOfPages | awk '{print $2}')


# # prepare stamp tex
# cp default_stamp.tex report_stamp.tex
# cp default_stamp.tex ${report_dir}/report_stamp.tex

# sed -i -e "s/ORDER_ID_HERE/${order_id}/" ${report_dir}/report_stamp.tex
# sed -i -e "s/PDF_PAGE_HERE/${report_page_n}/" ${report_dir}/report_stamp.tex

# stamp_name=${order_id}_stamp

# pdflatex -output-format "pdf_document" -jobname ${stamp_name} -output-directory ${report_dir} ${report_dir}/report_stamp.tex

# # stamping 
# pdftk ${cover_rm_report_pdf} multistamp ${report_dir}/${stamp_name}.pdf output ${report_dir}/paged_sample_report.pdf

# # cover merge 
# cover_pdf=${report_dir}/wts_analysis_report_cover_0.pdf

# pdftk ${cover_pdf} ${report_dir}/paged_sample_report.pdf cat output ${report_dir}/cover_merged_report.pdf

# ps2pdf ${report_dir}/cover_merged_report.pdf ${report_dir}/${order_id}_WTS_Analysis_Sample_Report.pdf 



