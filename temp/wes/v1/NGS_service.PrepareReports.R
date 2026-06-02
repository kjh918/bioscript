
#---| PACKAGES |-----------------------------------------------------------------------------------#
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("RMySQL"))
    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("parallel"))
#--------------------------------------------------------------------------------------------------#

#---| ARGUMENTS |----------------------------------------------------------------------------------#
    option_list = list(
        make_option(c("--BASE_DIR"  ), action="store", default="/data/wes", type="character", help="data base dir"),
        make_option(c("--SEQ_FOLDER"), action="store", default=NA, type="character", help="data seq-folder")
    )
    #--------------------------------------------------------------------------#
    ARGS = parse_args(OptionParser(option_list=option_list))
    #--------------------------------------------------------------------------#
    BASE_DIR   <- ARGS$BASE_DIR
    SEQ_FOLDER <- ARGS$SEQ_FOLDER
#--------------------------------------------------------------------------------------------------#

#---| FUNCTIONS |----------------------------------------------------------------------------------#
    indexing = function(X) { sapply(unique(X), function(y) list(y)) }
#--------------------------------------------------------------------------------------------------#

#---| CHECK VALUES |-------------------------------------------------------------------------------#
    message(" >>> check and prepare report files...")
    if( is.na(SEQ_FOLDER) ){ stop("No SEQ_FOLDER.")}
#--------------------------------------------------------------------------------------------------#

#---| FIXED VALUES |-------------------------------------------------------------------------------#
    REPORT_BASE_DIR="/storage/home/kangsm/shinyWeb"
#--------------------------------------------------------------------------------------------------#
    
#---| RUN |----------------------------------------------------------------------------------------#
    dbCon    <- dbConnect(dbDriver("MySQL"), host='192.168.0.34', user='gcx', port=3306, password='gencurix!!', db='gcx_ngs_service' )
    ORDER_ID <- dbGetQuery(dbCon, sprintf("SELECT * FROM seqfolder_info WHERE seq_folder = '%s'", SEQ_FOLDER))
    #sinfo    <- dbGetQuery(dbCon, sprintf("SELECT seq_folder,seq_id FROM seqid_info WHERE seq_folder = '%s'", SEQ_FOLDER))
    dbDisconnect(dbCon)
    #--------------------------------------------------------------------------#   
    tumorOnly_List <- read.table(sprintf("%s/%s/meta/sample.list", BASE_DIR, SEQ_FOLDER), header=F, sep=" ")[,1]
    
    TS_N_analysisList = sprintf("%s/%s/meta/matched_analysis.list", BASE_DIR, SEQ_FOLDER)
    if( file.exists(TS_N_analysisList) )
    {
        nt_TS_List <- read.table(sprintf("%s/%s/meta/matched_analysis.list", BASE_DIR, SEQ_FOLDER), header=F, sep=" ") %>% 
            dplyr::rename(seqid=2, mode=4) %>% filter(mode == "TS") %>% .$seqid
    }else{
        nt_TS_List <- NULL
    }

    ORG_N_analysisList = sprintf("%s/%s/meta/matched_analysis.list", BASE_DIR, SEQ_FOLDER)
    if( file.exists(TS_N_analysisList) )
    {
        nt_ORG_List <- read.table(sprintf("%s/%s/meta/matched_analysis.list", BASE_DIR, SEQ_FOLDER), header=F, sep=" ") %>% 
            dplyr::rename(seqid=2, mode=4) %>% filter(mode == "ORG") %>% .$seqid
    }else{
        nt_ORG_List <- NULL
    }
    #--------------------------------------------------------------------------#
    message(" >>> create tumor-only unfiltered-raw-variants data...")
    raw_variants_Tonly <- mclapply(indexing(tumorOnly_List), function(idx)
    {
        read.delim(sprintf("%s/%s/%s/vcf/%s.mutect2.filtered.vep.annotated.maf", BASE_DIR, SEQ_FOLDER, idx, idx)) %>%
            dplyr::select(c("seq_folder","seq_id","vkey","TUMOR_ALT_COUNT","TUMOR_DEPTH","DEPTH"))
    }, mc.cores = length(tumorOnly_List))
    #--------------------------------------------------------------------------#
    message(" >>> create matched-normal TS unfiltered-raw-variants data...")
    if( length(nt_TS_List) != 0 )
    {
        raw_variants_TS_N <- mclapply(indexing(nt_TS_List), function(idx)
        {
            read.delim(sprintf("%s/%s/%s/vcf/%s.mutect2.NT.filtered.vep.annotated.maf", BASE_DIR, SEQ_FOLDER, idx, idx)) %>%
                dplyr::select(c("seq_folder","seq_id","vkey","TUMOR_ALT_COUNT","TUMOR_DEPTH","DEPTH"))
        }, mc.cores = length(nt_TS_List))
    }else{
        raw_variants_TS_N <- mclapply(indexing(nt_TS_List), function(idx) data.frame() )
    }
    #--------------------------------------------------------------------------#
    message(" >>> create matched-normal ORG unfiltered-raw-variants data...")
    if( length(nt_ORG_List) != 0 )
    {
        raw_variants_ORG_N <- mclapply(indexing(nt_ORG_List), function(idx)
        {
            read.delim(sprintf("%s/%s/%s/vcf/%s.mutect2.ORG.NT.filtered.vep.annotated.maf", BASE_DIR, SEQ_FOLDER, idx, idx)) %>%
                dplyr::select(c("seq_folder","seq_id","vkey","TUMOR_ALT_COUNT","TUMOR_DEPTH","DEPTH"))
        }, mc.cores = length(nt_ORG_List))
    }else{
        raw_variants_ORG_N <- mclapply(indexing(nt_ORG_List), function(idx) data.frame() )
    }
    #--------------------------------------------------------------------------#
    save(raw_variants_Tonly,raw_variants_TS_N, raw_variants_ORG_N,
        file=sprintf("%s/%s/meta/%s.unfiltered.raw.variants.vaf.Rdata", BASE_DIR, SEQ_FOLDER, SEQ_FOLDER) )
#--------------------------------------------------------------------------------------------------#

#---| REPORT FOLDER & FILES |----------------------------------------------------------------------#
    message(" >>> create report folder and copy files to folder...")
    STANDARD_REPORT_DIR = sprintf("%s/STANDARD_REPORT/%s_%s", REPORT_BASE_DIR, ORDER_ID$ngs_order_id, SEQ_FOLDER)
    if( !dir.exists(STANDARD_REPORT_DIR) ){ system(sprintf("mkdir -p %s", STANDARD_REPORT_DIR)) }
    ADVANCED_REPORT_DIR = sprintf("%s/ADVANCED_REPORT/%s_%s", REPORT_BASE_DIR, ORDER_ID$ngs_order_id, SEQ_FOLDER)
    if( !dir.exists(ADVANCED_REPORT_DIR) ){ system(sprintf("mkdir -p %s", ADVANCED_REPORT_DIR)) }
    save(raw_variants_Tonly,raw_variants_TS_N, raw_variants_ORG_N,
        file=sprintf("%s/%s.unfiltered.raw.variants.vaf.Rdata", ADVANCED_REPORT_DIR, SEQ_FOLDER ) )
    #--------------------------------------------------------------------------#
    system(sprintf("ln -s %s/resources %s/resources", REPORT_BASE_DIR, STANDARD_REPORT_DIR))
    system(sprintf("ln -s %s/resources %s/resources", REPORT_BASE_DIR, ADVANCED_REPORT_DIR))
    #--------------------------------------------------------------------------#
    system(sprintf("cp %s/resources/DEFAULT_standard_analysis_report.Rmd %s/%s_standard_analysis_report.Rmd", REPORT_BASE_DIR, STANDARD_REPORT_DIR, ORDER_ID$ngs_order_id))
    system(sprintf("cp %s/resources/DEFAULT_advanced_analysis_report.Rmd %s/%s_advanced_analysis_report.Rmd", REPORT_BASE_DIR, ADVANCED_REPORT_DIR, ORDER_ID$ngs_order_id))

    CLIENT_ID   = paste0( "C", ORDER_ID$facility_id )
    CLIENT_NAME = gsub("_", "", ORDER_ID$client_name)
    system(sprintf("sed -i 's/DEFAULT_ORDER_ID/%s/g' %s/%s_standard_analysis_report.Rmd",    ORDER_ID$ngs_order_id, STANDARD_REPORT_DIR, ORDER_ID$ngs_order_id))
    system(sprintf("sed -i 's/DEFAULT_SEQ_FOLDER/%s/g' %s/%s_standard_analysis_report.Rmd",  SEQ_FOLDER,            STANDARD_REPORT_DIR, ORDER_ID$ngs_order_id))
    system(sprintf("sed -i 's/DEFAULT_CLIENT_ID/%s/g' %s/%s_standard_analysis_report.Rmd",   CLIENT_ID,             STANDARD_REPORT_DIR, ORDER_ID$ngs_order_id))
    system(sprintf("sed -i 's/DEFAULT_CLIENT_NAME/%s/g' %s/%s_standard_analysis_report.Rmd", CLIENT_NAME,           STANDARD_REPORT_DIR, ORDER_ID$ngs_order_id))

    system(sprintf("sed -i 's/DEFAULT_ORDER_ID/%s/g' %s/%s_advanced_analysis_report.Rmd",    ORDER_ID$ngs_order_id, ADVANCED_REPORT_DIR, ORDER_ID$ngs_order_id))
    system(sprintf("sed -i 's/DEFAULT_SEQ_FOLDER/%s/g' %s/%s_advanced_analysis_report.Rmd",  SEQ_FOLDER,            ADVANCED_REPORT_DIR, ORDER_ID$ngs_order_id))
    system(sprintf("sed -i 's/DEFAULT_CLIENT_ID/%s/g' %s/%s_advanced_analysis_report.Rmd",   CLIENT_ID,             ADVANCED_REPORT_DIR, ORDER_ID$ngs_order_id))

    # if( !dir.exists(STANDARD_REPORT_DIR) )
    # {
    #     system(sprintf("mkdir -p %s", STANDARD_REPORT_DIR))
    #     #system(sprintf("ln -s %s/resources %s/resources", REPORT_BASE_DIR, STANDARD_REPORT_DIR))
    #     #system(sprintf("cp %s/resources/default_standard_analysis_report.Rmd %s/", REPORT_BASE_DIR, STANDARD_REPORT_DIR))
    #     system(sprintf("mv %s/default_standard_analysis_report.Rmd %s/%s_standard_analysis_report.Rmd", STANDARD_REPORT_DIR, STANDARD_REPORT_DIR, ORDER_ID$ngs_order_id))
    # }
    # #--------------------------------------------------------------------------#
    # if( !dir.exists(ADVANCED_REPORT_DIR) )
    # {
    #     system(sprintf("mkdir -p %s", ADVANCED_REPORT_DIR))
    #     system(sprintf("ln -s %s/resources %s/resources", REPORT_BASE_DIR, ADVANCED_REPORT_DIR))
    #     system(sprintf("cp %s/resources/default_advanced_analysis_report.Rmd %s/", REPORT_BASE_DIR, ADVANCED_REPORT_DIR))
    #     system(sprintf("mv %s/default_advanced_analysis_report.Rmd %s/%s_advanced_analysis_report.Rmd", STANDARD_REPORT_DIR, STANDARD_REPORT_DIR, ORDER_ID$ngs_order_id))
    #     save(raw_variants_Tonly,raw_variants_TS_N,raw_variants_TS_N,
    #         file=sprintf("%s/%s.unfiltered.raw.variants.vaf.Rdata", ADVANCED_REPORT_DIR, SEQ_FOLDER) )
    # }
    message(" >>> REPORT FOLDERS AND FILES PREPARED.")
#--------------------------------------------------------------------------------------------------#


