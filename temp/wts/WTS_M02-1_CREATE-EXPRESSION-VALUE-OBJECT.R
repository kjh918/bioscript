
# RNAseq BAM-QC RESULT INTEGRATION AND SAVE SCRIPT

#---| PACKAGES |-------------------------------------------------------------------------------------------------------------------------------------#
    options(stringsAsFactors=FALSE)   
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("plyr"))
    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("Hmisc"))
    suppressPackageStartupMessages(library("RMySQL"))
    suppressPackageStartupMessages(library("reshape2"))
    suppressPackageStartupMessages(library("yaml"))
#----------------------------------------------------------------------------------------------------------------------------------------------------#

#---| PARSE ARGS |-----------------------------------------------------------------------------------------------------------------------------------#
    option_list = list( 
        make_option(c("--BASE_DIR"),        action="store", default=NA,    type="character", help="BASE_DIR"),   
        make_option(c("--SEQ_FOLDER"),      action="store", default=NA,    type="character", help="SEQ_FOLDER"),
        make_option(c("--SEQ_ID"),          action="store", default=NA,    type="character", help="SEQ_ID"),
        make_option(c("--ANALYSIS_CONFIG"), action="store", default=NA,    type="character", help="analysis config YAML file. REQUIRED."),
        make_option(c("--DB_IMPORT"),       action="store", default=FALSE, type="logical",   help="DB_IMPORT")
    )
    #---------------------------------------------------------------------------  
    ARGS = parse_args(OptionParser(option_list=option_list))
    #---------------------------------------------------------------------------
    BASE_DIR        <- ARGS$BASE_DIR
    SEQ_FOLDER      <- ARGS$SEQ_FOLDER
    SEQ_ID          <- ARGS$SEQ_ID
    ANALYSIS_CONFIG <- ARGS$ANALYSIS_CONFIG
    DB_IMPORT       <- ARGS$DB_IMPORT
    #---------------------------------------------------------------------------
    if( is.na(ANALYSIS_CONFIG) )
    { stop("|---!!! wts analysis config yaml file is not found. REQUIRD. please check again. STOPPED.") }
#----------------------------------------------------------------------------------------------------------------------------------------------------#

#---| DATABASE CONNECTION |--------------------------------------------------------------------------------------------------------------------------#
    source("/data/wts/params/ruo_wts_db.R")
#----------------------------------------------------------------------------------------------------------------------------------------------------#

#---| READ ANALYSIS CONFIG YAML FILE |---------------------------------------------------------------------------------------------------------------#
    CONFIG <- read_yaml(ANALYSIS_CONFIG)
    if( is.null( CONFIG$GLOBAL$RDS_DIR ) )
    { stop("|---!!! no global RDS path found. REQUIRED. please check again. STOPPED.") }
#----------------------------------------------------------------------------------------------------------------------------------------------------#

#---| PRE-DEFINED FOLDERS |--------------------------------------------------------------------------------------------------------------------------#
    readcount_rds_dir <- sprintf("%s/RDS_ReadCount", CONFIG$GLOBAL$RDS_DIR)
    tpm_rds_dir       <- sprintf("%s/RDS_TPM",       CONFIG$GLOBAL$RDS_DIR)
    fpkm_rds_dir      <- sprintf("%s/RDS_FPKM",      CONFIG$GLOBAL$RDS_DIR)
#----------------------------------------------------------------------------------------------------------------------------------------------------#

#---| READ EXPRESSION VALUES |-----------------------------------------------------------------------------------------------------------------------#
    rsem_res_dir  <- sprintf("%s/%s/%s/quant/rsem", BASE_DIR, SEQ_FOLDER, SEQ_ID)
    rsem_gene_res <- read.delim(sprintf("%s/%s.genes.results", rsem_res_dir, SEQ_ID), header=TRUE)
    #rsem_gene_res$ens_geneid <- sapply( rsem_gene_res$gene_id, function(y) unlist(strsplit(y, "\\."))[1] ) 
    #---------------------------------------------------------------------------
    rsem_gene_res <- rsem_gene_res %>% mutate( 
        named_count = paste(gene_id, expected_count, sep=","),
        named_tpm   = paste(gene_id, TPM, sep=","),
        named_fpkm  = paste(gene_id, FPKM, sep=","),
        named_eff_length  = paste(gene_id, effective_length, sep=","),
        named_gene_length = paste(gene_id, length, sep=","),
        named_tx_id = paste(gene_id, transcript_id.s., sep=":")
    )
    #---------------------------------------------------------------------------
    Values_Counts <- paste(rsem_gene_res$named_count, collapse=";")
    Values_TPM    <- paste(rsem_gene_res$named_tpm,   collapse=";")
    Values_TPKM   <- paste(rsem_gene_res$named_fpkm,  collapse=";")
    Values_EffLength  <- paste(rsem_gene_res$named_eff_length,   collapse=";")
    Values_GeneLength <- paste(rsem_gene_res$named_gene_length,  collapse=";")
    Values_Tx         <- paste(rsem_gene_res$named_tx_id,        collapse=";")

    # expression values table (for database table) -----------------------------
    EXPR_VALUES_TABLE <- data.frame(
        seq_folder          = SEQ_FOLDER,
        seq_id              = SEQ_ID,
        swtools             = 'rsem',
        read_count_rds_path = sprintf("%s/%s.ReadCount.rds", readcount_rds_dir, SEQ_ID),
        tpm_rds_path        = sprintf("%s/%s.TPM.rds",       tpm_rds_dir, SEQ_ID      ),
        fpkm_rds_path       = sprintf("%s/%s.FPKM.rds",      fpkm_rds_dir, SEQ_ID     ),
        read_count          = Values_Counts,
        tpm                 = Values_TPM,
        fpkm                = Values_TPKM,
        eff_length          = Values_EffLength,
        gene_length         = Values_GeneLength,
        tx_id               = Values_Tx
    )

    # read_count object --------------------------------------------------------
    Obj_ReadCount        <- list( setNames( object = rsem_gene_res$expected_count, nm = rsem_gene_res$gene_id ) )
    names(Obj_ReadCount) <- SEQ_ID

    # TPM object ---------------------------------------------------------------
    Obj_TPM        <- list( setNames( object = rsem_gene_res$TPM, nm = rsem_gene_res$gene_id ) )
    names(Obj_TPM) <- SEQ_ID

    # read_count object --------------------------------------------------------
    Obj_FPKM        <- list( setNames( object = rsem_gene_res$FPKM, nm = rsem_gene_res$gene_id ) )
    names(Obj_FPKM) <- SEQ_ID
#----------------------------------------------------------------------------------------------------------------------------------------------------#

#---| SAVE RDS |-------------------------------------------------------------------------------------------------------------------------------------#
    saveRDS( Obj_ReadCount, file = sprintf("%s/%s.ReadCount.rds", readcount_rds_dir, SEQ_ID))
    saveRDS( Obj_TPM,       file = sprintf("%s/%s.TPM.rds",       tpm_rds_dir,       SEQ_ID))
    saveRDS( Obj_FPKM,      file = sprintf("%s/%s.FPKM.rds",      fpkm_rds_dir,      SEQ_ID))
#----------------------------------------------------------------------------------------------------------------------------------------------------#


#---| IMPORT INTO DATABASE |-------------------------------------------------------------------------------------------------------------------------#
    if( DB_IMPORT )
    {
        dbCon <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db=DB_NAME_DATA )
        clearPreExistData <- dbGetQuery(dbCon, sprintf("DELETE FROM expr_values WHERE seq_folder = '%s' AND seq_id = '%s'", SEQ_FOLDER, SEQ_ID))
        writeNewData      <- dbWriteTable(dbCon, name='expr_values', value=EXPR_VALUES_TABLE, row.names=FALSE, append=TRUE)
        dbDisconnect(dbCon)
    }
#----------------------------------------------------------------------------------------------------------------------------------------------------#
