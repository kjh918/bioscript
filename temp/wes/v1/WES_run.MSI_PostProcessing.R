
# SCRIPT VERSION : v0.1
# GUIDELINE ID   : GMC-NGS-B02-A03
# DATE           : 2024-03-11
# AUTHOR         : kangsm@gencurix.com

#---| PACKAGES |-----------------------------------------------------------------------------------#
options(stringsAsFactors=FALSE)   
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Hmisc"))
suppressPackageStartupMessages(library("RMySQL"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("lubridate"))
#--------------------------------------------------------------------------------------------------#

#---| Set Defualt Values |-------------------------------------------------------------------------#
    baseDir <- "/data/wes"
#--------------------------------------------------------------------------------------------------#

#---| ARGUMENTS |----------------------------------------------------------------------------------#
    option_list = list( 
        make_option(c("--seqFolder"), action="store", default=NA, type="character", help="SEQ_FOLDER"),   
        make_option(c("--seqID"),     action="store", default=NA, type="character", help="SEQ_ID"),
        make_option(c("--baseDir")  , action="store", default=NA, type="character", help="$BASE_DIR"),
        make_option(c("--runMode"),   action="store", default=NA, type="character", help="MSI Run Mode. 'tonly' or 'nt'"),
        make_option(c("--nType"),     action="store", default=NA, type="character", help="type of normal sample"),
        make_option(c("--importDB"),  action="store", default=TRUE, type="logical", help="RESULT IMPORT INTO DB")
    )
    #--------------------------------------------------------------------------#
    ARGS <- parse_args(OptionParser(option_list=option_list))
    #--------------------------------------------------------------------------#
    seqFolder  <- ARGS$seqFolder
    seqID      <- ARGS$seqID
    baseDir    <- ARGS$baseDir
    runMode    <- ARGS$runMode
    normalType <- ARGS$nType
    IMPORT_DB  <- ARGS$importDB
#--------------------------------------------------------------------------------------------------#

#---| DATABASE CONNECTIONS |-----------------------------------------------------------------------#
    source("/data/wes/params/ruo_wes_db.R")
#--------------------------------------------------------------------------------------------------#

#---| FOLDERS |------------------------------------------------------------------------------------#
    DATA_DIR    <- sprintf("%s/%s", baseDir, seqFolder)
    LOG_DIR     <- sprintf("%s/%s/log", DATA_DIR, seqID)
    MSI_RES_DIR <- sprintf("%s/%s/msi", DATA_DIR, seqID)
#--------------------------------------------------------------------------------------------------#

#---| LOG FILE |-----------------------------------------------------------------------------------#
    LOG_FILE <- sprintf("%s/%s.MSI.post-processing.%s.log", LOG_DIR, seqID, today())
    if( ! file.exists(LOG_FILE) ){ system(sprintf("touch %s", LOG_FILE)) }
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    write(" ", LOG_FILE, append=T)
    write("WES MSI Calculation Post-Processing log", LOG_FILE, append=T)
    write("----------------------------------------------------------------------", LOG_FILE, append=T)
    write(sprintf(" RUN DATE   : %s", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T) 
    write(sprintf(" SEQ FOLDER : %s", seqFolder), LOG_FILE, append=T) 
    write(sprintf(" SEQ ID     : %s", seqID), LOG_FILE, append=T) 
    write(sprintf(" DB IMPOT   : %s", IMPORT_DB), LOG_FILE, append=T) 
    write("----------------------------------------------------------------------", LOG_FILE, append=T)
    write(" ", LOG_FILE, append=T)
#==================================================================================================# 

#==================================================================================================# 
    write(sprintf("%s | MSI Calculation Result Summary Start.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message("MSI Calculation Result Summary Start.")
#==================================================================================================# 

#---| RESULT TAG |---------------------------------------------------------------------------------#
    if(  runMode == "tonly" )
    {
        RUN_MODE_TAG = "Tumor-Only"
        RES_TAG      = "Tonly"
    }else{
        if( normalType == "ORG" )
        {
            RUN_MODE_TAG = "Matched-Normal-ORG"
            RES_TAG      = "ORG.NT"
        }else{
            RUN_MODE_TAG = "Matched-Normal"
            RES_TAG      = "NT"
        }
    }
#--------------------------------------------------------------------------------------------------#

#---| READ MSI RESULTS |---------------------------------------------------------------------------#
    msisensor_res <- read.delim(sprintf("%s/%s.msisensor2.msi",MSI_RES_DIR, seqID))
    msisesnor_msi <- ifelse( msisensor_res[1,3] >= 20, "MSI-H", "MSS" )
    #--------------------------------------------------------------------------#
    if( runMode == "tonly" )
    {
        mantis_res   <- data.frame(G="T-only", V="-")
        mantis_score <- "-"
        mantis_msi   <- "-"
        runModeTag   <- RUN_MODE_TAG
    }else{
        mantis_res   <- read.delim(sprintf("%s/%s.mantis.msi.status",MSI_RES_DIR, seqID))
        mantis_score <- mantis_res[1,2]
        mantis_msi   <- ifelse( mantis_res[1,2] >= 0.4, "MSI-H", "MSS" )
        runModeTag   <- RUN_MODE_TAG
    }
    #--------------------------------------------------------------------------#
    MSI_RESULT <- data.frame(
        seq_folder            = seqFolder,
        seq_id                = seqID,
        msi_run_mode          = runModeTag,
        mantis_score          = mantis_score,
        msisensor2_score      = msisensor_res[1,3],
        mantis_msi_status     = mantis_msi,
        msisensor2_msi_status = msisesnor_msi 
    )
#--------------------------------------------------------------------------------------------------#

#---| RESULT SAVE |--------------------------------------------------------------------------------#
    write(sprintf("%s | MSI Calculation Result Summary Finished.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    message("MSI Calculation Result Summary Finished.")
#==================================================================================================# 

#---| IMPORT INTO DB |-----------------------------------------------------------------------------#
    if( IMPORT_DB )
    {
        #==========================================================================================#    
        write(sprintf("%s | MSI-CALCULATION RESULTS IMPORT INTO DB START.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
        message(sprintf("\n  MSI-CALCULATION RESULTS IMPORT INTO DB START.", seqID))
        #==========================================================================================#
        
        dbCon        <- dbConnect(dbDriver("MySQL"), host=db_host, user=db_user, port=db_port, password=db_pw, db=db_name_data )
        clearPreData <- dbGetQuery(dbCon, sprintf("DELETE FROM msi_status WHERE seq_folder = '%s' AND seq_id = '%s' AND msi_run_mode = '%s'", seqFolder, seqID, runModeTag))
        writeNewData <- dbWriteTable(dbCon, name='msi_status', value=MSI_RESULT, row.names=F, append=T)
        dbDisconnect(dbCon) 

        #==========================================================================================#    
        write(sprintf("%s | MSI-CALCULATION RESULTS IMPORT INTO DB FINISHED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
        message(sprintf("\n  MSI-CALCULATION RESULTS IMPORT INTO DB FINISHED.", seqID))
        #==========================================================================================#
    }else{
        #==========================================================================================#    
        write(sprintf("%s | MSI-CALCULATION RESULTS DB IMPORT SKIPPED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
        message(sprintf("\n  MSI-CALCULATION RESULTS DB IMPORT SKIPPED."))
        #==========================================================================================#
    }
#--------------------------------------------------------------------------------------------------#

#---| SAVE AS FILE |-------------------------------------------------------------------------------#
    write.table( MSI_RESULT,
        sprintf("%s/%s.%s.MSI.status.calculation.results.txt", MSI_RES_DIR, seqID, RES_TAG ),
        quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t"
    )
#--------------------------------------------------------------------------------------------------#

#==================================================================================================#    
    write(sprintf("%s | MSI-CALCULATION RESULTS SAVED AS FILE.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message(sprintf("\n  MSI-CALCULATION RESULTS SAVED AS FILE."))
#==================================================================================================#

#==================================================================================================#
    write("   ", LOG_FILE, append=T)
    write(sprintf("%s | WES MSI Calculation Post-Processing DONE.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    write("----------------------------------------------------------------------", file = LOG_FILE, append = TRUE ) 
    write("   ", LOG_FILE, append=T)
#==================================================================================================# 

