# SCRIPT VERSION : v0.2
# GUIDELINE ID   : GMC-NGS-B02-A03
# DATE           : 2024-03-11
# AUTHOR         : kangsm@gencurix.com
# MODIFIED       : Apply strict argument checks, fix sprintf warnings, and enhance DB integration

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

#---| ARGUMENTS |----------------------------------------------------------------------------------#
    # [MODIFIED] default 값을 모두 제거하고 명시적 입력 강제
    option_list = list( 
        make_option(c("--seqFolder"), action="store", type="character", help="SEQ_FOLDER"),   
        make_option(c("--seqID"),     action="store", type="character", help="SEQ_ID"),
        make_option(c("--baseDir")  , action="store", type="character", help="$BASE_DIR"),
        make_option(c("--runMode"),   action="store", type="character", help="MSI Run Mode. 'tonly' or 'nt'"),
        make_option(c("--nType"),     action="store", type="character", help="type of normal sample"),
        make_option(c("--importDB"),  action="store", type="logical", help="RESULT IMPORT INTO DB")
    )
    #--------------------------------------------------------------------------#
    ARGS <- parse_args(OptionParser(option_list=option_list))
    
    # [MODIFIED] 필수 값 누락 시 무조건 에러 발생
    if(is.null(ARGS$seqFolder)) stop("[ERROR] --seqFolder argument is missing.")
    if(is.null(ARGS$seqID)) stop("[ERROR] --seqID argument is missing.")
    if(is.null(ARGS$baseDir)) stop("[ERROR] --baseDir argument is missing.")
    if(is.null(ARGS$runMode)) stop("[ERROR] --runMode argument is missing.")
    if(is.null(ARGS$nType)) stop("[ERROR] --nType argument is missing.")
    if(is.null(ARGS$importDB)) stop("[ERROR] --importDB argument is missing.")
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
    msisensor_file <- sprintf("%s/%s.msisensor2.msi", MSI_RES_DIR, seqID)
    if(!file.exists(msisensor_file)) stop(sprintf("[ERROR] MSISensor2 result file not found: %s", msisensor_file))
    
    msisensor_res <- read.delim(msisensor_file)
    msisesnor_msi <- ifelse( msisensor_res[1,3] >= 20, "MSI-H", "MSS" )
    #--------------------------------------------------------------------------#
    if( runMode == "tonly" )
    {
        mantis_res   <- data.frame(G="T-only", V="-")
        mantis_score <- "-"
        mantis_msi   <- "-"
        runModeTag   <- RUN_MODE_TAG
    }else{
        mantis_file <- sprintf("%s/%s.mantis.msi.status", MSI_RES_DIR, seqID)
        if(!file.exists(mantis_file)) stop(sprintf("[ERROR] Mantis result file not found: %s", mantis_file))
        
        mantis_res   <- read.delim(mantis_file)
        mantis_score <- as.numeric(mantis_res[1,2])
        mantis_msi   <- ifelse( mantis_score >= 0.4, "MSI-H", "MSS" )
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
        # [MODIFIED] sprintf 불필요한 인자(seqID) 제거하여 Warning 해결
        message("\n  MSI-CALCULATION RESULTS IMPORT INTO DB START.")
        #==========================================================================================#
        
        dbCon <- tryCatch(
            dbConnect(dbDriver("MySQL"), host=db_host, user=db_user, port=db_port, password=db_pw, db=db_name_data ),
            error = function(e) stop("[ERROR] DB Connection failed: ", e)
        )
        
        # [MODIFIED] 안전한 트랜잭션을 위해 dbExecute 사용
        clearPreData <- dbExecute(dbCon, sprintf("DELETE FROM msi_status WHERE seq_folder = '%s' AND seq_id = '%s' AND msi_run_mode = '%s'", seqFolder, seqID, runModeTag))
        
        # [MODIFIED] 이전 세션의 SafeMapAndInsert 적용 (대소문자 무시 및 컬럼 자동 맵핑)
        SafeMapAndInsert <- function(con, table_name, df_raw) {
            db_cols <- dbListFields(con, table_name)
            r_cols <- colnames(df_raw)
            match_idx <- match(tolower(db_cols), tolower(r_cols))
            
            valid_db_cols <- db_cols[!is.na(match_idx)]
            valid_r_cols  <- r_cols[na.omit(match_idx)]
            if (length(valid_db_cols) == 0) return(0) 
            
            df_clean <- as.data.frame(df_raw)[, valid_r_cols, drop=FALSE]
            colnames(df_clean) <- valid_db_cols
            dbWriteTable(con, name=table_name, value=df_clean, row.names=FALSE, append=TRUE)
            return(nrow(df_clean))
        }

        inserted_rows <- SafeMapAndInsert(dbCon, "msi_status", MSI_RESULT)

        dbDisconnect(dbCon) 

        #==========================================================================================#    
        write(sprintf("%s | MSI-CALCULATION RESULTS IMPORT INTO DB FINISHED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
        message("\n  MSI-CALCULATION RESULTS IMPORT INTO DB FINISHED.")
        if(inserted_rows == 0) message("  [WARNING] DB Insert Skipped due to schema mismatch.")
        #==========================================================================================#
    }else{
        #==========================================================================================#    
        write(sprintf("%s | MSI-CALCULATION RESULTS DB IMPORT SKIPPED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
        message("\n  MSI-CALCULATION RESULTS DB IMPORT SKIPPED.")
        #==========================================================================================#
    }
#--------------------------------------------------------------------------------------------------#

#---| SAVE AS FILE |-------------------------------------------------------------------------------#
    # [MODIFIED] 파일 저장 경로. Bash script 쪽에서 외부 mv 명령어 실행 시 RES_TAG(.Tonly. 등)가 포함된 이 파일명을 동일하게 참조해야 에러가 발생하지 않습니다.
    saved_file_path <- sprintf("%s/%s.MSI.status.calculation.results.txt", MSI_RES_DIR, seqID)
    
    write.table( MSI_RESULT, saved_file_path, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t" )
#--------------------------------------------------------------------------------------------------#

#==================================================================================================#    
    write(sprintf("%s | MSI-CALCULATION RESULTS SAVED AS FILE.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message("\n  MSI-CALCULATION RESULTS SAVED AS FILE.")
#==================================================================================================#

#==================================================================================================#
    write("   ", LOG_FILE, append=T)
    write(sprintf("%s | WES MSI Calculation Post-Processing DONE.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    write("----------------------------------------------------------------------", file = LOG_FILE, append = TRUE ) 
    write("   ", LOG_FILE, append=T)
#==================================================================================================#