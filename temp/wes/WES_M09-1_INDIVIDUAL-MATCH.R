
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
      make_option(c("--baseDir")  , action="store", default=NA, type="character", help="$BASE_DIR"),
      make_option(c("--assembly"),  action="store", default="hg19", type="character", help="Reference Genome Version"),
      make_option(c("--importDB"),  action="store", default=FALSE, type="logical", help="Result import database")
    )
    #--------------------------------------------------------------------------#  
    ARGS = parse_args(OptionParser(option_list=option_list))
    #--------------------------------------------------------------------------#
    seqFolder       <- ARGS$seqFolder
    baseDir         <- ARGS$baseDir
    genomeAessembly <- ARGS$assembly
    importDB        <- ARGS$importDB
    #loadFromDB      <- ARGS$sampleInfoFromDB
    #groupID         <- ARGS$sampleGroupID
#--------------------------------------------------------------------------------------------------#

#---| DATABASE CONNECTIONS |-----------------------------------------------------------------------#
    source("/data/wes/params/ruo_wes_db.R")
#--------------------------------------------------------------------------------------------------#

#---| FOLDERS |------------------------------------------------------------------------------------#
    DATA_DIR   <- sprintf("%s/%s", baseDir, seqFolder)
    LOG_DIR    <- sprintf("%s/meta/analysis_log", DATA_DIR)
    IM_RES_DIR <- sprintf("%s/meta/individual_matching", DATA_DIR)
    if( ! dir.exists(IM_RES_DIR) ){ system(sprintf("mkdir -p %s", IM_RES_DIR)) }
#--------------------------------------------------------------------------------------------------#

#---| LOGFILE |------------------------------------------------------------------------------------#
    LOG_FILE <- sprintf("%s/run.individual.matching.%s.log", LOG_DIR, today())
    if( ! file.exists(LOG_FILE) ){ system(sprintf("touch %s",LOG_FILE)) }
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    write(" ", LOG_FILE, append=T)
    write("WES Individual-Matching log", LOG_FILE, append=T)
    write("----------------------------------------------------------------------", LOG_FILE, append=T)
    write(sprintf(" RUN DATE   : %s", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T) 
    write(sprintf(" SEQ FOLDER : %s", seqFolder), LOG_FILE, append=T) 
    write(sprintf(" DB IMPORT  : %s", importDB), LOG_FILE, append=T) 
    write("----------------------------------------------------------------------", LOG_FILE, append=T)
    write(" ", LOG_FILE, append=T)
#==================================================================================================# 

#==================================================================================================# 
    write(sprintf("%s | Sample List of Group Creation Start.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message("Sample List of Group Creation Start.")
#==================================================================================================# 

#---| SAMPLE INFO |--------------------------------------------------------------------------------#       
        dbCon <- dbConnect(dbDriver("MySQL"), host=db_host, user=db_user, port=db_port, password=db_pw, db=db_name_info )
        sinfo <- dbGetQuery(dbCon, sprintf("SELECT * FROM seqid_info WHERE seq_folder = '%s'", seqFolder))
        dbDisconnect(dbCon)
        #----------------------------------------------------------------------#
        sinfo[is.na(sinfo)] <- ""
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    write(sprintf("%s | Run NGS-Checkmate Start.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message("Run NGS-Checkmate Start.")
#==================================================================================================# 

#---| RUN NGS-CHECKMATE by GROUP |-----------------------------------------------------------------#
    RUN_INDIVIDUAL_MATCHING <- "/storage/home/kangsm/runScripts/WES_run.NGSCheckMate_Rev.sh"
    SAMPLE_LIST_FILE        <- sprintf("%s/meta/sample.list", DATA_DIR)
    #--------------------------------------------------------------------------#
    RUN_CMD <- sprintf("%s --baseDir %s --seqFolder %s --sampleListFile %s --assembly %s --platform %s", RUN_INDIVIDUAL_MATCHING, baseDir, seqFolder, SAMPLE_LIST_FILE, genomeAessembly, "WES")
    #--------------------------------------------------------------------------#
    system(RUN_CMD)
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    write(" ", LOG_FILE, append=T)
    write(sprintf("%s | Run NGS-Checkmate Finished.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message("Run NGS-Checkmate Finished.")
#==================================================================================================# 

#==================================================================================================# 
    write(sprintf("%s | Individual-Matching Result Summary Start.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message("Individual-Matching Result Summary Start.")
#==================================================================================================# 

#---| INDIVIDUAL-MATCHING RESULT SUMMARY |---------------------------------------------------------#
    #MATCHING_RESULTS <- read.delim(sprintf("%s/%s_matched.txt", IM_RES_DIR, seqFolder), header=FALSE) 
    MATCHING_RESULTS <- read.delim(sprintf("%s/%s_all.txt", IM_RES_DIR, seqFolder), header=FALSE) 
    if( nrow(MATCHING_RESULTS) > 0 )
    {
        MATCHING_RESULTS <- MATCHING_RESULTS %>%
            dplyr::rename(seqid_1=1, match_res=2, seqid_2=3, cor=4, dp=5) %>%
            mutate( seqid_1 = gsub("\\.indiv.match.vcf", "", seqid_1) ) %>%
            mutate( seqid_2 = gsub("\\.indiv.match.vcf", "", seqid_2) ) %>%
            mutate( seq_folder = seqFolder, group="" )

        #==========================================================================================# 
            write(sprintf("%s | Individual-Matching Result Summary Finished.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
            message("Individual-Matching Result Summary Finished.")
        #==========================================================================================# 

        #---| INDIVIDUAL-MATCHING RESULT SAVE AS FILE |--------------------------------------------#
            write.table( MATCHING_RESULTS,
                sprintf("%s/%s.Individual.Matching.Results.txt", IM_RES_DIR, seqFolder),
                quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t"
            )
        #------------------------------------------------------------------------------------------#

        #==========================================================================================# 
            write(sprintf("%s | Individual-Matching Result Saved as File.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
            message("Individual-Matching Result Saved as File.")
        #==========================================================================================# 

        #---| RESULT IMPORT INTO DB |--------------------------------------------------------------#
            if( importDB )
            {
                dbCon <- dbConnect(dbDriver("MySQL"), host=db_host, user=db_user, port=db_port, password=db_pw, db=db_name_data )
                clearPreData <- dbGetQuery(dbCon, sprintf("DELETE FROM indiv_match WHERE seq_folder = '%s'", seqFolder))
                writeNewData <- dbWriteTable(dbCon, name='indiv_match', value=MATCHING_RESULTS, row.names=F, append=T)
                dbDisconnect(dbCon) 
                ##
                #==================================================================================#
                    write(sprintf("%s | Individual-Matching Result Import into DB Finished.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
                    message("Individual-Matching Result Import into DB Finished.")
                #==================================================================================#
            }else{
                #==========================================================================================#
                    write(sprintf("%s | Individual-Matching Result Import into DB Skipped.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
                    message("Individual-Matching Result Import into DB Skipped.")
                #==================================================================================#

            }
        #------------------------------------------------------------------------------------------#
    }
#--------------------------------------------------------------------------------------------------#




#==================================================================================================# 
    write(" ", LOG_FILE, append=T)
    write("WES Individual-Matching DONE.", LOG_FILE, append=T)
    write("----------------------------------------------------------------------", LOG_FILE, append=T)
    write(" ", LOG_FILE, append=T)
#==================================================================================================# 


