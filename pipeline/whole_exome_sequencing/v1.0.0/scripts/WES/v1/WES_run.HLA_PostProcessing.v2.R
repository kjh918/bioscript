

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
      make_option(c("--importDB"),  action="store", default=TRUE, type="logical", help="RESULT IMPORT INTO DB")
    )
    #--------------------------------------------------------------------------#  
    ARGS = parse_args(OptionParser(option_list=option_list))
    #--------------------------------------------------------------------------#
    seqFolder <- ARGS$seqFolder
    seqID     <- ARGS$seqID
    baseDir   <- ARGS$baseDir
    IMPORT_DB <- ARGS$importDB
#--------------------------------------------------------------------------------------------------#

#---| DATABASE CONNECTIONS |-----------------------------------------------------------------------#
    source("/data/wes/params/ruo_wes_db.R")
#--------------------------------------------------------------------------------------------------#

#---| FOLDERS |------------------------------------------------------------------------------------#
    DATA_DIR    <- sprintf("%s/%s", baseDir, seqFolder)
    LOG_DIR     <- sprintf("%s/%s/log", DATA_DIR, seqID)
    HLA_RES_DIR <- sprintf("%s/%s/hla", DATA_DIR, seqID)
#--------------------------------------------------------------------------------------------------#

#---| LOG FILE |-----------------------------------------------------------------------------------#
    LOG_FILE <- sprintf("%s/%s.HLA.post-processing.%s.log", LOG_DIR, seqID, today())
    if( ! file.exists(LOG_FILE) ){ system(sprintf("touch %s", LOG_FILE)) }
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    write(" ", LOG_FILE, append=T)
    write("WES HLA-TYPING Post-Processing log", LOG_FILE, append=T)
    write("----------------------------------------------------------------------", LOG_FILE, append=T)
    write(sprintf(" RUN DATE   : %s", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T) 
    write(sprintf(" SEQ FOLDER : %s", seqFolder), LOG_FILE, append=T) 
    write(sprintf(" SEQ ID     : %s", seqID), LOG_FILE, append=T) 
    write(sprintf(" DB IMPOT   : %s", IMPORT_DB), LOG_FILE, append=T) 
    write("----------------------------------------------------------------------", LOG_FILE, append=T)
    write(" ", LOG_FILE, append=T)
#==================================================================================================# 

#---| READ HLA-TYPING RESULTS |--------------------------------------------------------------------#
    # OptiType
    optitypeResult <- sprintf("%s/optitype/%s.optitype_result.tsv",  HLA_RES_DIR, seqID)
    opti_res       <- read.delim(optitypeResult)
    # HLA-LA
    hlalaResult <- sprintf("%s/hla-la/%s/hla/R1_bestguess_G.txt", HLA_RES_DIR, seqID)
    hla_res     <- read.delim(hlalaResult)
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    write(sprintf("%s | HLA-Typing Results Read.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message("HLA-Typing Results Read.")
#==================================================================================================# 

#---| MHC-I RESULT SUMMARY |-----------------------------------------------------------------------#
    RES_major = data.frame()
    ##
    for( loc in c("A","B","C") ) {
        loc_opti <- sort(c(opti_res[1, sprintf("%s1", loc)], opti_res[1, sprintf("%s2", loc)]))
        loc_hla  <- sort(hla_res %>% filter( Locus == loc ) %>% .$Allele)
        loc_hla2 <- sapply(loc_hla, function(y) {
            z = unlist(strsplit(y, ";"))
            v = sapply(z, function(w) paste(unlist(strsplit(w, "\\:"))[1:2], collapse=":"))
            k = paste(v, collapse=";")
            return(k)
        })
        ##
        res_major <- data.frame(
            seq_folder       = seqFolder, 
            seq_id           = seqID, 
            concordance      = "", 
            report           = 1, 
            locus            = loc,
            allele1_optitype = loc_opti[1],  
            allele2_optitype = loc_opti[2],
            allele1_hla_la   = loc_hla2[1],  
            allele2_hla_la   = loc_hla2[2],
            allele1_hla_full = loc_hla[1],   
            allele2_hla_full = loc_hla[2]
        )
        ##
        RES_major = rbind(RES_major, res_major)
    }
    ##
    RES_major$concordance = apply(RES_major, 1, function(y) ifelse( y[6] == y[8] & y[7] == y[9], "Yes", "No" ) )
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    write(sprintf("%s | MHC-I Result Summary Done.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message("MHC-I Result Summary Done.")
#==================================================================================================# 

#---| MHC-II RESULT SUMMARY |----------------------------------------------------------------------#
    RES_minor = data.frame()
    ##
    for(loc2 in unique(hla_res$Locus[ hla_res$Locus %nin%  c("A","B","C") ])) {
        loc2_hla  <- sort(hla_res %>% filter( Locus == loc2 ) %>% .$Allele)
        loc2_hla2 <- sapply(loc2_hla, function(y) {
            z = unlist(strsplit(y, ";"))
            v = sapply(z, function(w) paste(unlist(strsplit(w, "\\:"))[1:2], collapse=":"))
            k = paste(v, collapse=";")
            return(k)
        })
        ##
        res_minor = data.frame(
            seq_folder       = seqFolder, 
            seq_id           = seqID, 
            concordance      = "", 
            report           = 0, 
            locus            = loc2, 
            allele1_optitype = "",  
            allele2_optitype = "",
            allele1_hla_la   = loc2_hla2[1],  
            allele2_hla_la   = loc2_hla2[2],
            allele1_hla_full = loc2_hla[1],   
            allele2_hla_full = loc2_hla[2]
        )
        ##
        RES_minor = rbind(RES_minor, res_minor)
    }
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    write(sprintf("%s | MHC-II Result Summary Done.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message("MHC-II Result Summary Done.")
#==================================================================================================# 

#---| COMBINE MHC-I & MHC-II RESULTS |-------------------------------------------------------------#
    HLA_RESULTS = rbind(RES_major, RES_minor)
#--------------------------------------------------------------------------------------------------#

#---| CHECK CONCORDANCE between SWTOOL's RESULTS |-------------------------------------------------#
    matchResult = HLA_RESULTS %>% filter( report == 1 ) %>% .$concordance %>% unique 
    if( length(matchResult) == 1 && matchResult == 'Yes' )
    {
        #==========================================================================================# 
        message(sprintf("\n  SEQ_ID = %s : ALL MHC-I HLA-TYPE matched between OptiType and HLA-LA)", seqID)) 
        write("----------------------------------------------------------------------", file = LOG_FILE, append = TRUE )
        write(sprintf("SEQ_ID = %s : ALL MHC-I HLA-TYPE matched between OptiType and HLA-LA", seqID), file = LOG_FILE, append = TRUE)
        write("----------------------------------------------------------------------", file = LOG_FILE, append = TRUE ) 
        #==========================================================================================# 
    }else{
        #==========================================================================================# 
        message(sprintf("\n  SEQ_ID = %s : more than 1 allele(s) of MHC-I HLA-TYPE NOT matched", seqID))
        write("----------------------------------------------------------------------", file = LOG_FILE, append = TRUE )
        write(sprintf("SEQ_ID = %s : more than 1 allele(s) of MHC-I HLA-TYPE NOT matched", seqID), file = LOG_FILE,  append = TRUE)
        write("----------------------------------------------------------------------", file = LOG_FILE, append = TRUE )
        #==========================================================================================# 
    }
#--------------------------------------------------------------------------------------------------#

#---| IMPORT INTO DB |-----------------------------------------------------------------------------#
    if( IMPORT_DB )
    {
        #==========================================================================================#    
        write(sprintf("%s | HLA-TYPING RESULTS IMPORT INTO DB START.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
        message(sprintf("\n  HLA-TYPING RESULTS IMPORT INTO DB START."))
        #==========================================================================================#
        
        dbCon <- dbConnect(dbDriver("MySQL"), host=db_host, user=db_user, port=db_port, password=db_pw, db=db_name_data )
        clearPreData <- dbGetQuery(dbCon, sprintf("DELETE FROM hla_type WHERE seq_folder = '%s' AND seq_id = '%s'", seqFolder, seqID))
        writeNewData <- dbWriteTable(dbCon, name='hla_type', value=HLA_RESULTS, row.names=F, append=T)
        dbDisconnect(dbCon) 

        #==========================================================================================#    
        write(sprintf("%s | HLA-TYPING RESULTS IMPORT INTO DB FINISHED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
        message(sprintf("\n  HLA-TYPING RESULTS IMPORT INTO DB FINISHED."))
        #==========================================================================================#
    }else{
        #==========================================================================================#    
        write(sprintf("%s | HLA-TYPING RESULTS DB IMPORT SKIPPED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
        message(sprintf("\n  HLA-TYPING RESULTS DB IMPORT SKIPPED."))
        #==========================================================================================#
    }
#--------------------------------------------------------------------------------------------------#

#---| SAVE AS FILE |-------------------------------------------------------------------------------#
    write.table( HLA_RESULTS,
        sprintf("%s/%s.HLA.typing.results.txt", HLA_RES_DIR, seqID ),
        quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t"
    )
#--------------------------------------------------------------------------------------------------#

#==================================================================================================#    
    write(sprintf("%s | HLA-TYPING RESULTS SAVED AS FILE.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message(sprintf("\n  HLA-TYPING RESULTS SAVED AS FILE."))
#==================================================================================================#
       
#==================================================================================================#
    write("   ", LOG_FILE, append=T)
    write(sprintf("%s | WES HLA-TYPING Post-Processing DONE.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    write("----------------------------------------------------------------------", file = LOG_FILE, append = TRUE ) 
    write("   ", LOG_FILE, append=T)
#==================================================================================================# 
    
    
