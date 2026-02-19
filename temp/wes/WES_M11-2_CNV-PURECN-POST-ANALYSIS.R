

# SCRIPT VERSION : v0.1
# GUIDELINE ID   : GMC-NGS-B02-A03
# DATE           : 2024-03-11
# AUTHOR         : kangsm@gencurix.com

#---| PACKAGES |-----------------------------------------------------------------------------------#
options(stringsAsFactors=FALSE)   
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("lubridate"))
#--------------------------------------------------------------------------------------------------#

#---| Set Defualt Values |-------------------------------------------------------------------------#
    baseDir     <- "/data/wes"
    #varCallMode <- "nt"
    resFolder   <- "Tonly"
#--------------------------------------------------------------------------------------------------#

#---| ARGUMENTS |----------------------------------------------------------------------------------#
    option_list = list( 
        make_option(c("--seqFolder"), action="store", default=NA, type="character", help="SEQ_FOLDER"),   
        make_option(c("--seqID"),     action="store", default=NA, type="character", help="SEQ_ID"),
        make_option(c("--baseDir")  , action="store", default=NA, type="character", help="$BASE_DIR"),
        make_option(c("--resultFolder"),    action="store", default="mutect2.NT", type="character", help="Result Folder Prefix")  
    )
    #--------------------------------------------------------------------------#  
    ARGS = parse_args(OptionParser(option_list=option_list))
    #--------------------------------------------------------------------------#
    seqFolder <- ARGS$seqFolder
    seqID     <- ARGS$seqID
    baseDir   <- ARGS$baseDir
    resFolder <- ARGS$resultFolder
#--------------------------------------------------------------------------------------------------#

#---| FOLDERS |------------------------------------------------------------------------------------#
    DATA_DIR        <- sprintf("%s/%s", baseDir, seqFolder)
    PURECN_RES_DIR  <- sprintf("%s/%s/cnv/purecn", DATA_DIR, seqID)
#--------------------------------------------------------------------------------------------------#

#---| READ PURECN RESULT |-------------------------------------------------------------------------#
    PURECN_DATA_DIR <- sprintf("%s/%s", PURECN_RES_DIR, resFolder)
    #--------------------------------------------------------------------------#
    purecnResultFile <- sprintf("%s/%s.csv", PURECN_DATA_DIR, seqID)
    if( file.exists(purecnResultFile) )
    {
        purecnRes  <- read.csv(purecnResultFile)
        PURITY     <- purecnRes$Purity
        PLOIDY     <- as.numeric(purecnRes$Ploidy)
        PLOIDY_INT <- round(as.numeric(purecnRes$Ploidy))
        #----------------------------------------------------------------------#
        purecnSummary <- data.frame(
            ID         = seqID,
            PURITY     = PURITY,
            PLOIDY     = round(PLOIDY, 2),
            PLOIDY_INT = PLOIDY_INT
        )
        #----------------------------------------------------------------------#
        purecnSummaryFileName <- sprintf("%s/%s.%s.PureCN.summary", PURECN_RES_DIR, seqID, resFolder)
        #----------------------------------------------------------------------#
        write.table( purecnSummary, sprintf("%s.tsv", purecnSummaryFileName), quote=F, col.names=T, row.names=F, sep="\t" )
        write.table( purecnSummary, sprintf("%s.txt", purecnSummaryFileName), quote=F, col.names=F, row.names=F, sep=" "  )
    }
#--------------------------------------------------------------------------------------------------#


