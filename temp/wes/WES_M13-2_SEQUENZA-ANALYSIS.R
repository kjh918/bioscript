
#---| PACKAGES |-----------------------------------------------------------------------------------#
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("sequenza"))
    suppressPackageStartupMessages(library("iotools"))
    suppressPackageStartupMessages(library("readr"))
    grDevices::X11.options(type='cairo')
    options(device='x11')
    options(stringsAsFactors=FALSE) 
    options(bitmapType='cairo')
#--------------------------------------------------------------------------------------------------#

#---| ARGUMENTS |----------------------------------------------------------------------------------#
    option_list = list( 
        make_option(c("--BASE_DIR"  ), action="store", default=NA, type="character", help="BASE DIR"),
        make_option(c("--SEQ_FOLDER"), action="store", default=NA, type="character", help="SEQ FOLDER"),   
        make_option(c("--SEQ_ID"    ), action="store", default=NA, type="character", help="TUMOR SEQ ID"),
        make_option(c("--RESULT_DIR"), action="store", default=NA, type="character", help="RESULT OUTPUT DIR"),
        make_option(c("--BIN_SIZE"  ), action="store", default=NA, type="character", help="SEQUENZA ANALYSIS BIN SIZE")
    )
    #--------------------------------------------------------------------------#  
    ARGS = parse_args(OptionParser(option_list=option_list))
    #--------------------------------------------------------------------------#
    BASE_DIR   <- ARGS$BASE_DIR
    SEQ_FOLDER <- ARGS$SEQ_FOLDER
    SEQ_ID     <- ARGS$SEQ_ID
    RESULT_DIR <- ARGS$RESULT_DIR
    BIN_SIZE   <- ARGS$BIN_SIZE
#--------------------------------------------------------------------------------------------------#

#---| INPUT CHECK |--------------------------------------------------------------------------------#
    ## path & sample-id check
    if( is.na(BASE_DIR) | is.na(SEQ_FOLDER) | is.na(SEQ_ID) ){ stop("Some Required Values are NOT FOUND.") }
    ## default result dir = data dir ( if not specified )
    if( is.na(RESULT_DIR) ){ RESULT_DIR = sprintf("%s/%s/%s/cnv/sequenza", BASE_DIR, SEQ_FOLDER, SEQ_ID) }
#--------------------------------------------------------------------------------------------------#

#---| RUN SEQUENZA |-------------------------------------------------------------------------------#
    seqzResFile <- sprintf("%s/%s/%s/cnv/sequenza/%s_bin%s.seqz.gz", BASE_DIR, SEQ_FOLDER, SEQ_ID, SEQ_ID, BIN_SIZE)
    if( !file.exists(seqzResFile) ){ stop("No SEQZ FILE FOUND.")}
    #--------------------------------------------------------------------------#
    # 1. read and processing
    message(">>>>> Running sequenza.extract....")
    seqzRes <- sequenza.extract(seqzResFile, verbose=FALSE)
    message(">>>>> Running sequenza.extract....DONE.")
    # 2. fit and analysis
    message(">>>>> Running sequenza.fit....")
    seqResFit <- sequenza.fit(seqzRes)
    message(">>>>> Running sequenza.fit....DONE.")
    # 3. save results
    message(">>>>> Running sequenza.result....")
    sequenza.results(
        sequenza.extract = seqzRes,
        cp.table         = seqResFit, 
        sample.id        = SEQ_ID,
        out.dir          = RESULT_DIR
    )
    # 4. save Robject additionally
    RDATA_NAME <- sprintf("%s/%s/%s/cnv/sequenza/%s_sequenza.AnalysisRobject.Rdata", BASE_DIR, SEQ_FOLDER, SEQ_ID, SEQ_ID)
    save(seqzRes,seqResFit, file=RDATA_NAME )
    message(">>>>> Running sequenza.result....DONE.")
#--------------------------------------------------------------------------------------------------#


# BASE_DIR   = "/storage/home/kangsm/analysis/gcx.MRD"
# SEQ_FOLDER = "test_sequenza"
# SEQ_ID     = "CRC01_T"
# RESULT_DIR = sprintf("%s/%s/%s/cnv/sequenza", BASE_DIR, SEQ_FOLDER, SEQ_ID) 


