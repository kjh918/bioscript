
options(stringsAsFactors=FALSE)   
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Hmisc"))
suppressPackageStartupMessages(library("RMySQL"))
suppressPackageStartupMessages(library("reshape2"))

#---| PARSE ARGS |
    option_list = list( make_option(c("--seqFolder"), action="store", default=NA, type="character", help="SEQ_FOLDER")    ),
    option_list = list( make_option(c("--seqID"),     action="store", default=NA, type="character", help="SEQ_ID")        ),
    option_list = list( make_option(c("--logFile"),   action="store", default=NA, type="character", help="LOG FILE NAME") )    
  ##  
    ARGS = parse_args(OptionParser(option_list=option_list))
    SEQ_FOLDER = ARGS$seqFolder
    SEQ_ID     = ARGS$seqID
    LOG_NAME   = ARGS$logFile
  ##
    DATA_DIR = sprintf("/data/wes/%s",SEQ_FOLDER)
    LOGFILE  = sprintf("/data/wes/%s/%s/%S", SEQ_FOLDER, SEQ_ID, LOG_NAME)
#---/
#---| OPTITYPE & HLA-LA RESULT ANALYSIS |
#---\
    # message(sprintf("[ %s ] HLA typing result summarization start.", date()))
    # write(sprintf("[ %s ] HLA typing result summarization start.", date()), file=LogFileName, append=T)
    # message("--------------------------------------------------------------------------------------")
    # write("--------------------------------------------------------------------------------------", file=LogFileName, append=T )
    # message(sprintf("  Seq_Folder = %s", SEQ_FOLDER))
    # write(sprintf("  Seq_Folder = %s", SEQ_FOLDER), file=LogFileName, append=T)
#---/
  #-| READ RESULTS
    optitypeResult = sprintf("%s/%s/hla/optitype/%s.optitype_result.tsv", DATA_DIR, SEQ_ID, SEQ_ID)
    hlalaResult    = sprintf("%s/%s/hla/hla-la/%s/hla/R1_bestguess_G.txt", DATA_DIR, SEQ_ID, SEQ_ID)
  ## 
    opti_res = read.delim(optitypeResult)
    hla_res  = read.delim(hlalaResult)
  #-| MHC-I
    RES_major = data.frame()
    for( loc in c("A","B","C") ) {
        loc_opti = sort(c(opti_res[1, sprintf("%s1", loc)], opti_res[1, sprintf("%s2", loc)]))
        loc_hla  = sort(hla_res %>% filter( Locus == loc ) %>% .$Allele)
        loc_hla2 = sapply(loc_hla, function(y) {
            z = unlist(strsplit(y, ";"))
            v = sapply(z, function(w) paste(unlist(strsplit(w, "\\:"))[1:2], collapse=":"))
            k = paste(v, collapse=";")
            return(k)
        })
        res_major = data.frame(
            seq_folder = SEQ_FOLDER, seq_id = SEQ_ID, concordance = "", report = 1, 
            locus = loc,
            allele1_optitype = loc_opti[1],  allele2_optitype = loc_opti[2],
            allele1_hla_la   = loc_hla2[1],  allele2_hla_la   = loc_hla2[2],
            allele1_hla_full = loc_hla[1],   allele2_hla_full = loc_hla[2]
        )
        RES_major = rbind(RES_major, res_major)
    }
    RES_major$concordance = apply(RES_major, 1, function(y) ifelse( y[6] == y[8] & y[7] == y[9], "Yes", "No" ) )
  #-| MHC-II
    RES_minor = data.frame()
    for(loc2 in unique(hla_res$Locus[ hla_res$Locus %nin%  c("A","B","C") ])) {
        loc2_hla  = sort(hla_res %>% filter( Locus == loc2 ) %>% .$Allele)
        loc2_hla2 = sapply(loc2_hla, function(y) {
            z = unlist(strsplit(y, ";"))
            v = sapply(z, function(w) paste(unlist(strsplit(w, "\\:"))[1:2], collapse=":"))
            k = paste(v, collapse=";")
            return(k)
        })
        res_minor = data.frame(
            seq_folder = SEQ_FOLDER, seq_id = SEQ_ID, concordance = "", report = 0, 
            locus = loc2, allele1_optitype = "",  allele2_optitype = "",
            allele1_hla_la   = loc2_hla2[1],  allele2_hla_la   = loc2_hla2[2],
            allele1_hla_full = loc2_hla[1],   allele2_hla_full = loc2_hla[2]
        )
        RES_minor = rbind(RES_minor, res_minor)
    }
  ##  
    RESULTS = rbind(RES_major, RES_minor)
#---/
#---| CHECK CONCORDANCE between SWTOOL's RESULTS
    matchResult = RESULTS %>% filter( report == 1 ) %>% .$concordance %>% unique 
    if( length(matchResult) == 1 & matchResult == 'Yes' )
    {
        message(sprintf("\n  SeqID = %s : all HLA-type results matched. (OptiType and HLA-LA)", SEQ_ID)) 
        write( "   ", file = LOGFILE, append = TRUE )
        write(sprintf("SeqID = %s : all HLA-type results matched. (OptiType and HLA-LA)", SEQ_ID), file = LOGFILE, append = TRUE)
        write( "   ", file = LOGFILE, append = TRUE ) 
    }else{
        message(sprintf("\n  SeqID = %s : more than 1 allele(s) of HLA-type NOT matched", SEQ_ID))
        write( "   ", file = LOGFILE, append = TRUE )
        write(sprintf("SeqID = %s : more than 1 allele(s) of HLA-type NOT matched", SEQ_ID), file = LOGFILE,  append = TRUE)
        write( "   ", file = LOGFILE, append = TRUE )
    }
#---/    
#---| RESULT DB IMPORT    
    dbCon = dbConnect(dbDriver("MySQL"), host='192.168.0.35', user='gcx', port=3306, password='gencurix!!', db='gcx_wes' )
    clearPreExistData = dbGetQuery(dbCon, sprintf("DELETE FROM hla_typing_result WHERE seq_folder = '%s' AND seq_id = '%s'", SEQ_FOLDER, SEQ_ID))
    writeNewData      = dbWriteTable(dbCon, name='hla_typing_result', value=RESULTS, row.names=F, append=T)
    dbDisconnect(dbCon) 
#---\    
    message(sprintf("  SeqID = %s : summarized results imported into DB", SEQ_ID))
    write(sprintf("SeqID = %s : HLA-typing Summarized Results Imported into DB, Done.", SEQ_ID), file = LOGFILE, append = TRUE)
#---/

