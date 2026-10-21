

options(stringsAsFactors=FALSE)   
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Hmisc"))
suppressPackageStartupMessages(library("RMySQL"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("openxlsx"))

#---| PARSE ARGS |
    option_list = list( make_option(c("--seqFolder"), action="store", default=NA, type="character", help="SEQ_FOLDER")    ),
    option_list = list( make_option(c("--seqID"),     action="store", default=NA, type="character", help="SEQ_ID")        )   
  ##  
    ARGS = parse_args(OptionParser(option_list=option_list))
    SEQ_FOLDER = ARGS$seqFolder
    SEQ_ID     = ARGS$seqID
  ##
    DATA_DIR = sprintf("/data/wes/%s",SEQ_FOLDER)
#---/
#---| HLA REPORT FILE GENERATION (Report)
    dbCon = dbConnect(dbDriver("MySQL"), host='192.168.0.35', user='gcx', port=3306, password='gencurix!!', db='gcx_wes' )
    HLA = dbGetQuery( dbCon, sprintf("SELECT * FROM hla_typing_result WHERE seq_folder = '%s' ", SEQ_FOLDER ))
    reportHLA = HLA %>% filter( report == 1 ) 
  ##
    ReportType1 = ldply(lapply( unique(reportHLA$seq_id), function(sq)
    {
        res = rbind(
            reportHLA %>% filter( seq_id == sq) %>% select(c('locus','allele1_optitype')) %>% mutate(allele=1) %>% rename(HLA_type=allele1_optitype),
            reportHLA %>% filter( seq_id == sq) %>% select(c('locus','allele2_optitype')) %>% mutate(allele=2) %>% rename(HLA_type=allele2_optitype)
        ) %>% arrange( locus, allele )
        res2 = reshape( res, idvar='allele', timevar='locus', direction='wide') %>% mutate( seq_id = sq ) %>% 
            rename( 'HLA-A'=HLA_type.A, 'HLA-B'=HLA_type.B, 'HLA-C'=HLA_type.C) %>% relocate( 5,1,2,3,4 )
        return(res2)
    }))
  ##
    ReportType2 = ldply(lapply(unique(reportHLA$seq_id), function(sq)
    {
        cbind(
            reportHLA %>% filter( seq_id == sq, locus == 'A' ) %>% select(c('seq_id','allele1_optitype','allele2_optitype')) %>% 
                rename( 'HLA-A allele1'=allele1_optitype, 'HLA-A allele2'=allele2_optitype ),
            reportHLA %>% filter( seq_id == sq, locus == 'B' ) %>% select(c('allele1_optitype','allele2_optitype')) %>% 
                rename( 'HLA-B allele1'=allele1_optitype, 'HLA-B allele2'=allele2_optitype ),
            reportHLA %>% filter( seq_id == sq, locus == 'C' ) %>% select(c('allele1_optitype','allele2_optitype')) %>% 
                rename( 'HLA-C allele1'=allele1_optitype, 'HLA-C allele2'=allele2_optitype )
        ) 
    }))
  ##
    ReportType3 = ldply(lapply( unique(reportHLA$seq_id), function(sq)
    {
        res = rbind(
            reportHLA %>% filter( seq_id == sq) %>% select(c('locus','allele1_optitype')) %>% mutate(allele=1) %>% rename(HLA_type=allele1_optitype),
            reportHLA %>% filter( seq_id == sq) %>% select(c('locus','allele2_optitype')) %>% mutate(allele=2) %>% rename(HLA_type=allele2_optitype)
        ) %>% arrange( locus, allele )
        res2 = res %>% group_by(locus) %>% summarise( "HLA-Type" = paste(HLA_type, collapse="\n") ) %>% mutate( 'Seq-ID' = sq) %>% 
            mutate(locus=gsub("^","HLA-", locus)) %>% rename('HLA-Class'= locus) %>% relocate(3,1,2)
        return(res2)
    }))
  ##
    supplHLA = HLA %>% filter( report == 0 ) 
    SupplReport = ldply(lapply( unique(supplHLA$seq_id), function(sq) 
    {
        res = rbind(
            supplHLA %>% filter( seq_id == sq) %>% select(c('locus','allele1_hla_full')) %>% mutate(allele=1) %>% rename(HLA_type=allele1_hla_full),
            supplHLA %>% filter( seq_id == sq) %>% select(c('locus','allele2_hla_full')) %>% mutate(allele=2) %>% rename(HLA_type=allele2_hla_full)
        ) %>% arrange( locus, allele )
        res2 = reshape( res, idvar='allele', timevar='locus', direction='wide') %>% 
            rename( 'HLA-E'=HLA_type.E, 'HLA-F'=HLA_type.F, 'HLA-G'=HLA_type.G ) 
        colnames(res2) = gsub("^HLA_type.", "", colnames(res2))
        res2 = rbind(res2, apply(res2, 2, function(y) paste(y, collapse="\n"))) %>% select(2:11) %>% mutate( 'Seq-ID' = sq ) %>% relocate( 11,1:10 ) %>%
            slice(3)
        return(res2)    
    }))
  #-| Save as XLSX-file
    HLA_REPORT = openxlsx::createWorkbook()
    openxlsx::addWorksheet(HLA_REPORT, "MHC-I")
    openxlsx::addWorksheet(HLA_REPORT, "MHC-II")
    openxlsx::writeDataTable(HLA_REPORT,"MHC-I", ReportType3)
    openxlsx::writeDataTable(HLA_REPORT,"MHC-II",SupplReport)
    hla_style = openxlsx::createStyle( fontSize = 11, halign = "center", valign = "center", wrapText = TRUE )
    openxlsx::addStyle(HLA_REPORT, sheet = 1, hla_style, rows = 1:100, cols = 1:3, gridExpand = TRUE)
    openxlsx::addStyle(HLA_REPORT, sheet = 2, hla_style, rows = 1:100, cols = 1:20, gridExpand = TRUE)
    openxlsx::setColWidths(HLA_REPORT, 1, cols = 1:3,  widths = 15)
    openxlsx::setColWidths(HLA_REPORT, 2, cols = 2:11,  widths = 30)
    openxlsx::removeFilter(HLA_REPORT, 1) ; openxlsx::removeFilter(HLA_REPORT, 2)
    openxlsx::saveWorkbook(HLA_REPORT, sprintf("%s/meta/%s_HLA_typing_results.xlsx", DATA_DIR, SEQ_FOLDER ), overwrite = TRUE)
    #---\     
        message("--------------------------------------------------------------------------------------")
        write("--------------------------------------------------------------------------------------", file=LogFileName, append=T)
        message(sprintf("  HLA-typing Result Saved as XLSX-File."))
        write(sprintf("  HLA-typing Result Saved as XLSX-File."), file=LogFileName, append=T)
        message(sprintf("\n  FileName = %s/meta/%s_HLA_typing_results.xlsx", DATA_DIR, SEQ_FOLDER))
        write(sprintf("\n  FileName = %s/meta/%s_HLA_typing_results.xlsx", DATA_DIR, SEQ_FOLDER), file=LogFileName, append=T)
        message("\n------------------------------------------------------------------------------------")
        write("\n------------------------------------------------------------------------------------", file=LogFileName, append=T)
        message(sprintf("[ %s ] HLA typing result summarization finished.", date()))
        write(sprintf("[ %s ] HLA typing result summarization finished.", date()), file=LogFileName, append=T)
       # close(LOG)
    #---/


