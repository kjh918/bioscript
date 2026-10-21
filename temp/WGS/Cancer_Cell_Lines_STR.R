

    options(stringsAsFactors=FALSE) 
    suppressPackageStartupMessages(library("Hmisc"))
    suppressPackageStartupMessages(library("plyr"))
    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("setwidth"))
    suppressPackageStartupMessages(library("RMySQL"))
    suppressPackageStartupMessages(library("Cairo"))
    grDevices::X11.options(type='cairo')
    options(device='x11')
    options(bitmapType='cairo') 

    # cell lines STR

    cclstr = openxlsx::read.xlsx("/data/cbNIPT/meta/Cell.Lines.STR.Values.xlsx")
    str_runlist = data.table::fread("/data/cbNIPT/meta/str.run.list", header=F) %>% as.data.frame()
    CodisSTRmarkers = read.delim("/storage/home/kangsm/myDB/STR_references/hg38.annotated.markers.raw.bed") %>% filter( !is.na(note))

    ## HTR8
         HTR8_STR = CodisSTRmarkers[,c(1,2,3,7,4,5,6)] %>% 
            left_join(., cclstr %>% filter( Cell_Name == "HTR8" ) %>% dplyr::rename(str_id=STR_ID, Cell=Cell_Name, STR_Db=STR_Values) , by="str_id") %>% 
            left_join(., read.delim("/data/cbNIPT/meta/HTR8_gcx_str_test_res.tsv"), by="str_id")

        for( j in c(36:46))
        {
            SeqID      = str_runlist[j, 1]
            SampleID   = str_runlist[j, 2]
            str_res = read.delim(sprintf("/data/cbNIPT/processed_dataset/%s/str/%s.STR.called.list.tsv", SeqID, SampleID)) %>% dplyr::select(c("str_id","allele_motif_repeats", "dp"))
            str_res$str = apply(str_res, 1, function(w) paste0(paste0(unique(unlist(strsplit(w[2], "\\|"))), collapse=","), " (", as.numeric(w[3]), ")") )

            HTR8_STR = HTR8_STR %>% left_join(., str_res[,c("str_id","str")], by="str_id")
            colnames(HTR8_STR)[ncol(HTR8_STR)] = SampleID
        }
        openxlsx::write.xlsx(HTR8_STR, "HTR8_STR_results.xlsx")
    #
    ## JEG3
        JEG3_STR = CodisSTRmarkers[,c(1,2,3,7,4,5,6)] %>% 
            left_join(., cclstr %>% filter( Cell_Name == "JEG3" ) %>% dplyr::rename(str_id=STR_ID, Cell=Cell_Name, STR_Db=STR_Values) , by="str_id") %>% 
            left_join(., read.delim("/data/cbNIPT/meta/JEG3_gcx_str_test_res.tsv"), by="str_id")
        for( j in c(24,32,33,34,35,47))
        {
            SeqID      = str_runlist[j, 1]
            SampleID   = str_runlist[j, 2]
            str_res = read.delim(sprintf("/data/cbNIPT/processed_dataset/%s/str/%s.STR.called.list.tsv", SeqID, SampleID)) %>% dplyr::select(c("str_id","allele_motif_repeats", "dp"))
            str_res$str = apply(str_res, 1, function(w) paste0(paste0(unique(unlist(strsplit(w[2], "\\|"))), collapse=","), " (", as.numeric(w[3]), ")") )

            JEG3_STR = JEG3_STR %>% left_join(., str_res[,c("str_id","str")], by="str_id")
            colnames(JEG3_STR)[ncol(JEG3_STR)] = SampleID
        }
        openxlsx::write.xlsx(JEG3_STR, "JEG3_STR_results.xlsx")
    #
    ## Cell Lines (A549)
        j=64
        SeqID      = str_runlist[j, 1]
        SampleID   = str_runlist[j, 2]
        SampleName = str_runlist[j, 3]

        A549_STR = CodisSTRmarkers[,c(1,2,3,7,4,5,6)] %>% 
            left_join(., cclstr %>% filter( Cell_Name == SampleName ) %>% dplyr::rename(str_id=STR_ID, Cell=Cell_Name, STR_Db=STR_Values) , by="str_id") 
        
        for( j in c(64,65) )
        {
            SeqID      = str_runlist[j, 1]
            SampleID   = str_runlist[j, 2]
            SampleName = str_runlist[j, 3]

            str_res = read.delim(sprintf("/data/cbNIPT/processed_dataset/%s/str/%s.STR.called.list.tsv", SeqID, SampleID)) %>% dplyr::select(c("str_id","allele_motif_repeats", "dp"))
            str_res$str = apply(str_res, 1, function(w) paste0(paste0(unique(unlist(strsplit(w[2], "\\|"))), collapse=","), " (", as.numeric(w[3]), ")") )
            A549_STR = A549_STR %>% left_join(., str_res[,c("str_id","str")], by="str_id")
            colnames(A549_STR)[ncol(A549_STR)] = SampleID
        }
        openxlsx::write.xlsx(A549_STR, "A549_STR_results.xlsx")

    #
    ## Cell Lines (Other)
    CCL_STR = lapply( c(55,56,58,59,60,61,62,63), function(j)
    {
        SeqID      = str_runlist[j, 1]
        SampleID   = str_runlist[j, 2]
        SampleName = str_runlist[j, 3]
        
        CCL_STR = CodisSTRmarkers[,c(1,2,3,7,4,5,6)] %>% 
            left_join(., cclstr %>% filter( Cell_Name == SampleName ) %>% dplyr::rename(str_id=STR_ID, Cell=Cell_Name, STR_Db=STR_Values) , by="str_id")
        str_res = read.delim(sprintf("/data/cbNIPT/processed_dataset/%s/str/%s.STR.called.list.tsv", SeqID, SampleID)) %>% dplyr::select(c("str_id","allele_motif_repeats", "dp"))
            str_res$str = apply(str_res, 1, function(w) paste0(paste0(unique(unlist(strsplit(w[2], "\\|"))), collapse=","), " (", as.numeric(w[3]), ")") )
            CCL_STR = CCL_STR %>% left_join(., str_res[,c("str_id","str")], by="str_id")
            colnames(CCL_STR)[ncol(CCL_STR)] = SampleID
        return(CCL_STR)
    })
        openxlsx::write.xlsx(CCL_STR, "CellLines_STR_results.xlsx")

    #
    ##
        RM_STR = CodisSTRmarkers[,c(1,2,3,7,4,5,6)] 
        
        for( j in c(66,67,68,69) )
        {
            SeqID      = str_runlist[j, 1]
            SampleID   = str_runlist[j, 2]
            SampleName = str_runlist[j, 3]

            str_res = read.delim(sprintf("/data/cbNIPT/processed_dataset/%s/str/%s.STR.called.list.tsv", SeqID, SampleID)) %>% dplyr::select(c("str_id","allele_motif_repeats", "dp"))
            str_res$str = apply(str_res, 1, function(w) paste0(paste0(unique(unlist(strsplit(w[2], "\\|"))), collapse=","), " (", as.numeric(w[3]), ")") )
            RM_STR = RM_STR %>% left_join(., str_res[,c("str_id","str")], by="str_id")
            colnames(RM_STR)[ncol(RM_STR)] = SampleID
        }
        openxlsx::write.xlsx(RM_STR, "Coriell_STR_results.xlsx")

    ## 

        WGA_STR = CodisSTRmarkers[,c(1,2,3,7,4,5,6)] 
        
        for( j in c(6,15,7,16,8,17,9,18,10,19,11,20,12,21,13,22,14,23) )
        {
            SeqID      = str_runlist[j, 1]
            SampleID   = str_runlist[j, 2]
            SampleName = str_runlist[j, 3]

            str_res = read.delim(sprintf("/data/cbNIPT/processed_dataset/%s/str/%s.STR.called.list.tsv", SeqID, SampleID)) %>% dplyr::select(c("str_id","allele_motif_repeats", "dp"))
            str_res$str = apply(str_res, 1, function(w) paste0(paste0(unique(unlist(strsplit(w[2], "\\|"))), collapse=","), " (", as.numeric(w[3]), ")") )
            WGA_STR = WGA_STR %>% left_join(., str_res[,c("str_id","str")], by="str_id")
            colnames(WGA_STR)[ncol(WGA_STR)] = SampleID
        }
        openxlsx::write.xlsx(WGA_STR, "WGA_WGS_STR_results.xlsx")


