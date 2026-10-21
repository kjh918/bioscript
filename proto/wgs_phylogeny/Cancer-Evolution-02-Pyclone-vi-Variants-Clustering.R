
###
    options(stringsAsFactors=FALSE) 
    suppressPackageStartupMessages({
        library("Hmisc")
        library("plyr")
        library("dplyr")
        library("setwidth")
        library("optparse")
        library("GenomicRanges")
        library("data.table")
        library("BiocParallel")
    })
    indexing = function(X) { sapply(unique(X), function(y) list(y)) }
    SysTime  = function(){format(Sys.time(), "%Y-%m-%d %H:%M:%S")}
    PYCLONE_VI_RUN_CMD <- "singularity exec -B /storage,/data /storage/images/pyclone-vi_v0.1.2.sif pyclone-vi"
    IncludeClass1 <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins",
        "In_Frame_Del","In_Frame_Ins","Splice_Region","Splice_Site"
    )
    IncludeClass2 <- c("Intron","Silent","5'UTR","3'UTR")
    Rounding <- function(x, digits = 0) {   
        posneg <- sign(x)
        z <- abs(x) * 10^digits
        z <- floor(z + 0.5)
        v <- z / 10^digits * posneg
        return(v)
    }
###
### 
    BaseDir      <- "/data/wgs/RUO_Service"
    BatchID      <- "WGS-26-02"
    AnalysisList <- sprintf("%s/%s/meta/clonality.evolution.analysis.list", BaseDir, BatchID )
    GenomeAssembly="hg19"

    VarDataDir <- sprintf("%s/%s/analysis/variant-sites-allele-count", BaseDir, BatchID)   
    WorkDir    <- sprintf("%s/%s/analysis/pyclone-vi-variant-clusters", BaseDir, BatchID)
    system(sprintf("mkdir -p %s", WorkDir))

    sinfo <- read.table(AnalysisList, header=F, sep=" ") %>% dplyr::rename(SeqID=1, Group=2, Order=3)
    sinfo$SampleID <- make.names(sinfo$SeqID)
    tumor_ids      <- sinfo %>% filter( Group != "normal" ) %>% .$SeqID
    tumor_samples  <- sinfo %>% filter( Group != "normal" ) %>% .$SampleID
###
### Variants Clusters : WGS

    wgs_var_count_list <- readRDS(sprintf("%s/%s.wgs.variants.read.count.list.RDS", VarDataDir, BatchID))
    
    VarClassMap <- wgs_var_count_list[[1]][,c("vkey","variant_class")] 
    for( i in 2:length(wgs_var_count_list) ){
        VarClassMap <- VarClassMap %>% left_join(.,
            wgs_var_count_list[[i]][,c("vkey","variant_class")],
            by="vkey"
        )
    } 
    colnames(VarClassMap)[2:ncol(VarClassMap)] <- paste0("v", 1:(ncol(VarClassMap)-1))
    VarClassMap$varclass <- apply(VarClassMap[,-1], 1, function(v) unique(v[!is.na(v)]) )

    wgs_pyclone_input_list <- lapply( indexing(tumor_ids), function(tid) {
   
        wgs_purple_data_dir <- sprintf("%s/%s/%s/hmftools/purple", BaseDir, BatchID, tid)
        wgs_cnv_data <- read.delim(sprintf("%s/%s.purple.cnv.somatic.tsv",wgs_purple_data_dir, tid))
        wgs_qc_data  <- read.delim(sprintf("%s/%s.purple.purity.tsv",wgs_purple_data_dir, tid))
        wgs_cnv_data$purity = wgs_qc_data$purity

        VARS <- data.table(wgs_var_count_list[[tid]])
        VARS$variant_class <- VarClassMap[match(VARS$vkey, VarClassMap$vkey), "varclass"]
        #VARS <- VARS %>% filter( variant_class %in% c(IncludeClass1, IncludeClass2), total_reads >= 10 )
        VARS <- VARS %>% filter( total_reads >= 10 )

        CNVS   <- data.table( wgs_cnv_data )

        VARS$majorAlleleCopyNumber = 1
        VARS$minorAlleleCopyNumber = 1
        VARS$purity = unique(CNVS$purity)
        VARS$normalCopyNumber = 2

        var_gr <- GRanges( VARS$CHROM, IRanges(start=VARS$POS, end=VARS$POS), sample_id=tid, index=paste(VARS$CHROM,VARS$POS,VARS$REF,VARS$ALT, sep="_") )
        cnv_gr <- GRanges( CNVS$chromosome, IRanges(start=CNVS$start, end=CNVS$end), sample_id=tid, 
            majorAlleleCopyNumber = CNVS$majorAlleleCopyNumber, minorAlleleCopyNumber = CNVS$minorAlleleCopyNumber
        )
        hits <- findOverlaps(var_gr, cnv_gr)
        if (length(hits) > 0) {
            hit_dt <- data.table( q = queryHits(hits), s = subjectHits(hits) )
            # if multiple CNV segments overlap a point, keep first
            hit_dt <- hit_dt[!duplicated(q)]
            VARS[hit_dt$q, majorAlleleCopyNumber := mcols(cnv_gr)$majorAlleleCopyNumber[hit_dt$s]]
            VARS[hit_dt$q, minorAlleleCopyNumber := mcols(cnv_gr)$minorAlleleCopyNumber[hit_dt$s]]
        } 

        PycloneInputData <- VARS %>% mutate(
            mutation_id = paste(CHROM,POS,REF,ALT,sep="_"),
            sample_id   = tid,
            majorAlleleCopyNumber = Rounding(majorAlleleCopyNumber, 0),
            minorAlleleCopyNumber = Rounding(minorAlleleCopyNumber, 0)
        ) %>% dplyr::select(
            c("mutation_id","sample_id","ref_reads","alt_reads","majorAlleleCopyNumber","minorAlleleCopyNumber","normalCopyNumber","purity")
        ) %>% dplyr::rename(
            ref_counts=ref_reads, alt_counts=alt_reads, major_cn=majorAlleleCopyNumber, minor_cn=minorAlleleCopyNumber, normal_cn=normalCopyNumber, tumour_content=purity
        )
        return(PycloneInputData)
    })

    all_wgs_pyclone_input <- ldply(wgs_pyclone_input_list) %>% dplyr::select(!c(".id"))

    total_wgs_mutations <- length(unique(all_wgs_pyclone_input$mutation_id))
    if( total_wgs_mutations < 200 ){ WEIGHT_WGS = 5 }else{ WEIGHT_WGS = 2 }

    message(sprintf("\n[ %s ] create pyclone-vi input : %s", SysTime(), BatchID))
    write.table( all_wgs_pyclone_input, sprintf("%s/%s.wgs.pyclone-vi.input.tsv", WorkDir, BatchID), quote=F, col.names=T, row.names=F, sep="\t")
    WGS_PYCLINE_VI_INPUT   <- sprintf("%s/%s.wgs.pyclone-vi.input.tsv",  WorkDir, BatchID)
    WGS_PYCLONE_VI_OUTPUT  <- sprintf("%s/%s.wgs.pyclone-vi.output.h5",  WorkDir, BatchID)
    WGS_PYCLONE_VI_RESULT  <- sprintf("%s/%s.wgs.pyclone-vi.result.tsv", WorkDir, BatchID)
    message(sprintf("[ %s ] run pyclone-vi fit", SysTime()))
    system(sprintf( "%s fit --in-file %s --out-file %s --num-restarts 20 --num-clusters 15 --density beta-binomial --mix-weight-prior %s --seed 123456", 
        PYCLONE_VI_RUN_CMD, WGS_PYCLINE_VI_INPUT, WGS_PYCLONE_VI_OUTPUT, WEIGHT_WGS
    ))
    message(sprintf("[ %s ] convert pyclone-vi result", SysTime()))
    system(sprintf( "%s write-results-file --in-file %s --out-file %s", PYCLONE_VI_RUN_CMD, WGS_PYCLONE_VI_OUTPUT, WGS_PYCLONE_VI_RESULT ))

###
### Variants Clusters : TSO
    tso_var_count_list <- readRDS(sprintf("%s/%s.tso.variants.read.count.list.RDS", VarDataDir, BatchID))

    VarClassMapTso <- tso_var_count_list[[1]][,c("vkey","variant_class")] 
    for( i in 2:length(tso_var_count_list) ){
        VarClassMapTso <- VarClassMapTso %>% left_join(.,
            tso_var_count_list[[i]][,c("vkey","variant_class")],
            by="vkey"
        )
    } 
    colnames(VarClassMapTso)[2:ncol(VarClassMapTso)] <- paste0("v", 1:(ncol(VarClassMapTso)-1))
    VarClassMapTso$varclass <- apply(VarClassMapTso[,-1], 1, function(v) unique(v[!is.na(v)]) )

    tso_pyclone_input_list <- lapply( indexing(tumor_ids), function(tid) {
   
        tso_purple_data_dir <- sprintf("%s/%s/%s/tso500/tso_hmftools/08.purple", BaseDir, BatchID, tid)
        tso_cnv_data <- read.delim(sprintf("%s/%s.purple.cnv.somatic.tsv",tso_purple_data_dir, tid))
        tso_qc_data  <- read.delim(sprintf("%s/%s.purple.purity.tsv",tso_purple_data_dir, tid))
        tso_cnv_data$purity = tso_qc_data$purity

        VARS <- data.table(tso_var_count_list[[tid]])
        VARS$variant_class <- VarClassMapTso[match(VARS$vkey, VarClassMapTso$vkey), "varclass"]
        #VARS <- VARS %>% filter( variant_class %in% c(IncludeClass1, IncludeClass2), vaf >= 0.03 )
        VARS <- VARS %>% filter( vaf >= 0.03 )

        CNVS   <- data.table( tso_cnv_data )

        VARS$majorAlleleCopyNumber = 1
        VARS$minorAlleleCopyNumber = 1
        VARS$purity = unique(CNVS$purity)
        VARS$normalCopyNumber = 2

        var_gr <- GRanges( VARS$CHROM, IRanges(start=VARS$POS, end=VARS$POS), sample_id=tid, index=paste(VARS$CHROM,VARS$POS,VARS$REF,VARS$ALT, sep="_") )
        cnv_gr <- GRanges( CNVS$chromosome, IRanges(start=CNVS$start, end=CNVS$end), sample_id=tid, 
            majorAlleleCopyNumber = CNVS$majorAlleleCopyNumber, minorAlleleCopyNumber = CNVS$minorAlleleCopyNumber
        )
        hits <- findOverlaps(var_gr, cnv_gr)
        if (length(hits) > 0) {
            hit_dt <- data.table( q = queryHits(hits), s = subjectHits(hits) )
            # if multiple CNV segments overlap a point, keep first
            hit_dt <- hit_dt[!duplicated(q)]
            VARS[hit_dt$q, majorAlleleCopyNumber := mcols(cnv_gr)$majorAlleleCopyNumber[hit_dt$s]]
            VARS[hit_dt$q, minorAlleleCopyNumber := mcols(cnv_gr)$minorAlleleCopyNumber[hit_dt$s]]
        } 

        PycloneInputData <- VARS %>% mutate(
            mutation_id = paste(CHROM,POS,REF,ALT,sep="_"),
            sample_id   = tid,
            majorAlleleCopyNumber = floor(majorAlleleCopyNumber+0.5),
            minorAlleleCopyNumber = floor(minorAlleleCopyNumber+0.5)
        ) %>% dplyr::select(
            c("mutation_id","sample_id","ref_reads","alt_reads","majorAlleleCopyNumber","minorAlleleCopyNumber","normalCopyNumber","purity")
        ) %>% dplyr::rename(
            ref_counts=ref_reads, alt_counts=alt_reads, major_cn=majorAlleleCopyNumber, minor_cn=minorAlleleCopyNumber, normal_cn=normalCopyNumber, tumour_content=purity
        )
        return(PycloneInputData)
    })
    all_tso_pyclone_input <- ldply(tso_pyclone_input_list) %>% dplyr::select(!c(".id"))

    total_tso_mutations <- length(unique(all_tso_pyclone_input$mutation_id))
    if( total_tso_mutations < 200 ){ WEIGHT_TSO = 5 }else{ WEIGHT_TSO = 2 }

    message(sprintf("\n[ %s ] create pyclone-vi input : %s", SysTime(), BatchID))
    write.table( all_tso_pyclone_input, sprintf("%s/%s.tso.pyclone-vi.input.tsv", WorkDir, BatchID), quote=F, col.names=T, row.names=F, sep="\t")
    TSO_PYCLINE_VI_INPUT   <- sprintf("%s/%s.tso.pyclone-vi.input.tsv", WorkDir, BatchID)
    TSO_PYCLONE_VI_OUTPUT  <- sprintf("%s/%s.tso.pyclone-vi.output.h5", WorkDir, BatchID)
    TSO_PYCLONE_VI_RESULT  <- sprintf("%s/%s.tso.pyclone-vi.result.tsv", WorkDir, BatchID)
    message(sprintf("[ %s ] run pyclone-vi fit", SysTime()))
    system(sprintf( "%s fit --in-file %s --out-file %s --num-clusters 15 --density beta-binomial --num-restarts 20 --mix-weight-prior %s --seed 123456", 
        PYCLONE_VI_RUN_CMD, TSO_PYCLINE_VI_INPUT, TSO_PYCLONE_VI_OUTPUT, WEIGHT_TSO
    ))
    message(sprintf("[ %s ] convert pyclone-vi result", SysTime()))
    system(sprintf( "%s write-results-file --in-file %s --out-file %s", PYCLONE_VI_RUN_CMD, TSO_PYCLONE_VI_OUTPUT, TSO_PYCLONE_VI_RESULT ))
###
### Variants Clusters : WGS-TSO-MERGED
    
    wgs_tso_merged_var_count_list <- readRDS(sprintf("%s/%s.wgs.tso.merged.variants.read.count.list.RDS", VarDataDir, BatchID))

    VarClassMapWgsTso <- wgs_tso_merged_var_count_list[[1]][,c("vkey","variant_class")] 
    for( i in 2:length(wgs_tso_merged_var_count_list) ){
        VarClassMapWgsTso <- VarClassMapWgsTso %>% left_join(.,
            wgs_tso_merged_var_count_list[[i]][,c("vkey","variant_class")],
            by="vkey"
        )
    } 
    colnames(VarClassMapWgsTso)[2:ncol(VarClassMapWgsTso)] <- paste0("v", 1:(ncol(VarClassMapWgsTso)-1))
    VarClassMapWgsTso$varclass <- apply(VarClassMapWgsTso[,-1], 1, function(v) unique(v[!is.na(v)]) )

    wgs_tso_pyclone_input_list <- lapply( indexing(tumor_ids), function(tid) {
   
        wgs_purple_data_dir <- sprintf("%s/%s/%s/hmftools/purple", BaseDir, BatchID, tid)
        wgs_cnv_data <- read.delim(sprintf("%s/%s.purple.cnv.somatic.tsv",wgs_purple_data_dir, tid))
        wgs_qc_data  <- read.delim(sprintf("%s/%s.purple.purity.tsv",wgs_purple_data_dir, tid))
        wgs_cnv_data$purity = wgs_qc_data$purity

        VARS <- data.table(wgs_tso_merged_var_count_list[[tid]])
        VARS$variant_class <- VarClassMapWgsTso[match(VARS$vkey, VarClassMapWgsTso$vkey), "varclass"]
        # VARS <- VARS %>% filter( variant_class %in% c(IncludeClass1, IncludeClass2), total_reads >= 10 )
        VARS <- VARS %>% filter( total_reads >= 10 )

        CNVS   <- data.table( wgs_cnv_data )

        VARS$majorAlleleCopyNumber = 1
        VARS$minorAlleleCopyNumber = 1
        VARS$purity = unique(CNVS$purity)
        VARS$normalCopyNumber = 2

        var_gr <- GRanges( VARS$CHROM, IRanges(start=VARS$POS, end=VARS$POS), sample_id=tid, index=paste(VARS$CHROM,VARS$POS,VARS$REF,VARS$ALT, sep="_") )
        cnv_gr <- GRanges( CNVS$chromosome, IRanges(start=CNVS$start, end=CNVS$end), sample_id=tid, 
            majorAlleleCopyNumber = CNVS$majorAlleleCopyNumber, minorAlleleCopyNumber = CNVS$minorAlleleCopyNumber
        )
        hits <- findOverlaps(var_gr, cnv_gr)
        if (length(hits) > 0) {
            hit_dt <- data.table( q = queryHits(hits), s = subjectHits(hits) )
            # if multiple CNV segments overlap a point, keep first
            hit_dt <- hit_dt[!duplicated(q)]
            VARS[hit_dt$q, majorAlleleCopyNumber := mcols(cnv_gr)$majorAlleleCopyNumber[hit_dt$s]]
            VARS[hit_dt$q, minorAlleleCopyNumber := mcols(cnv_gr)$minorAlleleCopyNumber[hit_dt$s]]
        } 

        PycloneInputData <- VARS %>% mutate(
            mutation_id = paste(CHROM,POS,REF,ALT,sep="_"),
            sample_id   = tid,
            majorAlleleCopyNumber = Rounding(majorAlleleCopyNumber, 0),
            minorAlleleCopyNumber = Rounding(minorAlleleCopyNumber, 0)
        ) %>% dplyr::select(
            c("mutation_id","sample_id","ref_reads","alt_reads","majorAlleleCopyNumber","minorAlleleCopyNumber","normalCopyNumber","purity")
        ) %>% dplyr::rename(
            ref_counts=ref_reads, alt_counts=alt_reads, major_cn=majorAlleleCopyNumber, minor_cn=minorAlleleCopyNumber, normal_cn=normalCopyNumber, tumour_content=purity
        )
        return(PycloneInputData)
    })
    all_merged_pyclone_input <- ldply(wgs_tso_pyclone_input_list) %>% dplyr::select(!c(".id"))

    total_tso_mutations <- length(unique(all_merged_pyclone_input$mutation_id))
    if( total_tso_mutations < 200 ){ WEIGHT_MERGED = 5 }else{ WEIGHT_MERGED = 2 }

    message(sprintf("\n[ %s ] create pyclone-vi input : %s", SysTime(), BatchID))
    write.table( all_merged_pyclone_input, sprintf("%s/%s.wgs.tso.merged.pyclone-vi.input.tsv", WorkDir, BatchID), quote=F, col.names=T, row.names=F, sep="\t")
    MERGED_PYCLINE_VI_INPUT   <- sprintf("%s/%s.wgs.tso.merged.pyclone-vi.input.tsv",  WorkDir, BatchID)
    MERGED_PYCLONE_VI_OUTPUT  <- sprintf("%s/%s.wgs.tso.merged.pyclone-vi.output.h5",  WorkDir, BatchID)
    MERGED_PYCLONE_VI_RESULT  <- sprintf("%s/%s.wgs.tso.merged.pyclone-vi.result.tsv", WorkDir, BatchID)
    message(sprintf("[ %s ] run pyclone-vi fit", SysTime()))
    system(sprintf( "%s fit --in-file %s --out-file %s --num-clusters 15 --density beta-binomial --num-restarts 20 --mix-weight-prior %s --seed 123456", 
        PYCLONE_VI_RUN_CMD, MERGED_PYCLINE_VI_INPUT, MERGED_PYCLONE_VI_OUTPUT, WEIGHT_MERGED
    ))
    message(sprintf("[ %s ] convert pyclone-vi result", SysTime()))
    system(sprintf( "%s write-results-file --in-file %s --out-file %s", PYCLONE_VI_RUN_CMD, MERGED_PYCLONE_VI_OUTPUT, MERGED_PYCLONE_VI_RESULT ))

###




