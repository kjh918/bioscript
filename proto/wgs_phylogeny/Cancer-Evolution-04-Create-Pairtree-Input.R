### packages
    options(stringsAsFactors=FALSE) 
    suppressPackageStartupMessages({
        library("Hmisc")
        library("plyr")
        library("dplyr")
        library("setwidth")
        library("optparse")
        library("GenomicRanges")
        library("data.table")
    })
    indexing = function(X) { sapply(unique(X), function(y) list(y)) }
###
### 
    BaseDir      <- "/data/wgs/RUO_Service"
    BatchID      <- "WGS-26-02"
    AnalysisList <- sprintf("%s/%s/meta/clonality.evolution.analysis.list", BaseDir, BatchID )
    
    VarDataDir <- sprintf("%s/%s/analysis/variant-sites-allele-count", BaseDir, BatchID)
    PycloneDir <- sprintf("%s/%s/analysis/pyclone-vi-variant-clusters", BaseDir, BatchID)
    WorkDir    <- sprintf("%s/%s/analysis/cancer-evolution-pairtree2", BaseDir, BatchID)
    system(sprintf("mkdir -p %s", WorkDir))

    sinfo <- read.table(AnalysisList, header=F, sep=" ") %>% dplyr::rename(SeqID=1, Group=2, Order=3)
    sinfo$SampleID <- make.names(sinfo$SeqID)
    tumor_ids      <- sinfo %>% filter( Group != "normal" ) %>% .$SeqID
    tumor_samples  <- sinfo %>% filter( Group != "normal" ) %>% .$SampleID
###
### 
    wgs_var_count_list <- readRDS(sprintf("%s/%s.wgs.variants.read.count.list.RDS", VarDataDir, BatchID))
    wgs_pyclone_res    <- fread(sprintf("%s/%s.wgs.pyclone-vi.result.tsv", PycloneDir, BatchID))
    variant_list       <- data.frame(
        mutation_id = Reduce(union, lapply(wgs_var_count_list, function(x) x %>% mutate( var_name=paste(x$CHROM,x$POS, x$REF, x$ALT, sep="_")) %>% .$var_name ))
    ) 
    cluster_variant_data <- variant_list %>% left_join(., wgs_pyclone_res[,c("mutation_id","cluster_id")] %>% unique(), by="mutation_id") %>% 
        filter( !is.na(cluster_id)) %>% mutate( var_id = paste0("s", 0:(nrow(.)-1)))
###

### 
    wgs_tso_var_count_list <- readRDS(sprintf("%s/%s.wgs.tso.merged.variants.read.count.list.RDS", VarDataDir, BatchID))
    wgs_tso_pyclone_res    <- fread(sprintf("%s/%s.wgs.tso.merged.pyclone-vi.result.tsv", PycloneDir, BatchID))
    variant_list       <- data.frame(
        mutation_id = Reduce(union, lapply(wgs_tso_var_count_list, function(x) x %>% mutate( var_name=paste(x$CHROM,x$POS, x$REF, x$ALT, sep="_")) %>% .$var_name ))
    ) 
    cluster_variant_data <- variant_list %>% left_join(., wgs_tso_pyclone_res[,c("mutation_id","cluster_id")] %>% unique(), by="mutation_id") %>% 
        filter( !is.na(cluster_id)) %>% mutate( var_id = paste0("s", 0:(nrow(.)-1)))
###

### create pairtree input
    VAR_DATA_LIST = wgs_var_count_list
    # ssm.txt
    pairtree_input_list <- lapply( indexing(tumor_ids), function(tid) {
        # cnv
        purple_dir <- sprintf("%s/%s/%s/hmftools/purple", BaseDir, BatchID, tid)
        cnv_data   <- read.delim(sprintf("%s/%s.purple.cnv.somatic.tsv", purple_dir, tid))
        qc_data    <- read.delim(sprintf("%s/%s.purple.purity.tsv",      purple_dir, tid))
        cnv_data$purity = qc_data$purity       
       
        VARS   <- data.table(VAR_DATA_LIST[[tid]])
        VARS[which(is.na(VARS$alt_reads)), "alt_reads"] <- 0
        VARS   <- VARS %>% mutate( total_reads = ref_reads + alt_reads )
        CNVS   <- data.table(cnv_data)
        PURITY <- unique(CNVS$purity)

        v_gr <- GRanges( VARS$CHROM, IRanges(start=VARS$POS, end=VARS$POS), sample_id=tid, var_key=paste(VARS$CHROM,VARS$POS, VARS$REF, VARS$ALT, sep=":"))
        c_gr <- GRanges( CNVS$chromosome, IRanges(start=CNVS$start, end=CNVS$end), sample_id=tid, copyNumber = CNVS$copyNumber, minorAlleleCopyNumber = CNVS$minorAlleleCopyNumber)

        hits <- findOverlaps(v_gr, c_gr)

        VARS[, copyNumber := 2]
        VARS[, minorAlleleCopyNumber := 1]
        VARS[, purity := PURITY]

        if (length(hits) > 0) {
            hit_dt <- data.table(
                q = queryHits(hits),
                s = subjectHits(hits)
            )
            # if multiple CNV segments overlap a point, keep first
            hit_dt <- hit_dt[!duplicated(q)]

            VARS[hit_dt$q, copyNumber := mcols(c_gr)$copyNumber[hit_dt$s]]
            VARS[hit_dt$q, minorAlleleCopyNumber := mcols(c_gr)$minorAlleleCopyNumber[hit_dt$s]]
        } 

        VARS[, adjusted_total_cn := 2 + purity * (copyNumber - 2)]
        VARS[adjusted_total_cn <= 0 | is.na(adjusted_total_cn), adjusted_total_cn := 2]
        VARS[, var_read_prob := 1 / adjusted_total_cn]
        VARS[var_read_prob > 1, var_read_prob := 1]
        VARS[var_read_prob < 1e-4, var_read_prob := 1e-4]

        RES <- VARS %>% mutate( 
            mutation_id = paste(VARS$CHROM,VARS$POS, VARS$REF, VARS$ALT, sep="_")
        ) %>% left_join( 
            ., cluster_variant_data, by="mutation_id"
        ) %>% filter( 
            !is.na(var_id) 
        ) %>% dplyr::select(
            c("var_id","mutation_id","alt_reads","total_reads","var_read_prob")
        ) %>% dplyr::rename(
            id=1, name=2, var_reads=3, total_reads=4, var_read_prob=5 
        )
        return(RES)
    })
    ssm_data <- pairtree_input_list[[1]]
    for( j in 2:length(pairtree_input_list) ){
        ssm_data <- ssm_data %>% left_join(
            ., pairtree_input_list[[j]][,c("name","var_reads","total_reads","var_read_prob")] %>% 
                dplyr::rename( vr=var_reads, tr=total_reads, vrp=var_read_prob),
            by="name"
        ) %>% mutate( 
            var_reads     = paste(var_reads, vr, sep=","),
            total_reads   = paste(total_reads, tr, sep=","),
            var_read_prob = paste(var_read_prob, vrp, sep=",")
        ) %>% dplyr::select(
            c("id","name","var_reads","total_reads","var_read_prob")
        )
    }

    # param.json
    pairtree_param_var_list <- cluster_variant_data %>% filter( mutation_id %in% ssm_data$name )
    Variant_cluster_list    <- lapply( indexing(pairtree_param_var_list$cluster_id), function(idx) 
        pairtree_param_var_list %>% filter( cluster_id==idx ) %>% .$var_id 
    )
    params <- list(
        samples  = tumor_ids,
        clusters = Variant_cluster_list,
        garbage  = list()
    )

    # save as proper format
    write.table(ssm_data, sprintf('%s/%s.pairtree.ssm.txt', WorkDir, BatchID), quote=F, col.names=T, row.names=F, sep="\t")
    jsonlite::write_json(params, sprintf('%s/%s.pairtree.params.json', WorkDir, BatchID), pretty = TRUE, auto_unbox = TRUE)
####




 
    
