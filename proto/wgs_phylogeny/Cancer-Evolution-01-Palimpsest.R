
### packages & resources
    options(stringsAsFactors=FALSE) 
    suppressPackageStartupMessages({
        library("Hmisc")
        library("plyr")
        library("dplyr")
        library("setwidth")
        library("optparse")
        library("RMySQL")
        library("vcfR")
    })
    suppressWarnings(suppressPackageStartupMessages(library("Palimpsest")))
    indexing = function(X) { sapply(unique(X), function(y) list(y)) }

    No.NA.Merge <- function(x) {   
        x = x[!is.na(x)]
        x = x[ x != "NA" ]
        x = x[ x != "" ]
        x2 = unique(x)
        if( length(x2) == 0 ){ v = "" }else{ v = paste0(x2, collapse=";") }
        return(v)
    }

###
### params
    BaseDir        <- "/data/wgs/RUO_Service"
    BatchID        <- "WGS-26-02"
    TumorSeqIdList <- "WGS-26-02_03,WGS-26-02_02"
    GenomeAssembly <- "hg19"
    CancerType     <- "BRCA"

    VarDataDir <- sprintf("%s/%s/analysis/variant-sites-allele-count",  BaseDir, BatchID)
    PycloneDir <- sprintf("%s/%s/analysis/pyclone-vi-variant-clusters", BaseDir, BatchID)
    WorkDir    <- sprintf("%s/%s/analysis/clonality-palimpsest", BaseDir, BatchID)
    if( !dir.exists(WorkDir) ){ system(sprintf("mkdir -p %s", WorkDir)) }

    AnalysisList <- sprintf("%s/%s/meta/clonality.evolution.analysis.list", BaseDir, BatchID)

    sinfo <- read.table(AnalysisList, header=F, sep=" ") %>% dplyr::rename(SeqID=1, Group=2, Order=3, SampleName=4)
    sinfo$SampleID <- make.names(sinfo$SeqID)
    tumor_ids      <- sinfo %>% filter( Group != "normal" ) %>% .$SeqID
    tumor_samples  <- sinfo %>% filter( Group != "normal" ) %>% .$SampleID
###
###
    dbCon = dbConnect(dbDriver("MySQL"), host='192.168.0.34', user='gcx', port=3306, password='gencurix!!', db='public_database' )
    DirverGenes <- dbGetQuery(dbCon, "SELECT * FROM intogen_cancer_driver_catalog")
    dbDisconnect(dbCon)
    DirverGenes <- DirverGenes %>% filter( CANCER_TYPE == CancerType ) %>% group_by(SYMBOL) %>% 
        mutate( ACT_MODE=paste0(unique(ROLE), collapse=",") ) %>% dplyr::select(
        c("SYMBOL","CANCER_TYPE","ACT_MODE","IS_DRIVER")
    ) %>% unique() %>% as.data.frame()

    GENOME_FASTA="/storage/references_and_index/hg19/hmftools/Homo_sapiens.GRCh37.GATK.illumina.fasta"

    if( GenomeAssembly == "hg19" ){
        suppressPackageStartupMessages({library("BSgenome.Hsapiens.UCSC.hg19")})
        GENOME=BSgenome.Hsapiens.UCSC.hg19
        GRCh37_TYPE = TRUE
    }else{
        suppressPackageStartupMessages({library("BSgenome.Hsapiens.UCSC.hg38")})
        GENOME=BSgenome.Hsapiens.UCSC.hg38
        GRCh37_TYPE = FALSE
    }
    IncludeClass <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins",
        "In_Frame_Del","In_Frame_Ins","Splice_Region","Splice_Site","Intron","Silent","5'UTR","3'UTR"
    )

###
### run
       
    # input-variant-data
    wgs_var_count_list <- readRDS(sprintf("%s/%s.wgs.variants.read.count.list.RDS", VarDataDir, BatchID))

    VarInfo <- ldply(wgs_var_count_list) %>% group_by(vkey) %>% reframe(
        Gene          = No.NA.Merge(gene),
        VariantClass  = No.NA.Merge(variant_class),
        VariantType   = No.NA.Merge(variant_type),
        HGVSp         = No.NA.Merge(hgvsp),
        Tier          = No.NA.Merge(TIER)
    )

    var_data_list <- lapply( indexing(tumor_ids), function(tid) {
        VARS <- wgs_var_count_list[[tid]] %>% left_join(., VarInfo, by="vkey") %>% 
            filter( alt_reads > 0, VariantClass %in% IncludeClass ) %>% 
            mutate( Sample = tid ) %>% dplyr::select(
                c("Sample","VariantType","CHROM","POS","REF","ALT","alt_reads","total_reads","matched_normal_dp","Gene","HGVSp")
            ) %>% dplyr::rename(
                Sample=1, Type=2, CHROM=3, POS=4, REF=5, ALT=6, Tumor_Varcount=7, Tumor_Depth=8, Normal_Depth=9, Gene=10, HGVSp=11
            )
        VARS$Type <- sapply( VARS$Type, function(vt) {
            if( vt %nin% c("INS","DEL") ){ vtype = "SNV" }else{ vtype = vt }
            return(vtype)
        })
        return(VARS)
    })
        
    # input-cnv-data
    cnv_data_list <- lapply( indexing(tumor_ids), function(tid) {
        purple_dir <- sprintf("%s/%s/%s/hmftools/purple", BaseDir, BatchID, tid)
        seg_read   <- read.delim(sprintf("%s/%s.purple.segment.tsv",purple_dir,tid)) %>% mutate(seg_id=1:nrow(.))
        cnv_read   <- read.delim(sprintf("%s/%s.purple.cnv.somatic.tsv",purple_dir,tid)) %>% mutate(cnv_id=1:nrow(.))
        qc_read    <- read.delim(sprintf("%s/%s.purple.purity.tsv",purple_dir,tid))
        #
        gr_seg <- GRanges(seg_read$chromosome, IRanges(start=seg_read$start, end=seg_read$end), ratio=seg_read$observedTumorRatio)
        gr_cnv <- GRanges(cnv_read$chromosome, IRanges(start=cnv_read$start, end=cnv_read$end), cnv_id=cnv_read$cnv_id)
        #
        hits  <- findOverlaps(gr_cnv, gr_seg)
        ov    <- pintersect( gr_cnv[queryHits(hits)], gr_seg[subjectHits(hits)] )
        ov_dt <- data.table::data.table(
            cnv_id = mcols(gr_cnv)$cnv_id[queryHits(hits)],
            ratio  = mcols(gr_seg)$ratio[subjectHits(hits)],
            overlap_bp = width(ov)
        )
        ratio_by_cnv <- ov_dt %>% group_by(cnv_id) %>% reframe(
            tumor_ratio=weighted.mean(ratio, overlap_bp, na.rm = TRUE),
            ov_length=sum(overlap_bp)
        )
        cnv_seg <- merge(cnv_read, ratio_by_cnv, by = "cnv_id", all.x = TRUE, sort = FALSE)
        cnv_res <- cnv_seg %>% mutate( sample=tid, ploidy=qc_read$ploidy ) %>% 
            dplyr::select(c( "sample","chromosome","start","end","tumor_ratio","copyNumber","minorAlleleCopyNumber","majorAlleleCopyNumber","ploidy")) %>% 
            dplyr::rename(Sample=1, CHROM=2, POS_START=3, POS_END=4, LogR=5, ntot=6, Nmin=7, Nmaj=8, Ploidy=9) %>% 
            filter( LogR != 0 ) %>% mutate( LogR = log2(LogR) )
        return(cnv_res)      
    })      

    # annot-data
    annot_data_list <- lapply( indexing(tumor_ids), function(tid) {
        purple_dir <- sprintf("%s/%s/%s/hmftools/purple", BaseDir, BatchID, tid)
        qc_read    <- read.delim(sprintf("%s/%s.purple.purity.tsv",purple_dir,tid))
        annot_res  <- data.frame(
            Sample = tid,
            Gender = ifelse(qc_read$gender == "FEMALE", "F","M"),
            Purity = qc_read$purity
        )
        return(annot_res)
    })

### 
    VCF_DATA   <- var_data_list   %>% ldply() %>% dplyr::select(!c(".id"))
    CNV_DATA   <- cnv_data_list   %>% ldply() %>% dplyr::select(!c(".id"))
    ANNOT_DATA <- annot_data_list %>% ldply() %>% dplyr::select(!c(".id"))
 
    VCF_DATA$Type <- gsub("SNP", "SNV", VCF_DATA$Type)
    VCF_DATA$Type <- gsub("MNP", "SNV", VCF_DATA$Type)
    VCF_DATA$Type <- gsub("TNP", "SNV", VCF_DATA$Type)
    VCF_DATA <- VCF_DATA %>% filter( Type %in% c("SNV","INS","DEL") )

    ANNOT_VCF_DATA <- suppressWarnings(annotate_VCF(
        vcf          = VCF_DATA, 
        add_DBS_cats = FALSE, 
        add_ID_cats  = TRUE, 
        ref_fasta    = GENOME_FASTA,
        ref_genome   = GENOME,
        GRCh37_fasta = GRCh37_TYPE
    ))
###
### COSMIC SBS
    # Breast cancer-focused SBS set
    #selected_sigs <- c("SBS1", "SBS2", "SBS3", "SBS5", "SBS13", "SBS18")
    #SBS_ref_sigs <- SBS_cosmic[rownames(SBS_cosmic) %in% selected_sigs, ]
    
    SBS_SigExp <- deconvolution_fit(
        input_matrices    = palimpsest_input(vcf = ANNOT_VCF_DATA, Type = "SBS"), 
        input_signatures  = SBS_cosmic,
        resdir            = WorkDir, 
        signature_colours = signature_colour_generator(rownames(SBS_cosmic)), 
        input_vcf         = ANNOT_VCF_DATA,
        doplot            = FALSE
    )
    ID_SigExp <- deconvolution_fit(
        input_matrices    = palimpsest_input(vcf = ANNOT_VCF_DATA, Type = "ID"), 
        input_signatures  = ID_cosmic,
        resdir            = WorkDir, 
        signature_colours = signature_colour_generator(rownames(ID_cosmic)), 
        input_vcf         = ANNOT_VCF_DATA,
        doplot            = FALSE
    )

    SBS_Origins <- signature_origins( 
        input                  = ANNOT_VCF_DATA %>% filter( Type == "SNV"), 
        Type                   = "SBS", 
        input_signatures       = SBS_cosmic, 
        signature_contribution = SBS_SigExp
    )
    ID_Origins <- signature_origins( 
        input                  = ANNOT_VCF_DATA %>% filter( Type != "SNV"), 
        Type                   = "ID", 
        input_signatures       = ID_cosmic, 
        signature_contribution = ID_SigExp
    )

    openxlsx::write.xlsx( 
        list( COSMIC_SBS = SBS_Origins, COSMIC_IDS = ID_Origins ),
        sprintf("%s/%s.variants.mutsig.origin.assign.result.xlsx", WorkDir, BatchID ),
        rowNames=FALSE, overwrite=TRUE
    )

    MutSigOrigins <- rbind(
        SBS_Origins %>% mutate( index = paste(Sample, CHROM, POS, REF, ALT, sep="_") ) %>% 
            dplyr::select(c("index","SBS.Sig.max","SBS.Sig.max.prob")) %>% dplyr::rename( MutSig.Max=2, MutSig.Max.Prob=3 ),
        ID_Origins %>% mutate( index = paste(Sample, CHROM, POS, REF, ALT, sep="_") ) %>% 
            dplyr::select(c("index","ID.Sig.max","ID.Sig.max.prob")) %>% dplyr::rename( MutSig.Max=2, MutSig.Max.Prob=3 )
    ) %>% unique() 

    ANNOT_VCF_DATA <- ANNOT_VCF_DATA %>% 
        mutate( index = paste(Sample, CHROM, POS, REF, ALT, sep="_") ) %>% 
        left_join( ., MutSigOrigins, by = "index" ) %>% dplyr::select(!c("index"))

###
### Clonality Analysis
    #VCF_CNA <- cnaCCF_annot(vcf= VCF_DATA, annot_data = ANNOT_DATA, cna_data = CNV_DATA, CCF_boundary = 0.95)

    VCF_CNA <- cnaCCF_annot(
        vcf          = ANNOT_VCF_DATA, 
        annot_data   = ANNOT_DATA, 
        cna_data     = CNV_DATA %>% mutate(CHROM=gsub("^","chr", CHROM)), 
        CCF_boundary = 0.95
    )

    ChromDupTime <- chrTime_annot(
        vcf      = VCF_CNA, 
        cna_data = CNV_DATA %>% mutate(CHROM=gsub("^","chr", CHROM)), 
        cyto     = cytoband_hg19
    )

    PALIMPSEST_RESULT_LIST = list(
        VARIANT_RES  = ChromDupTime$vcf,
        MUT_TIME_RES = ChromDupTime$point.mut.time,
        CNV_RES      = ChromDupTime$cna_data
    )

    CLONALITY_RES <- ChromDupTime$vcf %>% mutate( index = paste( Sample, dup.ID, sep="_") ) %>% left_join(
        ., 
        ChromDupTime$point.mut.time %>% mutate( index = paste( Sample, ID, sep="_") ) %>% dplyr::select(c("index","Moltime")),
        by = "index"
    ) 

    write.table( 
        CLONALITY_RES, 
        sprintf("%s/%s.palimpsest.clonality.analysis.result.tsv", WorkDir, BatchID), 
        quote=F, col.names=T, row.names=F, sep="\t"
    )
    saveRDS( CLONALITY_RES, file=sprintf("%s/%s.palimpsest.clonality.analysis.result.RDS", WorkDir, BatchID ) )

    save( CLONALITY_RES, PALIMPSEST_RESULT_LIST, file=sprintf("%s/%s.palimpsest.clonality.analysis.result.objects.RData", WorkDir, BatchID ) )
###

