
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
    })
    indexing = function(X) { sapply(unique(X), function(y) list(y)) }
###
### 
    BaseDir="/data/wgs"
    BatchID="WGS-26-01"
    AnalysisList="/data/wgs/WGS-26-01/meta/analysis.list"
    
    VarDataDir <- sprintf("%s/%s/analysis/variant-sites-allele-count", BaseDir, BatchID)
    WorkDir    <- sprintf("%s/%s/analysis/phylowgs", BaseDir, BatchID)
    system(sprintf("mkdir -p %s", WorkDir))

    sinfo <- read.table(AnalysisList, header=F, sep=" ") %>% dplyr::rename(SeqID=1, Group=2, Order=3)
    sinfo$SampleID <- make.names(sinfo$SeqID)
    tumor_ids      <- sinfo %>% filter( Group != "normal" ) %>% .$SeqID
    tumor_samples  <- sinfo %>% filter( Group != "normal" ) %>% .$SampleID
###
###
    var_data <- read.delim(sprintf("%s/%s.tumors.variants.vaf.tsv", VarDataDir, BatchID)) %>% 
        mutate( vkey=paste(CHROM,POS,REF,ALT, sep="_") ) 

    var_data_list <- lapply(indexing(tumor_ids), function(tid) {
        id <- sinfo %>% filter( SeqID == tid ) %>% .$SampleID
        maf_file <- sprintf("%s/%s/%s/hmftools/purple/%s.purple.somatic.vep.annotation.maf", BaseDir, BatchID, tid, tid)
        maf_data <- read.delim(maf_file, skip=1)
        cnt_data <- read.delim(sprintf("%s/%s.variant.sites.bam.counts.summary.tsv", VarDataDir, tid)) %>% mutate(index=paste0(chrom, "_", pos))

        vd <- var_data[,c("CHROM","POS","REF","ALT", id)] %>% mutate(index=paste0(CHROM, "_", POS)) %>% 
            left_join(., cnt_data[,c("index","ref_count")], by="index")
        
        vd$ref_reads = apply( vd[,c(id, "ref_count")], 1, function(ad){
            if( ad[1] == "." ){ 
                rd <- as.numeric(ad[2])
            }else{
                rd <- as.numeric(unlist(strsplit(ad[1], ","))[1])
            }
            return(rd)
        })
        vd$alt_reads = sapply( vd[,id], function(ad){
            if( ad == "." ){ 
                rd <- 0
            }else{
                rd <- as.numeric(unlist(strsplit(ad, ","))[2])
            }
            return(rd)
        })
        vd <- vd %>% left_join(
            ., maf_data %>% mutate(index=paste0(Chromosome, "_", Start_Position)) %>% dplyr::select(c("index","Variant_Classification")), by="index"
        )
        return(vd)
    })
###
### cnv-data
    cnv_data_list <- lapply( indexing(tumor_ids), function(tid) {
        purple_dir = sprintf("%s/%s/%s/hmftools/purple", BaseDir, BatchID, tid)
        cnv_data = read.delim(sprintf("%s/%s.purple.cnv.somatic.tsv",purple_dir,tid))
        qc_data   = read.delim(sprintf("%s/%s.purple.purity.tsv",purple_dir,tid))
        cnv_data$purity = qc_data$purity
        return(cnv_data)
    })
###
### create pyclone-vi input
    # c("VarIndex","sample_id","ref_counts","var_counts","major_cn","minor_cn","normal_cn","tumour_content")
    pyc_input_data_list <- lapply( indexing(tumor_ids), function(tid) {
        VARS   <- data.table(var_data_list[[tid]])
        CNVS   <- data.table(cnv_data_list[[tid]])
        PURITY <- unique(CNVS$purity)
        v_gr   <- GRanges( VARS$CHROM, IRanges(start=VARS$POS, end=VARS$POS), sample_id=tid, var_key=paste(VARS$CHROM,VARS$POS, VARS$REF, VARS$ALT, sep=":"))
        c_gr   <- GRanges( CNVS$chromosome, IRanges(start=CNVS$start, end=CNVS$end), sample_id=tid, majorAlleleCopyNumber = CNVS$majorAlleleCopyNumber, minorAlleleCopyNumber = CNVS$minorAlleleCopyNumber)
        hits   <- findOverlaps(v_gr, c_gr)
        VARS[, majorAlleleCopyNumber := 1]
        VARS[, minorAlleleCopyNumber := 1]
        VARS[, purity := PURITY]
        VARS[, normalCopyNumber := 2]
        if (length(hits) > 0) {
            hit_dt <- data.table(
                q = queryHits(hits),
                s = subjectHits(hits)
            )
            # if multiple CNV segments overlap a point, keep first
            hit_dt <- hit_dt[!duplicated(q)]

            VARS[hit_dt$q, majorAlleleCopyNumber := mcols(c_gr)$majorAlleleCopyNumber[hit_dt$s]]
            VARS[hit_dt$q, minorAlleleCopyNumber := mcols(c_gr)$minorAlleleCopyNumber[hit_dt$s]]
        } 

        RES <- VARS %>% mutate( 
            mutation_id = paste(VARS$CHROM,VARS$POS, VARS$REF, VARS$ALT, sep="_"),
            sample_id   = tid,
            majorAlleleCopyNumber = round(majorAlleleCopyNumber, 0),
            minorAlleleCopyNumber = round(minorAlleleCopyNumber, 0)
        ) %>% filter( 
            ref_reads+alt_reads > 10
        ) %>% dplyr::select(
            c("mutation_id","sample_id","ref_reads","alt_reads","majorAlleleCopyNumber","minorAlleleCopyNumber","normalCopyNumber","purity")
        ) %>% dplyr::rename(
            ref_counts=ref_reads, alt_counts=alt_reads, major_cn=majorAlleleCopyNumber, minor_cn=minorAlleleCopyNumber, normal_cn=normalCopyNumber, tumour_content=purity
        ) 
        return(RES)
    })
    PYC_INPUT_DATA <- pyc_input_data_list %>% ldply() %>% dplyr::select(!c(".id"))
    write.table(PYC_INPUT_DATA, sprintf("%s/%s.pyclone-vi.input.tsv", WorkDir, BatchID), quote=F, col.names=T, row.names=F, sep="\t")
###
### run pyclone-vi
    PYCLONE_VI_RUN_CMD <- "singularity exec -B /storage,/data /storage/images/pyclone-vi_v0.1.2.sif pyclone-vi"
    PYCLINE_VI_INPUT   <- sprintf("%s/%s.pyclone-vi.input.tsv", WorkDir, BatchID)
    PYCLONE_VI_OUTPUT  <- sprintf("%s/%s.pyclone-vi.output.h5", WorkDir, BatchID)
    PYCLONE_VI_RESULT  <- sprintf("%s/%s.pyclone-vi.result.tsv", WorkDir, BatchID)
    PYC_MAX_CLUSTERS   <- 20

    system(sprintf(
        "%s fit --in-file %s --out-file %s --num-clusters %s --density beta-binomial --num-restarts 10 --seed 1234", 
        PYCLONE_VI_RUN_CMD, PYCLINE_VI_INPUT, PYCLONE_VI_OUTPUT, PYC_MAX_CLUSTERS
    ))
    system(sprintf(
        "%s write-results-file --in-file %s --out-file %s", 
        PYCLONE_VI_RUN_CMD, PYCLONE_VI_OUTPUT, PYCLONE_VI_RESULT
    ))
    pyc_data    <- read.delim(PYCLONE_VI_RESULT)
    pyc_cluster <- pyc_data[,c("mutation_id","cluster_id")] %>% unique()
###
### 
    RemoveClasses <- c("Silent","Intron","5'Flank","3'Flank","Targeted_Region","IGR","3'UTR","5'UTR","RNA")

    filtered_var_data_list <- lapply(indexing(tumor_ids), function(tid) {
        pid <- pyc_input_data_list[[tid]]
        vd  <- var_data_list[[tid]] %>% mutate( mutation_id=paste(CHROM,POS, REF, ALT, sep="_") )       
        VARS <- pid %>% left_join(., vd[,c("mutation_id","Variant_Classification")], by="mutation_id") %>% 
            left_join(., pyc_cluster, by="mutation_id") %>% 
            filter( Variant_Classification %nin% RemoveClasses ) %>% 
            filter( !is.na(Variant_Classification) )
        var_bed <- data.frame(
            chrom=sapply(VARS$mutation_id, function(x) unlist(strsplit(x, "_"))[1]),
            start=sapply(VARS$mutation_id, function(x) as.numeric(unlist(strsplit(x, "_"))[2])-1),
            end=sapply(VARS$mutation_id, function(x) as.numeric(unlist(strsplit(x, "_"))[2]))
        )
        write.table(var_bed, sprintf("%s/%s.filtered.var.sites.bed", WorkDir, tid), quote=F, col.names=F, row.names=F, sep="\t")
        return(VARS)
    })     

    for( tid in tumor_ids ){
        input_vcf      <- sprintf("%s/%s.purple.somatic.filtered.vcf.gz", WorkDir, tid)
        output_vcf     <- sprintf("%s/%s.phylowgs.input.filtered.vcf.gz", WorkDir, tid)
        var_site_bed   <- sprintf("%s/%s.filtered.var.sites.bed", WorkDir, tid)
        VCF_FILTER_CMD <- sprintf("bcftools view -R %s -Oz -o %s %s", var_site_bed, output_vcf, input_vcf )
        system(VCF_FILTER_CMD)
    }


### cnv_data for 'create_phylowgs_inputs.py'
    for( tid in tumor_samples ){
        
        purple_dir = sprintf("%s/%s/%s/hmftools/purple", BaseDir, BatchID, tid)
        
        cnv_data = read.delim(sprintf("%s/%s.purple.cnv.somatic.tsv",purple_dir,tid))
        qc_data   = read.delim(sprintf("%s/%s.purple.purity.tsv",purple_dir,tid))

        cnv_res <- cnv_data %>% mutate( minor_cn=round(minorAlleleCopyNumber, 0), major_cn=round(majorAlleleCopyNumber, 0)) %>% 
            mutate( copy_number = minor_cn+major_cn, cellular_prevalence=qc_data$purity ) %>% 
            dplyr::select(
                c("chromosome","start","end","copy_number","minor_cn","major_cn","cellular_prevalence")
            ) 
        write.table(cnv_res,
            sprintf("%s/%s.cnv.data.tsv", work_dir, tid),
            quote=F, col.names=T, row.names=F, sep="\t"
        )
    }
###


        



    }
