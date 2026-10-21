
# 3. Merge WGS-TSO500 VARIANTS COUNTS

### packages & functions
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
    Rounding <- function(x, digits = 0) {   
        posneg <- sign(x)
        z <- abs(x) * 10^digits
        z <- floor(z + 0.5)
        v <- z / 10^digits * posneg
        return(v)
    }
###

### 
    BaseDir <- "/data/wgs/RUO_Service"
    BatchID <- "WGS-26-02"
    GenomeAssembly <- "hg19"
    AnalysisList   <- sprintf("%s/%s/meta/clonality.evolution.analysis.list", BaseDir, BatchID)
    VarDataDir     <- sprintf("%s/%s/analysis/variant-sites-allele-count", BaseDir, BatchID)

    sinfo <- read.table(AnalysisList, header=F, sep=" ") %>% dplyr::rename(SeqID=1, Group=2, Order=3)
    sinfo$SampleID <- make.names(sinfo$SeqID)
    tumor_ids      <- sinfo %>% filter( Group != "normal" ) %>% .$SeqID
    tumor_samples  <- sinfo %>% filter( Group != "normal" ) %>% .$SampleID
###

###
    if( GenomeAssembly == "hg19" ){
        suppressPackageStartupMessages({library("BSgenome.Hsapiens.UCSC.hg19")})
        GENOME=BSgenome.Hsapiens.UCSC.hg19
    }else{
        suppressPackageStartupMessages({library("BSgenome.Hsapiens.UCSC.hg38")})
        GENOME=BSgenome.Hsapiens.UCSC.hg38
    }
###

### var-data
    var_data_wgs <- read.delim(sprintf("%s/%s.wgs.tumors.variants.vaf.tsv", VarDataDir, BatchID)) %>% 
        mutate( vkey=paste(CHROM,POS,REF,ALT, sep="_") ) 
    var_data_tso <- read.delim(sprintf("%s/%s.tso.tumors.variants.vaf.tsv", VarDataDir, BatchID)) %>% 
        mutate( vkey=paste(CHROM,POS,REF,ALT, sep="_") ) 
###

### allele-count-data
    allele_count_wgs <- lapply(indexing(tumor_ids), function(tid) 
        return(read.delim(sprintf("%s/%s.wgs.variant.sites.bam.counts.summary.tsv", VarDataDir, tid)) %>% mutate(index=paste0(chrom, "_", pos)))
    )
    allele_count_tso <- lapply(indexing(tumor_ids), function(tid) 
        return(read.delim(sprintf("%s/%s.tso.variant.sites.bam.counts.summary.tsv", VarDataDir, tid)) %>% mutate(index=paste0(chrom, "_", pos)))
    )
###

### FIX INDEL ALLELES
    # wgs
    vep_maf_wgs <- lapply( indexing(tumor_ids), function(tid) {
        maf_wgs <- sprintf("%s/%s/%s/hmftools/purple/%s.purple.somatic.vep.annotation.maf", BaseDir, BatchID, tid, tid)
        vep_maf <- read.delim(maf_wgs, skip=1) 
        vep_maf$REF <- apply(vep_maf[,c("Chromosome","vcf_pos","Variant_Type","Reference_Allele")], 1, function(v) {
            if( v[3] == "INS" ){
                ref_allele = as.character(GENOME[[paste0("chr",v[1])]][as.numeric(v[2]):as.numeric(v[2])])
            }else if( v[3] == "DEL" ){
                ref_allele = paste0(as.character(GENOME[[paste0("chr",v[1])]][as.numeric(v[2]):as.numeric(v[2])]), v[4])
            }else{
                ref_allele = v[4]
            }
            return(ref_allele)
        } )
        vep_maf$ALT <- apply(vep_maf[,c("Chromosome","vcf_pos","Variant_Type","Tumor_Seq_Allele2")], 1, function(v) {
            if( v[3] == "INS" ){
                alt_allele = paste0(as.character(GENOME[[paste0("chr",v[1])]][as.numeric(v[2]):as.numeric(v[2])]), v[4])
            }else if( v[3] == "DEL" ){
                alt_allele = as.character(GENOME[[paste0("chr",v[1])]][as.numeric(v[2]):as.numeric(v[2])])
            }else{
                alt_allele = v[4]
            }
            return(alt_allele)
        } )
        vep_maf_rev <- vep_maf %>% mutate( vkey=paste(Chromosome,vcf_pos,REF,ALT, sep="_") ) %>% 
            mutate( vaf = Rounding(t_alt_count/(t_alt_count+t_ref_count), 3)) %>% dplyr::select(
            c("vkey","Hugo_Symbol","vcf_pos","REF","ALT","Variant_Classification","Variant_Type","vaf","HGVSp_Short","n_depth","FILTER","TIER")
        ) %>% dplyr::rename(
            gene=Hugo_Symbol, pos=vcf_pos, ref=REF, alt=ALT, variant_class=Variant_Classification, variant_type=Variant_Type,
            hgvsp=HGVSp_Short, matched_normal_dp=n_depth
        )
        return(vep_maf_rev)
    })

    # tso
    vep_maf_tso <- lapply( indexing(tumor_ids), function(tid) {
        maf_tso <- sprintf("%s/%s/%s/tso500/tso_hmftools/08.purple/%s.purple.somatic.vep.annotation.maf", BaseDir, BatchID, tid, tid)
        vep_maf <- read.delim(maf_tso, skip=1) 
        vep_maf$REF <- apply(vep_maf[,c("Chromosome","vcf_pos","Variant_Type","Reference_Allele")], 1, function(v) {
            if( v[3] == "INS" ){
                ref_allele = as.character(GENOME[[paste0("chr",v[1])]][as.numeric(v[2]):as.numeric(v[2])])
            }else if( v[3] == "DEL" ){
                ref_allele = paste0(as.character(GENOME[[paste0("chr",v[1])]][as.numeric(v[2]):as.numeric(v[2])]), v[4])
            }else{
                ref_allele = v[4]
            }
            return(ref_allele)
        } )
        vep_maf$ALT <- apply(vep_maf[,c("Chromosome","vcf_pos","Variant_Type","Tumor_Seq_Allele2")], 1, function(v) {
            if( v[3] == "INS" ){
                alt_allele = paste0(as.character(GENOME[[paste0("chr",v[1])]][as.numeric(v[2]):as.numeric(v[2])]), v[4])
            }else if( v[3] == "DEL" ){
                alt_allele = as.character(GENOME[[paste0("chr",v[1])]][as.numeric(v[2]):as.numeric(v[2])])
            }else{
                alt_allele = v[4]
            }
            return(alt_allele)
        } )
        vep_maf_rev <- vep_maf %>% mutate( vkey=paste(Chromosome,vcf_pos,REF,ALT, sep="_") ) %>% 
            mutate( vaf = Rounding(t_alt_count/(t_alt_count+t_ref_count), 3)) %>% dplyr::select(
            c("vkey","Hugo_Symbol","vcf_pos","REF","ALT","Variant_Classification","Variant_Type","vaf","HGVSp_Short","FILTER","TIER")
        ) %>% dplyr::rename(
            gene=Hugo_Symbol, pos=vcf_pos, ref=REF, alt=ALT, variant_class=Variant_Classification, variant_type=Variant_Type,
            hgvsp=HGVSp_Short
        )
        return(vep_maf_rev)
    })
###

### WGS VARIANTS
    VarCount_WGS_List <- lapply( indexing(tumor_ids), function(tid) {
        
        VD <- var_data_wgs[,c("CHROM","POS","REF","ALT", make.names(tid))]
        VD$index <- apply(VD, 1, function(v) {
            if( nchar(v[3]) > 1 ){ 
                idx <- paste0( v[1], "_", as.numeric(v[2])+1 )
            }else{
                idx <- paste0( v[1], "_", as.numeric(v[2]) )
            }
            return(idx)
        })
        VD <- VD %>% left_join(
            ., 
            allele_count_wgs[[tid]] %>% unique(),
            by="index"
        )

        VD$ref_reads = apply( VD[,c(make.names(tid), "ref_count")], 1, function(ad) {
            if( ad[1] == "." ){ rd <- as.numeric(ad[2]) }else{ rd <- as.numeric(unlist(strsplit(ad[1], ","))[1]) }
            return( as.numeric(rd) )
        })
        VD$alt_reads <- apply(VD, 1, function(ad) {
            if( ad[5] == "." ){
                if( nchar(ad[3]) == nchar(ad[4]) ){ 
                    vt = 'snv' 
                    va = ad[4]
                }else if( nchar(ad[3]) > nchar(ad[4]) ){
                    vt = 'del'
                    va = paste0("-", gsub("^[A-Z]", "", ad[3]))
                }else if( nchar(ad[3]) < nchar(ad[4]) ){
                    vt = 'ins'
                    va = paste0("+", gsub("^[A-Z]", "", ad[4]))
                }
                if( vt == 'snv' ) { 
                    allele_cols <- c("A"=11, "T"=12, "G"=13, "C"=14)
                    rd <- ad[ allele_cols[ ad[4] ] ]
                }
                if( vt %in% c('ins','del') ) {  
                    alt_a = unlist(strsplit(ad[19], ","))
                    if( length(alt_a) > 0 ){
                        alt_df = data.frame(
                            allele = sapply(alt_a, function(x) unlist(strsplit(x, ":"))[1]),
                            count  = sapply(alt_a, function(x) unlist(strsplit(x, ":"))[2])
                        )
                        alt_c = alt_df[which(alt_df$allele == va), "count"]
                        if( length(alt_c) == 0 ){ alt_c = 0 }
                        if( alt_c == "" ){ rd = 0 }else{ rd = alt_c }
                    }else{
                        rd = 0
                    }                  
                }
            }else{
                rd = unlist(strsplit(ad[5], ","))[2]
            }
            return( as.numeric(rd) )
        })

        VD_RES <- VD %>% mutate( total_reads = ref_reads+alt_reads ) %>% 
            mutate( vkey=paste(CHROM,POS,REF,ALT, sep="_") ) %>% 
            left_join(., vep_maf_wgs[[tid]], by="vkey") %>% dplyr::select(
                c("vkey","CHROM","POS","REF","ALT","ref_reads","alt_reads","total_reads",
                "gene","variant_class","variant_type","vaf","hgvsp","matched_normal_dp","FILTER","TIER"
        ))

        return(VD_RES)
    })
###

### TSO VARIANTS
    VarCount_TSO_List <- lapply( indexing(tumor_ids), function(tid) {
        VD <- var_data_tso[,c("CHROM","POS","REF","ALT", make.names(tid))]
        VD$index <- apply(VD, 1, function(v) {
            if( nchar(v[3]) > 1 ){ 
                idx <- paste0( v[1], "_", as.numeric(v[2])+1 )
            }else{
                idx <- paste0( v[1], "_", as.numeric(v[2]) )
            }
            return(idx)
        })
        VD <- VD %>% left_join(
            ., 
            allele_count_tso[[tid]] %>% unique(),
            by="index"
        )
        VD$ref_reads = apply( VD[,c(make.names(tid), "ref_count")], 1, function(ad) {
            if( ad[1] == "." ){ rd <- as.numeric(ad[2]) }else{ rd <- as.numeric(unlist(strsplit(ad[1], ","))[1]) }
            return( as.numeric(rd) )
        })

        VD$alt_reads <- apply(VD, 1, function(ad) {
            if( ad[5] == "." ){
                if( nchar(ad[3]) == nchar(ad[4]) ){ 
                    vt = 'snv' 
                    va = ad[4]
                }else if( nchar(ad[3]) > nchar(ad[4]) ){
                    vt = 'del'
                    va = paste0("-", gsub("^[A-Z]", "", ad[3]))
                }else if( nchar(ad[3]) < nchar(ad[4]) ){
                    vt = 'ins'
                    va = paste0("+", gsub("^[A-Z]", "", ad[4]))
                }
                if( vt == 'snv' ) { 
                    allele_cols <- c("A"=11, "T"=12, "G"=13, "C"=14)
                    rd <- ad[ allele_cols[ad[4]]  ]
                }
                if( vt %in% c('ins','del') ) {  
                    alt_a = unlist(strsplit(ad[19], ","))
                    alt_df = data.frame(
                        allele = sapply(alt_a, function(x) unlist(strsplit(x, ":"))[1]),
                        count  = sapply(alt_a, function(x) unlist(strsplit(x, ":"))[2])
                    )
                    alt_c = alt_df[which(alt_df$allele == va), "count"]
                    if( length(alt_c) == 0 ){ alt_c = 0 }
                    rd = alt_c
                }
            }else{
                rd = unlist(strsplit(ad[5], ","))[2]
            }
            return( as.numeric(rd) )
        })

        VD_RES <- VD %>% mutate( total_reads = ref_reads+alt_reads ) %>% 
            mutate( vkey=paste(CHROM,POS,REF,ALT, sep="_") ) %>% 
            left_join(., vep_maf_tso[[tid]], by="vkey") %>% dplyr::select(
                c("vkey","CHROM","POS","REF","ALT","ref_reads","alt_reads","total_reads",
                "gene","variant_class","variant_type","vaf","hgvsp","FILTER","TIER"
        ))
        return(VD_RES)
    })
###

### WGS VARIANTS AT TSO500 SITES 

    allele_count_wgs_tso_site <- lapply(indexing(tumor_ids), function(tid) 
        return(read.delim(sprintf("%s/%s.wgs.variant.tso.sites.bam.counts.summary.tsv", VarDataDir, tid)) %>% mutate(index=paste0(chrom, "_", pos)))
    )
    wgs_vkey = Reduce( union, lapply(VarCount_WGS_List, function(y) y$vkey))

    VarCount_WGS_TSO_Merged <- lapply( indexing(tumor_ids), function(tid) {
        VD <- var_data_tso[,c("CHROM","POS","REF","ALT", make.names(tid))]
        VD$index <- apply(VD, 1, function(v) {
            if( nchar(v[3]) > 1 ){ 
                idx <- paste0( v[1], "_", as.numeric(v[2])+1 )
            }else{
                idx <- paste0( v[1], "_", as.numeric(v[2]) )
            }
            return(idx)
        }) 
        VD <- VD %>% left_join(
            ., 
            allele_count_wgs_tso_site[[tid]] %>% unique(),
            by="index"
        )

        VD$ref_reads = VD$ref_count
        
        VD$alt_reads <- unlist(apply(VD, 1, function(ad) {
            if( ad[5] == "." ){
                if( nchar(ad[3]) == nchar(ad[4]) ){ 
                    vt = 'snv' 
                    va = ad[4]
                }else if( nchar(ad[3]) > nchar(ad[4]) ){
                    vt = 'del'
                    va = paste0("-", gsub("^[A-Z]", "", ad[3]))
                }else if( nchar(ad[3]) < nchar(ad[4]) ){
                    vt = 'ins'
                    va = paste0("+", gsub("^[A-Z]", "", ad[4]))
                }
                if( vt == 'snv' ) { 
                    allele_cols <- c("A"=11, "T"=12, "G"=13, "C"=14)
                    rd <- ad[ allele_cols[ad[4]]  ]
                }
                if( vt %in% c('ins','del') ) {  
                    alt_a = unlist(strsplit(ad[19], ","))
                    alt_df = data.frame(
                        allele = sapply(alt_a, function(x) unlist(strsplit(x, ":"))[1]),
                        count  = sapply(alt_a, function(x) unlist(strsplit(x, ":"))[2])
                    )
                    alt_c = alt_df[which(alt_df$allele == va), "count"]
                    if( length(alt_c) == 0 ){ alt_c = 0 }
                    rd = alt_c
                }
            }else{
                allele_cols <- c("A"=11, "T"=12, "G"=13, "C"=14)
                rd <- ad[ allele_cols[ad[4]]  ]
            }
            return( as.numeric(rd) )
        }))

        VD <- VD %>% mutate( vkey=paste(CHROM,POS,REF,ALT, sep="_") ) %>% filter( vkey %nin% wgs_vkey ) %>% mutate(keep = TRUE)
        VD[which(VD[, make.names(tid)] != "." & VD$alt_reads == 0), "keep"] <- FALSE 

        VD_RES <- VD %>% mutate( total_reads = ref_reads+alt_reads ) %>%
            left_join(., vep_maf_tso[[tid]], by="vkey") %>% dplyr::select(
                c("vkey","CHROM","POS","REF","ALT","ref_reads","alt_reads","total_reads",
                "gene","variant_class","variant_type","vaf","hgvsp","FILTER","TIER","keep"
        ))
        VD_RES[which(VD_RES$alt_reads > 0 & VD_RES$vaf < 0.03), "keep"] <- FALSE       
        VD_RES <- VD_RES %>% filter( keep ) %>% mutate( vaf = Rounding( alt_reads/total_reads, 3) ) %>% dplyr::select(!c("keep"))

        VD_MERGED <- rbind(VarCount_WGS_List[[tid]], VD_RES)
        
        return(VD_MERGED)
    })
###
### SAVE as RDS
    saveRDS( VarCount_WGS_List, file=sprintf("%s/%s.wgs.variants.read.count.list.RDS", VarDataDir, BatchID))
    saveRDS( VarCount_TSO_List, file=sprintf("%s/%s.tso.variants.read.count.list.RDS", VarDataDir, BatchID))
    saveRDS( VarCount_WGS_TSO_Merged, file=sprintf("%s/%s.wgs.tso.merged.variants.read.count.list.RDS", VarDataDir, BatchID))
###

```

###