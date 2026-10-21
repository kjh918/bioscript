
###
    options(stringsAsFactors=FALSE) 
    suppressPackageStartupMessages({
        library("Hmisc")
        library("plyr")
        library("dplyr")
        library("setwidth")
        library("optparse")
        suppressWarnings(library("dndscv"))
        suppressWarnings(library("cancereffectsizeR"))
    })
###
    BaseDir="/data/wgs"
    BatchID="WGS-26-01"
    AnalysisList="/data/wgs/WGS-26-01/meta/analysis.list"
    WorkDir="/data/wgs/WGS-26-01/analysis/cancer-effect-size"

    sinfo <- read.table(sprintf("%s",AnalysisList), header=F, sep=" ") %>% dplyr::rename(SeqID=1, Group=2, Order=3)
    sinfo$SampleID <- make.names(sinfo$SeqID)
    tumor_ids      <- sinfo %>% filter( Group != "normal" ) %>% .$SeqID
    tumor_samples  <- sinfo %>% filter( Group != "normal" ) %>% .$SampleID

###
###
    CESA_RemoveClasses <- c("Intron","5'Flank","3'Flank","Targeted_Region","IGR","3'UTR","5'UTR","RNA")



    CESA_RES_List <- indexing( tumor_ids )
    CESA_RANK_LIST <- CESA_RES_List

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
            filter( alt_reads > 0, VariantClass %nin% CESA_RemoveClasses ) %>% 
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


    #merged_maf = data.frame()
    FILTERED_MERGED_MAF = data.frame()
    for( tid in tumor_ids ){    
        
        purple_dir <- sprintf("%s/%s/%s/hmftools/purple", BaseDir, BatchID, tid)
        maf_data <- read.delim(sprintf("%s/%s.purple.somatic.vep.annotation.maf",purple_dir,tid), skip=1)
        filtered_maf <- maf_data %>% filter( FILTER == "PASS" ) %>% 
            filter( Variant_Classification %nin% CESA_RemoveClasses ) %>% 
            filter( t_alt_count >= 2 ) %>% 
            filter( t_ref_count+t_alt_count >= 10 ) 
        
        FILTERED_MERGED_MAF = rbind(FILTERED_MERGED_MAF, filtered_maf)
    }

    write.table(FILTERED_MERGED_MAF,
        sprintf("%s/%s.cancereffectsizeR.input.filtered.merged.maf", WorkDir, BatchID),
        quote=F, col.names=T, row.names=F, sep="\t"
    )
    
        maf_read     <- preload_maf(sprintf("%s/%s.cancereffectsizeR.input.filtered.merged.maf", WorkDir, BatchID), refset = "ces.refset.hg19")
        maf_filtered <- maf_read[germline_variant_site == F][repetitive_region == F | cosmic_site_tier %in% 1:3]
    
        cesa_res <- CESAnalysis(refset = "ces.refset.hg19")
        cesa_res <- load_maf( cesa = cesa_res, 
            maf      = maf_filtered, 
            maf_name = tid, 
            coverage = "genome"
        )
        cesa_res <- gene_mutation_rates( cesa = cesa_res, 
            covariates = ces.refset.hg19$covariates$breast, 
            save_all_dndscv_output = TRUE
        )

        SigExclusions <- suggest_cosmic_signature_exclusions( cancer_type = "BRCA")
        cesa_res <- trinuc_mutation_rates( cesa = cesa_res, 
            signature_set = ces.refset.hg19$signatures$COSMIC_v3.2,
            signature_exclusions = SigExclusions,
            cores = 4
        )
        SelectedVariants <- select_variants( cesa = cesa_res )
        
        cesa_res <- ces_variant( cesa = cesa_res , 
            variants = SelectedVariants,  
            run_name="VarAnalysis", 
            conf = 0.95, 
            cores = 4
        )
        csea_res_ranked <- snv_results(cesa_res)[["VarAnalysis"]] %>% arrange(-selection_intensity) 

        CESA_RES_List[[ tid ]]  <- cesa_res
        CESA_RANK_LIST[[ tid ]] <- csea_res_ranked
    }


    save( CESA_RES_List, CESA_RANK_LIST, file=sprintf("%s/%s.cancereffectsizeR.RData", WorkDir, BatchID ) )
    write.table( CESA_RANK_LIST[[tumor_ids[1]]], sprintf("%s/%s.cancereffectsizeR.ranked.result.tsv", WorkDir, tumor_ids[1] ), quote=F, col.names=T, row.names=F, sep="\t")
    write.table( CESA_RANK_LIST[[tumor_ids[2]]], sprintf("%s/%s.cancereffectsizeR.ranked.result.tsv", WorkDir, tumor_ids[2] ), quote=F, col.names=T, row.names=F, sep="\t")
    
    