
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
        library("MutationTimeR")
        library("BiocParallel")
    })
    indexing = function(X) { sapply(unique(X), function(y) list(y)) }
    No.NA.Merge <- function(x) {   
        x = x[!is.na(x)]
        x = x[ x != "NA" ]
        x = x[ x != "" ]
        x2 = unique(x)
        if( length(x2) == 0 ){ v = "" }else{ v = paste0(x2, collapse=";") }
        return(v)
    }
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
    RunCores = 24

    VarDataDir <- sprintf("%s/%s/analysis/variant-sites-allele-count", BaseDir, BatchID)
    PycloneDir <- sprintf("%s/%s/analysis/pyclone-vi-variant-clusters", BaseDir, BatchID)
    WorkDir    <- sprintf("%s/%s/analysis/clonality-mutationTimeR", BaseDir, BatchID)
    system(sprintf("mkdir -p %s", WorkDir))

    sinfo <- read.table(AnalysisList, header=F, sep=" ") %>% dplyr::rename(SeqID=1, Group=2, Order=3)
    sinfo$SampleID <- make.names(sinfo$SeqID)
    tumor_ids      <- sinfo %>% filter( Group != "normal" ) %>% .$SeqID
    tumor_samples  <- sinfo %>% filter( Group != "normal" ) %>% .$SampleID

    IncludeClass <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins",
        "In_Frame_Del","In_Frame_Ins","Splice_Region","Splice_Site","Intron","Silent","5'UTR","3'UTR"
    )
    CodingClass <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Splice_Region","Splice_Site")
    Tag_Clonal = c("clonal [late]"="clonal","clonal [early]"="clonal","subclonal"="subclonal","clonal [NA]"="clonal", "none"="" )
    Tag_Timing = c("clonal [late]"="after_CN_gain","clonal [early]"="before_CN_gain","subclonal"="","clonal [NA]"="", "none"="" )

###
###
    wgs_var_count_list <- readRDS(sprintf("%s/%s.wgs.variants.read.count.list.RDS", VarDataDir, BatchID))
    wgs_pyclone_res    <- fread(sprintf("%s/%s.wgs.pyclone-vi.result.tsv", PycloneDir, BatchID))
    
    VarInfo <- ldply(wgs_var_count_list) %>% group_by(vkey) %>% reframe(
        Gene          = No.NA.Merge(gene),
        VariantClass  = No.NA.Merge(variant_class),
        VariantType   = No.NA.Merge(variant_type),
        HGVSp         = No.NA.Merge(hgvsp),
        Tier          = No.NA.Merge(TIER)
    ) %>% as.data.frame()

###
    for( tid in tumor_ids )
    {
        # tid <- tumor_ids[1]

        ## VCF
        {
            vcf_file <- sprintf("%s/%s.wgs.purple.somatic.filtered.vcf.gz", VarDataDir, tid)
            vcf_data <- readVcf(vcf_file)
            # fix AD
            ad_raw    <- geno(vcf_data)$AD
            ad_matrix <- matrix(unlist(ad_raw), ncol = 2, byrow = TRUE)
            ad_array  <- array(ad_matrix, dim = c(nrow(ad_matrix), 1, 2))
            dimnames(ad_array) <- list(NULL, colnames(vcf_data), c("ref", "alt"))
            geno(vcf_data)$AD  <- ad_array

            usable_vars <- wgs_var_count_list[[tid]] %>% mutate( 
                gene          = VarInfo[match(vkey, VarInfo$vkey), "Gene"],
                variant_class = VarInfo[match(vkey, VarInfo$vkey), "VariantClass"]
            ) %>% filter( alt_reads > 0, variant_class %in% IncludeClass ) %>% mutate(
                index = paste0(CHROM,":",POS,"_",REF,"/",ALT)
            )
        }
        MTR_INPUT_VCF_DATA <- vcf_data[intersect(rownames(vcf_data), usable_vars$index), ]
        
        ## CNV
        {
            purple_dir <- sprintf("%s/%s/%s/hmftools/purple", BaseDir, BatchID, tid)
            purple_cnv <- read.delim(sprintf("%s/%s.purple.cnv.somatic.tsv", purple_dir, tid))
            qc_data    <- read.delim(sprintf("%s/%s.purple.purity.tsv",      purple_dir, tid))
            purple_cnv$purity <- qc_data$purity

            # PURPLE 데이터를 mutationTimeR 포맷으로 정밀하게 변환
            refined_cnv <- purple_cnv %>% rowwise() %>% do({
                row   <- .
                major <- row$majorAlleleCopyNumber
                minor <- row$minorAlleleCopyNumber
                total_purity <- row$purity # 샘플 전체 순도
                
                # 정수부와 소수부 분리
                m_int <- floor(major)
                m_frac <- major - m_int
                
                if (m_frac > 0.05 && m_frac < 0.95) {
                    # 애매한 구간: 두 개의 행으로 분리
                    # Clone 1: 내림한 정수값, Clone 2: 올림한 정수값
                    data.frame(
                        chromosome = row$chromosome,
                        start      = row$start,
                        end        = row$end,
                        major_cn   = c(as.integer(m_int), as.integer(m_int + 1)),
                        minor_cn   = Rounding(minor, 0),
                        clonal_frequency = c(total_purity * (1 - m_frac), total_purity * m_frac) # 소수점 비율에 따라 전체 purity를 배분
                    )
                } else {
                    # 깨끗한 구간: 그대로 반올림
                    data.frame(
                        chromosome       = row$chromosome,
                        start            = row$start,
                        end              = row$end,
                        major_cn         = Rounding(major, 0),
                        minor_cn         = Rounding(minor, 0),
                        clonal_frequency = total_purity
                    )
                }
            }) %>% ungroup()
        }
        MTR_INPUT_CNV_DATA <- GRanges(
            seqnames         = refined_cnv$chromosome,
            ranges           = IRanges(refined_cnv$start, refined_cnv$end),
            major_cn         = refined_cnv$major_cn,
            minor_cn         = refined_cnv$minor_cn,
            clonal_frequency = refined_cnv$clonal_frequency
        )

        ## CLUSTERS
        pyc_data <- wgs_pyclone_res %>% filter( sample_id == tid ) %>% filter( mutation_id %in% usable_vars$vkey ) 
        
        MTR_INPUT_CLUSTERS <- pyc_data %>% group_by(cluster_id) %>% reframe(
            n_ssms = length(mutation_id)
        ) %>% mutate(
            proportion = Rounding( n_ssms/sum(n_ssms), 5 )
        ) %>% dplyr::rename(cluster=1) %>% as.data.frame()

        ## Run mutationTimeR
        register(MulticoreParam(workers = RunCores)) 
        
        mutTimeR <- mutationTime(
            vcf      = MTR_INPUT_VCF_DATA, 
            cn       = MTR_INPUT_CNV_DATA, 
            clusters = MTR_INPUT_CLUSTERS,
            purity   = qc_data$purity,
            gender   = tolower(qc_data$gender),
            isWgd    = as.logical(qc_data$wholeGenomeDuplication),
            n.boot   = 1000
        )

        # Results 
        vcf_res <- addMutTime(MTR_INPUT_VCF_DATA, mutTimeR$V)
        cnv_res <- cbind(as.data.frame(MTR_INPUT_CNV_DATA), mutTimeR$T) %>% as.data.frame() %>% dplyr::select(!c("timing_param"))

        vcf_res_raw <- cbind(as.data.table(rowRanges(vcf_res)), as.data.table(info(vcf_res)))
        vcf_res_raw$REF <- as.character(vcf_res_raw$REF)
        vcf_res_raw$ALT <- sapply(vcf_res_raw$ALT, function(x) paste(as.character(x), collapse = ","))
        vcf_res_raw <- vcf_res_raw %>% mutate( vkey = paste(as.character(seqnames), start, REF, ALT, sep="_") )

        MTR_RES <- usable_vars %>% filter( index %in% rownames(vcf_data) ) %>% 
            left_join(.,
                vcf_res_raw %>% mutate( vkey = paste(as.character(seqnames), start, REF, ALT, sep="_") ) %>% 
                    dplyr::select( c("vkey", colnames(mutTimeR$V)) ),
                by="vkey"           
            ) %>% dplyr::select(!c("index")) %>% mutate( CLS=as.character(CLS) ) 

        MTR_RES[which(is.na(MTR_RES$CLS)), "CLS"] = "none"
        MTR_RES$clonality <- Tag_Clonal[ MTR_RES$CLS ]
        MTR_RES$timing    <- Tag_Timing[ MTR_RES$CLS ]

        MTR_RES <- MTR_RES %>% left_join(., 
            pyc_data %>% as.data.frame() %>% dplyr::rename( vkey=1 ) %>% dplyr::select(c(1,3,4)), 
            by="vkey" 
        ) 

        vcf_gr <- GRanges(MTR_RES$CHROM,    IRanges(MTR_RES$POS, width=1) )
        cnv_gr <- GRanges(cnv_res$seqnames, IRanges(cnv_res$start, cnv_res$end) )
        ovlp   <- findOverlaps(vcf_gr, cnv_gr)

        cnv_res_cols <- which(colnames(cnv_res)=="clonal_frequency") : ncol(cnv_res)
        MTR_RES <- cbind(MTR_RES[queryHits(ovlp), ], cnv_res[subjectHits(ovlp), cnv_res_cols])

        col_types <- sapply(MTR_RES, class)
        list_cols <- names(col_types[col_types == "SimpleList" | col_types == "list"])
        for (col in list_cols) {
            MTR_RES[, col] <- sapply(MTR_RES[, col], function(x) paste(as.character(x), collapse = ","))
        }

        # cluster summary
        MTR_CLUSTERS <- list(
            all_variants = MTR_RES %>% group_by(cluster_id,clonality,timing) %>% reframe(
                n_mut      = length(vkey),
                mean_pGain = Rounding(mean(pGain, na.rm = TRUE), 3),
                mean_pSub  = Rounding(mean(pSub,  na.rm = TRUE), 3),
                mean_MutCN = Rounding(mean(MutCN, na.rm = TRUE), 3),
                mean_cellular_prevalence = mean(cellular_prevalence, na.rm=TRUE)
            ) %>% group_by(cluster_id) %>% mutate( mut_fraction = Rounding(n_mut/sum(n_mut),3) ) %>% as.data.frame(),
            coding_variants = MTR_RES %>% filter( variant_class %in% CodingClass ) %>% group_by(cluster_id,clonality,timing) %>% reframe(
                n_mut      = length(vkey),
                mean_pGain = Rounding(mean(pGain, na.rm = TRUE), 3),
                mean_pSub  = Rounding(mean(pSub,  na.rm = TRUE), 3),
                mean_MutCN = Rounding(mean(MutCN, na.rm = TRUE), 3),
                mean_cellular_prevalence = mean(cellular_prevalence, na.rm=TRUE)
            ) %>% group_by(cluster_id) %>% mutate( mut_fraction = Rounding(n_mut/sum(n_mut),3) ) %>% as.data.frame()
        )

        # save results
        write.table( MTR_RES, 
            sprintf("%s/%s.mutationTimeR.result.tsv", WorkDir, tid), 
            quote=F, col.names=T, row.names=F, sep="\t"
        )
        openxlsx::write.xlsx( MTR_CLUSTERS,
            sprintf("%s/%s.mutationTimeR.result.cluster.summary.xlsx", WorkDir, tid), 
            rowNames=F, overwrite=T
        )
        save(MTR_RES, vcf_res, cnv_res, MTR_CLUSTERS, file=sprintf("%s/%s.mutationTimeR.analysis.result.objects.RData", WorkDir, tid))

    }

###
    # "MutCN" : 해당 mutation이 몇 개의 allele에 존재하는지
    # "MutDeltaCN" :  mutation이 존재하는 allele의 CN 변화량
    # "MajCN"/"MinCN" :  major / minor allele copy number
    # "MajDerCN"/"MinDerCN" :  tumor-derived CN (normal 제외)
    # "CNF" : copy number fraction / segment fraction
    # "CNID" :  copy number segment ID
    # "pMutCN" : MutCN 추정의 posterior probability
    # "pGain" : mutation이 CN gain 이전에 발생했을 확률, pGain 크면 early
    # "pSingle" : duplication 이후 mutation일 가능성
    # "pSub" : subclonal mutation일 확률
    # "pAllSubclones" : 모든 subclone에서 존재할 확률
    # "pMutCNTail"  불확실성/이상값 indicator , 값 크면 unreliable 가능성
    #  "CLS" 
###
    # plotSample(vcf_data,cnv_data)
###
