
    options(stringsAsFactors=FALSE) 
    suppressPackageStartupMessages(library("Hmisc"))
    suppressPackageStartupMessages(library("plyr"))
    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("setwidth"))
    suppressPackageStartupMessages(library("RMySQL"))
    suppressPackageStartupMessages(library("optparse"))
    #suppressPackageStartupMessages(library("intansv"))
    suppressPackageStartupMessages(library("VariantAnnotation"))
    suppressPackageStartupMessages(library("GenomicRanges"))

    hg38_Genes_GR = readRDS("/storage/home/kangsm/myDB/rds/hg38_all.genes.GRanges.RDS")

    SeqID       = "cbNIPT_PubDB_04_07_DS2"
    DataDir     = sprintf("/data/cbNIPT/processed_dataset/%s/sv", SeqID)
    MinSize_INV = 5000000 # 5 Mb
    MinSize_ITX = 1000000 # 1 Mb
    MinSize_INS = 100000 # 100 kb

    WorkDir = DataDir
    setwd(WorkDir)

    BreakdancerResultFile=sprintf("%s/breakdancer.ctx", DataDir)
    DellyInvVcfFile=sprintf("%s/delly.INV.filtered.vcf", DataDir)
    DellyBndVcfFile=sprintf("%s/delly.BND.filtered.vcf", DataDir)

    # breakdancer
    {
        SkipRows = as.numeric(system(sprintf("cat %s | grep ^# | wc -l", BreakdancerResultFile), intern=TRUE))-1
        bdctx    = read.delim(BreakdancerResultFile, skip=SkipRows) %>% dplyr::rename(Chr1=1) %>% 
            mutate( Size_Mb = round(as.numeric(Size)/1000000, 3) ) %>% 
            dplyr::select( c("Chr1","Pos1","Chr2","Pos2","Type","Size_Mb","Size","Score","num_Reads") ) %>% 
            dplyr::rename( Chrom1=1, Start=2, Chrom2=3, End=4)

        bdres_inv = bdctx %>% filter( abs(Size) >= MinSize_INV, Score >= 60 ) %>% filter( Type %in% c("INV") )
        bdres_itx = bdctx %>% filter( abs(Size) >= MinSize_ITX, Score >= 60 ) %>% filter( Type %in% c("ITX") )
        bdres_ctx = bdctx %>% filter( Score >= 95, num_Reads >= 2 ) %>% filter( Type %in% c("CTX") )
        # inversion
        bdres_inv_GR = GRanges(bdres_inv$Chrom1, IRanges(start=bdres_inv$Start, end=bdres_inv$End)) 
        BDRES_INV = as.data.frame(reduce(bdres_inv_GR)) %>% mutate( type="INV", Size_Mb=round(width/1000000, 3)) %>% dplyr::select(!c("strand"))
        # intra-chrom translocation
        bdres_itx_GR = GRanges(bdres_itx$Chrom1, IRanges(start=bdres_itx$Start, end=bdres_itx$End)) 
        BDRES_ITX = as.data.frame(reduce(bdres_itx_GR)) %>% mutate( type="ITX", Size_Mb=round(width/1000000, 3)) %>% dplyr::select(!c("strand"))
        # inter-chrom trnaslocation
        ctx_pairs = unique( bdres_ctx[,c("Chrom1","Chrom2")])
        chr1_list = gtools::mixedsort(unique(ctx_pairs$Chrom1))
        CTX_PAIRS = data.frame()
        for( ctg in chr1_list )
        {
            ctx_pair_ctg = ctx_pairs %>% filter( Chrom1 == ctg )
            chr2_index = ctx_pair_ctg$Chrom2
            names(chr2_index) = 1:nrow(ctx_pair_ctg)
            CTX_PAIRS = rbind(CTX_PAIRS, ctx_pair_ctg[names(gtools::mixedsort(chr2_index)), ])
        }        
        BDRES_CTX = data.frame()
        for( j in 1:nrow(CTX_PAIRS) )
        {
            ctx_j = bdres_ctx %>% filter( Chrom1 == CTX_PAIRS[j, "Chrom1"], Chrom2 == CTX_PAIRS[j, "Chrom2"] ) %>% arrange( Start )

            start_pos = ctx_j$Start
            names(start_pos) = 1:nrow(ctx_j)
            start_pos_groups = split(start_pos, c(0, cumsum(diff(start_pos) > 1000000 )))

            end_pos = ctx_j$End
            names(end_pos) = 1:nrow(ctx_j)
            end_pos_groups = split(end_pos, c(0, cumsum(diff(end_pos) > 1000000 )))

            if( identical( lapply(start_pos_groups, length), lapply(end_pos_groups, length) ) )
            {
                ctx_merged = lapply( start_pos_groups, function(w) 
                {
                    ctx_j[ as.numeric(names(w)), ] %>% group_by( Chrom1, Chrom2 ) %>% reframe(
                        Start = min(Start), End = min(End), Type="CTX", Score=min(Score), num_Reads=sum(num_Reads)
                    ) %>% dplyr::relocate(c(1,3,2,4,5,6,7)) %>% as.data.frame()
                }) %>% ldply() %>% dplyr::select(!c(".id"))
            }else{
                ctx_merged = ctx_j %>% dplyr::select(c(1:5,8,9))
            }
            BDRES_CTX = rbind(BDRES_CTX, ctx_merged)
        }
    }
    # delly2
    {
        delly_vcf_inv = readVcf(DellyInvVcfFile)
        delly_res_inv = data.frame(
            as.data.frame(rowRanges(delly_vcf_inv))[,c("seqnames","start","REF","ALT","QUAL","FILTER")], 
            as.data.frame(info(delly_vcf_inv))[,c("END","MAPQ","SVTYPE")],
            geno(delly_vcf_inv)[["DR"]],
            geno(delly_vcf_inv)[["DV"]]
        ) %>% 
            dplyr::rename( Chrom=seqnames, Start=start, End=END, DR=10,DV=11) %>% 
            mutate( LENGTH=End-Start+1 ) %>% 
            filter(
                MAPQ   >= 30,
                DV     >= 1,
                LENGTH >= MinSize_INV
            )
        delly_inv_GR = GRanges(delly_res_inv$Chrom, IRanges(start=delly_res_inv$Start, end=delly_res_inv$End)) 
        DELLY_INV = as.data.frame(reduce(delly_inv_GR)) %>% mutate( type="INV", Size_Mb=round(width/1000000, 3)) %>% dplyr::select(!c("strand"))

        delly_vcf_bnd = readVcf(DellyBndVcfFile)
        delly_res_bnd = data.frame(
            as.data.frame(rowRanges(delly_vcf_bnd))[,c("seqnames","start","REF","ALT","QUAL","FILTER")], 
            as.data.frame(info(delly_vcf_bnd))[,c("CHR2","POS2","MAPQ")],
            geno(delly_vcf_bnd)[["DR"]],
            geno(delly_vcf_bnd)[["DV"]]
        ) %>% 
            dplyr::rename( Chrom1=seqnames, Pos1=start, Chrom2=CHR2, Pos2=POS2, DR=10,DV=11 )
        delly_res_bnd$Type = apply(delly_res_bnd[,c("Chrom1","Chrom2")], 1, function(w) ifelse( w[1] == w[2], "ITX", "CTX"))
        delly_res_bnd$ALT  = unlist(delly_res_bnd$ALT)
        delly_ITX = delly_res_bnd %>% filter( Type == 'ITX' ) %>% dplyr::relocate(c(1,2,7,8,12,3,4,5,6,9,10,11)) %>% filter( MAPQ >= 30, DR >= 1, DV >= 1 )
        delly_CTX = delly_res_bnd %>% filter( Type == 'CTX' ) %>% dplyr::relocate(c(1,2,7,8,12,3,4,5,6,9,10,11)) %>% filter( MAPQ >= 30, DR >= 1, DV >= 1 )
    }

    BDRES_INV
    BDRES_ITX
    BDRES_CTX
    DELLY_INV
    delly_ITX
    delly_CTX


    filter( MAPQ >= 20 ) %>% 

    
    filter( MAPQ >= 20 ) 


    readDelly("/data/cbNIPT/processed_dataset/cbNIPT_PubDB_04_07_DS/sv/delly.INV.filtered.vcf")    