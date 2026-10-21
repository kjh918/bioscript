
    library(karyoploteR)
    library(GenomicRanges)
    library(dplyr)
    library(RMySQL)
    library("Cairo")
    grDevices::X11.options(type='cairo')
    options(device='x11')
    options(bitmapType='cairo')

    BaseDir = "/data/wgs/RUO_Service"
    BatchID = "WGS-26-02"
    SeqID   = "WGS-26-02_02"
    SampleName = "BRCA-02-BM"

    res_dir = sprintf("%s/%s/analysis/summary_and_report", BaseDir, BatchID)


    cnv_gene_data = read.delim("//data/wgs/RUO_Service/WGS-26-02/WGS-26-02_02/hmftools/linx/WGS-26-02_02.linx.driver.catalog.tsv") %>% 
        filter( driver %in% c("AMP","DEL","DISRUPTION"))
    gene_list <- paste0(paste0("'", cnv_gene_data$gene, "'"), collapse=",")
    dbCon = dbConnect(dbDriver("MySQL"), host='192.168.0.34', user='gcx', port=3306, password='gencurix!!', db='public_database' )
    cnv_genes_info <- dbGetQuery(dbCon, sprintf("SELECT gene_name,chr,gene_start,gene_end,strand FROM hg19_gencode_v19_all_transcripts_seq WHERE gene_name IN (%s)", gene_list)) %>% 
        unique()
    dbDisconnect(dbCon)
    gene_gr = GRanges(cnv_genes_info$chr, IRanges(start=cnv_genes_info$gene_start, end=cnv_genes_info$gene_end), strand=cnv_genes_info$strand, gene=cnv_genes_info$gene_name)

    ## 예시 입력 형태
    ## cnv: chr, start, end, log2 혹은 cn
    ## sv : chr1, pos1, chr2, pos2, type

    cnv_data <- read.delim("/data/wgs/RUO_Service/WGS-26-02/WGS-26-02_02/hmftools/linx/WGS-26-02_02.linx.vis_copy_number.tsv")
    cnv_data <- cnv_data %>% mutate(
        CopyNumber = floor(CopyNumber+0.5),
        cnv_color = case_when(
            CopyNumber >= 9 ~ "#E41A1C",
            CopyNumber >= 5 & CopyNumber < 9 ~ "#dba9a7",
            CopyNumber <= 1 ~ "#6c8ebf",
            TRUE ~ "#eaeaea"
        )
    )

    cnv_gr <- makeGRangesFromDataFrame(
        cnv_data %>% mutate(Chromosome = paste0("chr", Chromosome)),
        seqnames.field = "Chromosome",
        start.field = "Start",
        end.field = "End",
        keep.extra.columns = TRUE
    )



    sv <- read.delim("/data/wgs/RUO_Service/WGS-26-02/WGS-26-02_02/hmftools/linx/WGS-26-02_02.linx.vis_sv_data.tsv") %>% 
        mutate( ChrStart = paste0("chr",ChrStart), ChrEnd = paste0("chr",ChrEnd) ) %>% 
        filter( Type %in% c("BND","INV","DEL") )

    sv_gr1 <- GRanges(sv$ChrStart, IRanges(start=sv$PosStart, width=1), sv_type=sv$Type)
    sv_gr2 <- GRanges(sv$ChrEnd,   IRanges(start=sv$PosEnd,   width=1), sv_type=sv$Type)

    genome <- getGenomeAndMask("hg19")$genome
    chrom_sizes <- seqlengths(genome)
    

    for( j in 1:nrow(sv)) {
        
        pdata_gr1 <- sv_gr1[j]
        pdata_gr2 <- sv_gr2[j]
        chrom_gr1 <- as.character(seqnames(pdata_gr1))
        chrom_gr2 <- as.character(seqnames(pdata_gr2))
        seqlevels(pdata_gr1) = seqlevels(pdata_gr2) = unique(c(chrom_gr1,chrom_gr2))
        pdata_cnv_gr  <- cnv_gr %>% dplyr::filter( seqnames %in% unique(c(chrom_gr1,chrom_gr2) ) )
        pdata_sv_type <- mcols(pdata_gr1)$sv_type

        message(sprintf(" [ Structural Variant Plot ] SV-TYPE = %s, CHROM = %s", pdata_sv_type, paste0(unique(c(chrom_gr1,chrom_gr2)), collapse=", ")))

        if( pdata_sv_type %in% c("INV","DEL") ){
            plot.params <- getDefaultPlotParams(plot.type=2)
        }else{
            plot.params <- getDefaultPlotParams(plot.type=3)
        }
        plot.params$ideogramheight <- 15
        plot.params$data2height    <- 100
        plot.params$data2inmargin  <- 0
        plot.params$ideogramlateralmargin <- 0.03
        
  
        png(file=sprintf("%s/%s.SV_%s.plot.png", res_dir, SeqID, j),  width=16, height=6, units="in", res=150 )

        if( pdata_sv_type %in% c("INV","DEL") ){
            kp <- plotKaryotype( genome="hg19", chromosomes = unique(c(chrom_gr1,chrom_gr2)), plot.type = 2, cex=2, plot.params = plot.params, labels.plotter = NULL )
            kpText( kp, chr = c(chrom_gr1), x = 0, y = 0.2, labels = chrom_gr1, pos=4, cex = 1.5, data.panel = 2)
        }else{
            kp <- plotKaryotype( genome="hg19", chromosomes = c(chrom_gr1,chrom_gr2), plot.type = 3, cex=2, plot.params = plot.params, labels.plotter = NULL )
            kpText( kp, chr = c(chrom_gr1), x = 0, y = 0.2, labels = chrom_gr1, pos=4, cex = 1.5, data.panel = 2)
            kpText( kp, chr = c(chrom_gr2), x = chrom_sizes[chrom_gr2], y = 0.2, labels = chrom_gr2, pos=2, cex = 1.5, data.panel = 2) # pos=4,
        }

        # cnv-panel
        kpRect(kp, data=pdata_cnv_gr, y0=0, y1=pdata_cnv_gr$CopyNumber, ymin=0, ymax=12, col=pdata_cnv_gr$cnv_color, border = NA)
        kpAxis( kp, data.panel = 1, ymin = 0, ymax = 12, r0 = 0, r1 = 1, numticks = 7, cex = 1.2 ) 
        kpAbline( kp, data.panel = 1, h = 0,  ymin = 0, ymax = 12)
        kpText( kp, chr = c(chrom_gr1), x = 0, y = 12, labels = "Copy Number", pos=4, cex = 1.5,  ymin = 0, ymax = 12 )

        # breakpoint marker
        kpPoints( kp, data = pdata_gr1, y = 0, data.panel = 2, r0 = 0.05, r1 = 0.1, col = "#1BA1E2", pch = 16, cex = 1.5 )
        kpPoints( kp, data = pdata_gr2, y = 0, data.panel = 2, r0 = 0.05, r1 = 0.1, col = "#1BA1E2", pch = 16, cex = 1.5 )
        kpPlotLinks( kp, data = pdata_gr1, data2 = pdata_gr2, data.panel = 2, r0 = 0.07, r1 = 0.55, col = "#1BA1E2", border = "#1BA1E2", lwd = 1.5)
        suppressWarnings(kpPlotMarkers( kp, data.panel = 2, data = gene_gr, labels = gene_gr$gene, r0 =0, r1 = 0.5, cex = 1))

        dev.off()
    }

    