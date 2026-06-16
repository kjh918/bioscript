
### packages usages 
    # tumor=runVarbin(
    #    dir,
    #    genome = c("hg38", "hg19"),
    #    resolution = c("220kb", "55kb", "110kb", "195kb", "280kb", "500kb", "1Mb", "2.8Mb"),
    #    remove_Y = FALSE, is_paired_end = FALSE, method = c("CBS", "multipcf"), vst = c("ft", "log"),
    #    seed = 17, min_bincount = 10, 
    #    alpha = 1e-05, # Significance level for changepoint detection in segmentation. Lower = fewer breakpoints; Higher = more breakpoints
    #    merge_levels_alpha = 1e-05, # Controls how aggressively similar CNV levels (e.g., 1.95 vs 2.05 copies) are merged. Larger values = less CN ; Smaller values more CN
    #    gamma = 40, # Regularization parameter controlling penalty for adding breakpoints. Larger : fewer segments, Lower : more segments
    #    name = "segment_ratios",
    #    BPPARAM = bpparam()
    # )
    # bincounts(tumor)
    # assay(tumor, 'ft')
    # ratios(tumor)
    # segment_ratios(tumor)
    # tumor <- runMetrics(tumor)
    # tumor <- findAneuploidCells(tumor)
    # tumor <- knnSmooth(tumor)
    # tumor <- calcInteger(tumor, method = 'scquantum', assay = 'smoothed_bincounts')
    # plotRatio(tumor)
###
#==================================================================================================#
### PACKAGES
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("copykit"))
    suppressPackageStartupMessages(library("miniUI"))
    suppressPackageStartupMessages(library("shiny"))
    suppressPackageStartupMessages(library("plyr"))
    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("parallel"))
    suppressPackageStartupMessages(library("RMySQL"))
    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("Cairo"))
    grDevices::X11.options(type='cairo')
    options(device='x11')
    options(bitmapType='cairo')
###
### ARGUMENTS
    option_list = list( 
        make_option(c("--SeqID"), action="store", default=NULL, type="character", help="sample-id. REQUIRED."),
        make_option(c("--AnalysisRunDir"), action="store", default=NULL, type="character", help="data folder. REQUIRED."),   
        make_option(c("--BinSize"), action="store", default=NULL, type="character", help="bin size. available choice : 220kb, 55kb, 110kb, 195kb, 280kb, 500kb, 1Mb, 2.8Mb. default = 110kb"),   
        make_option(c("--Ploidy"), action="store", default=NULL, type="character", help="sample ploidy. default = 2"),  
        make_option(c("--GenomeVersion"), action="store", default=NULL, type="character", help="reference genome version. default = hg38")
    )
    #--------------------------------------------------------------------------#  
    ARGS = parse_args(OptionParser(option_list=option_list))
    #--------------------------------------------------------------------------#
    SeqID          = ARGS$SeqID
    AnalysisRunDir = ARGS$AnalysisRunDir
    BinSize        = ARGS$BinSize
    Ploidy         = ARGS$Ploidy
    GenomeVersion  = ARGS$GenomeVersion

###
### MODIFIED FUNCTIONS
    plotRatio2 = function ( scCNA, PlotPloidy )
    {
        ploidy <- segment_ratio <- ratio <- start <- xstart <- xend <- NULL
        chr_ranges <- as.data.frame(SummarizedExperiment::rowRanges(scCNA))
        chr_lengths <- rle(as.numeric(chr_ranges$seqnames))$lengths
        chr_ranges_start <- chr_ranges %>% dplyr::group_by(seqnames) %>%
            dplyr::arrange(seqnames, start) %>% dplyr::filter(dplyr::row_number() == 1) %>% dplyr::ungroup()
        chr_ranges_end <- chr_ranges %>% dplyr::group_by(seqnames) %>%
            dplyr::arrange(seqnames, start) %>% dplyr::filter(dplyr::row_number() == dplyr::n()) %>% dplyr::ungroup()
        chrom_rects <- data.frame(chr = chr_ranges_start$seqnames, xstart = as.numeric(chr_ranges_start$abspos), xend = as.numeric(chr_ranges_end$abspos))
        xbreaks <- rowMeans(chrom_rects %>% dplyr::select(xstart, xend))
        if (nrow(chrom_rects)%%2 == 0) { 
            chrom_rects$colors <- c("#ffffff", "#e0e0e0")
        }else{
            chrom_rects$colors <- rep_len(c("#ffffff", "#e0e0e0"), nrow(chrom_rects))
        }
        ggchr_back <- list(geom_rect(data = chrom_rects, aes(xmin = xstart,xmax = xend, ymin = -Inf, ymax = Inf, fill = colors), alpha = 0.2), scale_fill_identity())
        sec_breaks <- c(0, 5e+08, 1e+09, 1.5e+09, 2e+09, 2.5e+09, 3e+09)
        sec_labels <- c(0, 0.5, 1, 1.5, 2, 2.5, 3)

        ggaes <- list(
            scale_x_continuous(
                breaks = xbreaks, labels = gsub("chr","", chrom_rects$chr), position = "top", expand = c(0,0), 
                sec.axis = sec_axis(~., breaks = sec_breaks, labels = sec_labels, name = "Continuous Genomic Position (Gb)")
            ), 
            theme_classic(), 
            xlab(""),
            ylab("CopyRatio"), 
            theme(
                axis.text.x = element_text(angle = 0, vjust = 0.5, size = 15), 
                axis.text.y = element_text(size = 15),
                legend.position = "none", 
                axis.ticks.x = element_blank(),
                axis.title = element_text(size = 15), 
                plot.title = element_text(size = 15),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
            )
        )
               
        abspos <- chr_ranges$abspos
        dat_seg <- copykit::segment_ratios(scCNA) %>% dplyr::mutate(abspos = abspos)
        dat_seg_l <- tidyr::gather(data = dat_seg, key = "sample", value = "segment_ratio", -abspos)
        dat_ratios <- copykit::ratios(scCNA) %>% dplyr::mutate(abspos = abspos)
        dat_ratios_l <- tidyr::gather(data = dat_ratios, key = "sample", value = "ratio", -abspos)
        
        if (nrow(dat_ratios_l) == nrow(dat_seg_l)) {
            df <- dat_ratios_l %>% dplyr::mutate(segment_ratio = dat_seg_l$segment_ratio)
        }else{
            stop("Nrow in copykit::segment_ratios() assay different than nrow in ratios().")
        }

        if (!is.null(colData(scCNA)$ploidy)) {
            df <- df %>% dplyr::mutate(
                integer = round(df$segment_ratio * (colData(scCNA)$ploidy[match(df$sample, colData(scCNA)$sample)]))
            )
        }
        # if( is.null(sample_name) )
        # {
        #     if( length(unique(df$sample)) != 1 )
        #     {
        #         stop("Multiple Samples Analyzed. Please Specify sample.")
        #     }else{
        #         df_plot = df
        #     }
        # }else{
        #     df_plot = df %>% filter( sample == sample_name )
        # }

        df_plot <- df 

        if (!is.null(colData(scCNA)$ploidy)) {
            cell_ploidy <- as.data.frame(colData(scCNA)) %>% dplyr::select(sample, ploidy) %>% pull(ploidy)
            color_ratio <- structure(ocean.balance(length(0:(2 * round(cell_ploidy)))), names = 0:(2 * round(cell_ploidy)))
            color_ratio[which(names(color_ratio) == round(cell_ploidy))] <- "gray"
            max_int_value <- round(max(df_plot$ratio) * cell_ploidy)
            if( Ploidy == 2 ){ if( max_int_value < 4 ){ max_int_value = 4 } }
            mean_bin_cell <- mean(df_plot$ratio)
            sec_axis_int <- list(
                scale_y_continuous(sec.axis = sec_axis(~. *cell_ploidy/mean_bin_cell, breaks = seq(0, max_int_value, 1), name = "Integer CopyNumber")))
            df_plot$integer[df_plot$integer > 2 * (round(mean(df_plot$integer)))] <- 2 * (round(mean(df_plot$integer)))
        }

        df_plot = df_plot %>% mutate( cn = integer * (mean_bin_cell/cell_ploidy))

        SampleName   = unique(df_plot$sample)
        plotMaxRatio = mean_bin_cell/PlotPloidy * round(PlotPloidy,0)*2.3
        #plotMaxRatio = df_plot %>% filter( integer <= Ploidy*2/2 ) %>% .$ratio %>% max()
        
        if( PlotPloidy == 2  ){ if( plotMaxRatio < 2.3 ){ plotMaxRatio =  2.3 } }
        #if( Ploidy  > 2.5){ if( mean_bin_cell*3 > plotMaxRatio ){ plotMaxRatio =  mean_bin_cell*3.1 } }
        
        p_cnv = ggplot(df_plot) + ggchr_back + ggaes + sec_axis_int + ggtitle( SampleName ) +
            geom_point(aes(abspos,ratio, color = as.character(integer)), pch = 19, size = 1.1, alpha = 0.6) + 
            geom_point(aes(abspos, cn), pch=15, col = "#2d2d2d", size = 1, alpha=0.3) + 
            scale_color_manual(values = color_ratio) +
            coord_cartesian(ylim=c(0,plotMaxRatio) )  

        return(p_cnv)
    }
    ##
    calcInteger2 = function (scCNA, assay = c("bincounts", "smoothed_bincounts",
        "segment_ratios"), method = "fixed", ploidy_value = NULL,
        name = "integer", penalty = 25, BPPARAM = bpparam())
    {
        assay = match.arg(assay)
        if ("smoothed_bincounts" %in% assayNames(scCNA) && assay ==
            "bincounts" && method == "scquantum") {
            warning("CopyKit detected that knnSmooth() has been performed.")
            warning("If working with knnSmooth datasets we recommend using the assay 'smoothed_bincounts'")
        }
        if (assay == "bincounts") {
            bin <- SummarizedExperiment::assay(scCNA, "bincounts")
        }
        if (assay %in% c("smoothed_bincounts", "segment_ratios")) {
            bin <- SummarizedExperiment::assay(scCNA, "segment_ratios")
        }
        seg <- SummarizedExperiment::assay(scCNA, "segment_ratios")
        if (!is.null(ploidy_value)) {
            if (method == "fixed") {
                if (is.null(ploidy_value) && !is.numeric(ploidy_value)) {
                    stop("Method fixed requires a numeric value for ploidy_value.")
                }
                message("Scaling ratio values by ploidy value ",
                    ploidy_value)
                SummarizedExperiment::colData(scCNA)$ploidy <- ploidy_value
                S4Vectors::metadata(scCNA)$ploidy_method <- "fixed"
            }
        }
        if (method == "metadata") {
            if (!is.null(colData(scCNA)$ploidy)) {
                message("Calculating integer values based on colData(scCNA)$ploidy info.")
            }
            else {
                stop("Method 'metadata' requires colData(scCNA)$ploidy information.")
            }
        }
        if (method == "scquantum") {
            rg <- as.data.frame(SummarizedExperiment::rowRanges(scCNA))
            if (assay %in% c("bincounts", "smoothed_bincounts")) {
                sc_quants <- BiocParallel::bplapply(assay(scCNA,
                    assay), scquantum::ploidy.inference, chrom = rg$seqnames,
                    start = rg$start, end = rg$end, penalty = penalty,
                    BPPARAM = BPPARAM)
            }
            if (assay == "segment_ratios") {
                sc_quants <- BiocParallel::bplapply(seq_along(seg),
                    function(z) {
                    segnums <- cumsum(c(TRUE, abs(diff(seg[, z])) >
                        1e-05))
                    seg_length <- rle(seg[, z])$lengths
                    seg_bins_mean <- tapply(bin[, z], segnums,
                        mean)
                    if (any(seg_length <= 3)) {
                        iod.est <- scquantum::timeseries.iod(bin[,
                        z])
                    }
                    else {
                        iod.est <- tapply(bin[, z], segnums, scquantum::timeseries.iod)
                    }
                    mean.est <- mean(bin[, z])
                    estimates <- scquantum::ploidy.inference(x = seg_bins_mean,
                        chrom = NULL, start = NULL, end = NULL, seg_length = seg_length,
                        iod = iod.est, mean_bincount = mean.est,
                        do_segmentation = FALSE)
                    })
            }
            sc_ploidies <- vapply(sc_quants, function(x) x$ploidy,
                numeric(1))
            sc_confidence <- vapply(sc_quants, function(x) x$confidence_ratio,
                numeric(1))
            ploidy_score <- abs(1 - sc_confidence)
            SummarizedExperiment::colData(scCNA)$ploidy <- sc_ploidies
            SummarizedExperiment::colData(scCNA)$confidence_ratio <- sc_confidence
            SummarizedExperiment::colData(scCNA)$ploidy_score <- ploidy_score
        }
        if (!identical(names(bin), colData(scCNA)$sample)) {
            stop("Order of cells in segment_ratios and colData() is not identical.")
        }
        if( dim(as.matrix(seg))[2] == 1 )
        {
            int_values <- data.frame(round(as.matrix(seg) %*% colData(scCNA)$ploidy))
            colnames(int_values) <- names(seg)
        }else{
            int_values <- round(as.matrix(seg) %*% diag(colData(scCNA)$ploidy)) %>%
            as.data.frame()
        }
        names(int_values) <- names(seg)
        SummarizedExperiment::assay(scCNA, name) <- int_values
        return(scCNA)
    }
    ##
    Call.CNV.from.CopyKit = function( scCNA )
    {
        chr_ranges <- as.data.frame(SummarizedExperiment::rowRanges(scCNA))
        chr_lengths <- rle(as.numeric(chr_ranges$seqnames))$lengths
        chr_ranges_start <- chr_ranges %>% dplyr::group_by(seqnames) %>%
            dplyr::arrange(seqnames, start) %>% dplyr::filter(dplyr::row_number() == 1) %>% dplyr::ungroup()
        chr_ranges_end <- chr_ranges %>% dplyr::group_by(seqnames) %>%
            dplyr::arrange(seqnames, start) %>% dplyr::filter(dplyr::row_number() == dplyr::n()) %>% dplyr::ungroup()
        abspos <- chr_ranges$abspos
        dat_seg <- copykit::segment_ratios(scCNA) %>% dplyr::mutate(abspos = abspos)
        dat_seg_l <- tidyr::gather(data = dat_seg, key = "sample", value = "segment_ratio", -abspos)
        dat_ratios <- copykit::ratios(scCNA) %>% dplyr::mutate(abspos = abspos)
        dat_ratios_l <- tidyr::gather(data = dat_ratios, key = "sample", value = "ratio", -abspos)

        seg_data <- dat_ratios_l %>% dplyr::mutate(segment_ratio = dat_seg_l$segment_ratio) %>% 
            dplyr::mutate( integer = round(segment_ratio * (colData(scCNA)$ploidy)) )
        
        cnv_raw = chr_ranges %>% left_join(., seg_data, by="abspos") %>% dplyr::rename(
            Chrom=seqnames, Start=start, End=end, Width=width, GC_content=gc_content,
            AbsGenomicPos=abspos, ChrArm=arm, SeqID=sample, Bin_logR=ratio, Segment_logR=segment_ratio, CN=integer
        ) %>% dplyr::select(!c("strand"))

        cn_values = cnv_raw$CN
        names(cn_values) = 1:nrow(cnv_raw)
        cn_groups = split(cn_values, c(0, cumsum(diff(cn_values) != 0)))

        cnv_result = mclapply( cn_groups, function(w) 
        {
            cn_group_data = cnv_raw[as.numeric(names(w)), ] 
            cn_group_logR = cn_group_data$Segment_logR
            names(cn_group_logR) = 1:nrow(cn_group_data)
            seg_groups = split(cn_group_logR, c(0, cumsum(diff(cn_group_logR) != 0)))
            
            seg_group_res = lapply( seg_groups, function(z) 
            {
                cn_group_data[as.numeric(names(z)), ] %>% group_by(Chrom, Segment_logR) %>% mutate(
                    Seg_Start=min(Start), Seg_End=max(End), Segment_CN=unique(CN)
                ) %>% mutate( Seg_Length = Seg_End-Seg_Start+1 ) %>% mutate( Seg_Length_Mb = round(Seg_Length/1000000,3) ) %>% dplyr::select(
                    c("Chrom","Seg_Start","Seg_End","Seg_Length","Seg_Length_Mb","Segment_logR","Segment_CN")
                ) %>% unique() %>% as.data.frame()
            }) %>% ldply() %>% dplyr::select(!c(".id"))
            
            return(seg_group_res)

        }, mc.cores=5) %>% ldply() %>% dplyr::select(!c(".id"))

        OutputResult = list(
            cnv_result = cnv_result,
            cnv_raw    = cnv_raw
        )
        return(OutputResult)
    }
    ##
###   
### MANUAL PARAMS 
    # SeqID   = "cbnipt.24.04.20_GM09888.cnv.ref.wgs"
    # BinSize = "110kb" # c("220kb", "55kb", "110kb", "195kb", "280kb", "500kb", "1Mb", "2.8Mb")
    # Ploidy  = 3.6
    # AnalysisRunDir = "/data/cbNIPT/copykit_analysis/cbnipt.24.03.11_SKBR3.10cells.ds.wga" 
    # GenomeVersion  = "hg38"
###
### PARAM CHECK
    if(is.null(SeqID) ){ stop(">> No SeqID. REQUIRED.") } 
    if(is.null(AnalysisRunDir) ){ stop(">> No Analysis-Run-Dir. REQUIRED.") }else{
        if(!dir.exists(AnalysisRunDir) ){ stop(">> Analysis-Run-Dir not found. No such folder.") }
    } 
    if(is.null(GenomeVersion)){ GenomeVersion = "hg38" }
    if(is.null(BinSize)      ){ BinSize = "110kb" }
    if(is.null(Ploidy)       ){ Ploidy = 2 }else{ Ploidy = as.numeric(Ploidy) }
###
### RUN CopyKit ANALYSIS
    # Read Data and Run Analysis
    CokyKitData = runVarbin(
        dir                = AnalysisRunDir,
        genome             = GenomeVersion,
        resolution         = BinSize,
        remove_Y           = FALSE, 
        min_bincount       = 1,
        alpha              = 0.01,
        merge_levels_alpha = 0.1,
        gamma              = 70,
        is_paired_end      = TRUE
    )
    CokyKitData = runMetrics(CokyKitData)
    # Add Sample Ploidy
    CokyKitData = calcInteger2(CokyKitData, method='fixed', ploidy_value = Ploidy )
    
    # CNV result as Plot
    p.cnv = plotRatio2(CokyKitData, PlotPloidy=Ploidy)
    ggsave(p.cnv, file=sprintf("%s/%s_binsize.%s_copykit.plot.CNA.png", AnalysisRunDir, SeqID, BinSize), width=24, height=6, unit='in', dpi=150 )

    # CNV result as Table
    cnv_result_list = Call.CNV.from.CopyKit(CokyKitData)
    write.table( 
        cnv_result_list$cnv_result, 
        sprintf("%s/%s_binsize.%s_copykit.called.CNA.tsv", AnalysisRunDir, SeqID, BinSize), 
        quote=F, col.names=T, row.names=F, sep="\t"
    )
    write.table(
        cnv_result_list$cnv_raw, 
        sprintf("%s/%s_binsize.%s_copykit.segmentation.raw.data.tsv", AnalysisRunDir, SeqID, BinSize), 
        quote=F, col.names=T, row.names=F, sep="\t"
    )

    copykit_stats = as.data.frame(colData(CokyKitData))
    write.table(
        copykit_stats,
        sprintf("%s/%s_binsize.%s_copykit.analysis.stats.tsv", AnalysisRunDir, SeqID, BinSize), 
        quote=F, col.names=T, row.names=F, sep="\t"
    )

    dbCon = dbConnect(dbDriver("MySQL"), host='192.168.0.34', user='gcx', port=3306, password='gencurix!!', db='gcx_cbnipt' )
    preclear  = dbGetQuery(dbCon, sprintf("DELETE FROM qc_copykit WHERE sample = '%s'", copykit_stats$sample))
    datawrite = dbWriteTable(dbCon, name="qc_copykit", value=copykit_stats, row.names=F, append=T)
    dbDisconnect(dbCon)

    # analysis object save
    saveRDS(CokyKitData, file=sprintf("%s/%s_binsize.%s_copykit.analysis.object.RDS", AnalysisRunDir, SeqID, BinSize))

###








