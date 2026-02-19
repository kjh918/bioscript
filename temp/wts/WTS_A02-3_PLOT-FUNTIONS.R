
#---| PACKAGES |-------------------------------------------------------------------------------------------------------------------------------------#
    options(stringsAsFactors=FALSE)   
    suppressPackageStartupMessages(library("Cairo"))
    grDevices::X11.options(type='cairo')
    options(device='x11')
    options(bitmapType='cairo')
#----------------------------------------------------------------------------------------------------------------------------------------------------#

#===| COLORS |=======================================================================================================================================#
    PredefineSampleGroupColors <- c(
        "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999","#66C2A5","#8DA0CB",
        "#A6D854","#FFD92F","#E5C494","#80B1D3","#FDB462","#B3DE69","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F","#FCCDE5"
    )
#====================================================================================================================================================#

#===| SAMPLE CLUSTERING PLOTS |======================================================================================================================#
#---/ UMAP result scatter plot (not ggplot2) /-----------------------------------------------------#
#' @description draw scatter-plot of umap clustering result 
#' @param UmapResult        result table of umap analysis
#' @param PlotDimension     2D-plot or 3D-plot. default = 3 (3D)
#' @param PlotName          plot filename. if NULL, default name will be used.
#' @param PlotTitle         plot title. if NULL, no title on plot.
#' @param Width             plot width. unit = 'inch'. 
#' @param Height            plot height unit = 'inch'. 
#' @param PlotSave          save plot as file or not. default = TRUE
#' @param DisplayOnScreen   display plot on screen with X11. default = FALSE
#' @param AxisColor         axis-line color for 2D-plot.
#' @param GridColor         axis & grid color of 3D-plot. box-border color in 2D-plot.
#' @param AxisLabelColor    axis-label font color. if NULL axis-color will be used as axis-label font color.
#' @param AxisFontSize      axis-label font size. default = 1.5
#' @param AxisTitleFontSize axis-title font size. default = 2
#' @param PointSize         scatter-plot point size. default = 3
#' @param PointColor        scatter-plot point colors. if NULL or less than the number of group, pre-defined colors will be used.
#' @param PointPch          scatter-plot point shape. either '19' or '21'. default = 19
#' @param Plot3D.Type       3D-plot type. either 'h' or 'p'. default = 'h'
#' @param Plot3D.Angle      3D-plot y-axis angle. default = 40
#' @export 
    WTS_plot.UMAP.Clustering.Scatter <- function( UmapResult, 
        PlotDimension=3, PlotName=NULL, PlotTitle=NULL, Width=8, Height=8, PlotSave=TRUE, DisplayOnScreen=FALSE, 
        AxisColor="#2d2d2d", GridColor="#d6d6d6", AxisLabelColor=NULL, AxisFontSize=1.5, AxisTitleFontSize=2,
        PointSize=1, PointColor=NULL, PointPch=19, Plot3D.Type="h", Plot3D.Angle=40 
    )
    {
        UMAP_COLORS <- c(
            "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999","#66C2A5","#8DA0CB",
            "#A6D854","#FFD92F","#E5C494","#80B1D3","#FDB462","#B3DE69","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F","#FCCDE5"
        )
        #-----------------------------------------------------------------------
        if( all(is.na(UmapResult$Z)) ){ CompN = 2 }else{ CompN = 3 }
        #-----------------------------------------------------------------------
        if( CompN < PlotDimension )
        {
            message("|---> Analysis Dimension is less than Plot-Dimension. Analysis Dimension will be used as Plot-Dimension.")
            message(sprintf("|---> Analysis Dimension = %s", CompN))
            PlotDimension <- CompN
        }
        #-----------------------------------------------------------------------
        if( any(is.null(UmapResult$group) | is.na(UmapResult$group)) ){ UmapResult$group <- UmapResult$seq_id }
        #-----------------------------------------------------------------------
        if( is.null(PointColor) )
        {
            GroupColors <- UMAP_COLORS
        }else{
            if( length(PointColor) < length(unique(UmapResult$group)) )
            {
                message("|---> provided colors are not enough to show groups. pre-defined default colors will be used.")
                GroupColors <- UMAP_COLORS
            }else{
                GroupColors <- PointColor[1:length(unique(UmapResult$group))]
            }
        }
        names(GroupColors) <- unique(UmapResult$group)
        UmapResult$colors  <- GroupColors[ UmapResult$group ]
        # plot filename --------------------------------------------------------
        if( is.null(PlotName) )
        { 
            message("|---> no plot filename. default filename will be used.")
            if( PlotDimension == 2 ){ PlotName <- "./UMAP.clustering.2D.scatterplot.png" }else{ PlotName <- "./UMAP.clustering.3D.scatterplot.png" }
        }
        # plot title and margin-------------------------------------------------
        TopMargin <- 5
        if( is.null(PlotTitle) )
        {
            message("|---> no plot title. ")
            TopMargin <- 2
            PlotTitle <- NULL
        }
        # axis label font color ------------------------------------------------
        if( is.null(AxisLabelColor) ){ AxisLabelColor <- AxisColor }
        # draw 2D plots --------------------------------------------------------
        if( PlotDimension == 2 )
        {
            if( CompN == 2 )
            {          
                # 2D SCATTER-PLOT AND SAVE AS FILE -----------------------------
                if( PlotSave )
                {
                    png( PlotName, width=Width, height=Height, units="in", res=150 )
                    par(mar=c(5,5,TopMargin,2))
                    plot( x=UmapResult$X, y=UmapResult$Y, xlab="D-1", ylab="D-2", main=PlotTitle, axes=FALSE, col=UmapResult$colors, bty="n", pch=PointPch, cex=PointSize, cex.lab=AxisTitleFontSize )
                    box( col=GridColor )
                    axis(1, col=AxisColor, col.ticks=AxisColor, col.axis=AxisLabelColor, cex.axis=AxisFontSize )
                    axis(2, col=AxisColor, col.ticks=AxisColor, col.axis=AxisLabelColor, cex.axis=AxisFontSize )
                    dev.off()
                }
                # Display on Screen --------------------------------------------
                if( DisplayOnScreen )
                {
                    x11( width=Width, height=Height, bg="white", family="Helvetica")
                    par(mar=c(5,5,TopMargin,2))
                    plot( x=UmapResult$X, y=UmapResult$Y, xlab="D-1", ylab="D-2", main=PlotTitle, axes=FALSE, col=UmapResult$colors, bty="n", pch=PointPch, cex=PointSize, cex.lab=AxisTitleFontSize )
                    box( col=GridColor )
                    axis(1, col=AxisColor, col.ticks=AxisColor, col.axis=AxisLabelColor, cex.axis=AxisFontSize )
                    axis(2, col=AxisColor, col.ticks=AxisColor, col.axis=AxisLabelColor, cex.axis=AxisFontSize )
                }
                #---------------------------------------------------------------
            }else{
                # 2D SCATTER-PLOT GROUPS AND SAVE AS FILE ----------------------
                if( PlotSave )
                {
                    png( PlotName, width=Width*3, height=Height, units="in", res=150 )
                    par(mfrow=c(1,3))
                    par(mar=c(5,5,TopMargin,2))
                    plot( x=UmapResult$X, y=UmapResult$Y, xlab="D-1", ylab="D-2", main=PlotTitle, axes=FALSE, col=UmapResult$colors, bty="n", pch=PointPch, cex=PointSize, cex.lab=AxisTitleFontSize )
                    box( col=GridColor )
                    axis(1, col=AxisColor, col.ticks=AxisColor, col.axis=AxisLabelColor, cex.axis=AxisFontSize )
                    axis(2, col=AxisColor, col.ticks=AxisColor, col.axis=AxisLabelColor, cex.axis=AxisFontSize )
                    plot( x=UmapResult$X, y=UmapResult$Z, xlab="D-1", ylab="D-3", main=PlotTitle, axes=FALSE, col=UmapResult$colors, bty="n", pch=PointPch, cex=PointSize, cex.lab=AxisTitleFontSize )
                    box( col=GridColor )
                    axis(1, col=AxisColor, col.ticks=AxisColor, col.axis=AxisLabelColor, cex.axis=AxisFontSize )
                    axis(2, col=AxisColor, col.ticks=AxisColor, col.axis=AxisLabelColor, cex.axis=AxisFontSize )
                    plot( x=UmapResult$Y, y=UmapResult$Z, xlab="D-2", ylab="D-3", main=PlotTitle, axes=FALSE, col=UmapResult$colors, bty="n", pch=PointPch, cex=PointSize, cex.lab=AxisTitleFontSize )
                    box( col=GridColor )
                    axis(1, col=AxisColor, col.ticks=AxisColor, col.axis=AxisLabelColor, cex.axis=AxisFontSize )
                    axis(2, col=AxisColor, col.ticks=AxisColor, col.axis=AxisLabelColor, cex.axis=AxisFontSize )
                    dev.off()
                }
                # Display on Screen --------------------------------------------
                if( DisplayOnScreen )
                {
                    x11( width=Width*3, height=Height, bg="white", family="Helvetica")
                    par(mfrow=c(1,3))
                    par(mar=c(5,5,TopMargin,2))
                    plot( x=UmapResult$X, y=UmapResult$Y, xlab="D-1", ylab="D-2", main=PlotTitle, axes=FALSE, col=UmapResult$colors, bty="n", pch=PointPch, cex=PointSize, cex.lab=AxisTitleFontSize )
                    box( col=GridColor )
                    axis(1, col=AxisColor, col.ticks=AxisColor, col.axis=AxisLabelColor, cex.axis=AxisFontSize )
                    axis(2, col=AxisColor, col.ticks=AxisColor, col.axis=AxisLabelColor, cex.axis=AxisFontSize )
                    plot( x=UmapResult$X, y=UmapResult$Z, xlab="D-1", ylab="D-3", main=PlotTitle, axes=FALSE, col=UmapResult$colors, bty="n", pch=PointPch, cex=PointSize, cex.lab=AxisTitleFontSize )
                    box( col=GridColor )
                    axis(1, col=AxisColor, col.ticks=AxisColor, col.axis=AxisLabelColor, cex.axis=AxisFontSize )
                    axis(2, col=AxisColor, col.ticks=AxisColor, col.axis=AxisLabelColor, cex.axis=AxisFontSize )
                    plot( x=UmapResult$Y, y=UmapResult$Z, xlab="D-2", ylab="D-3", main=PlotTitle, axes=FALSE, col=UmapResult$colors, bty="n", pch=PointPch, cex=PointSize, cex.lab=AxisTitleFontSize )
                    box( col=GridColor )
                    axis(1, col=AxisColor, col.ticks=AxisColor, col.axis=AxisLabelColor, cex.axis=AxisFontSize )
                    axis(2, col=AxisColor, col.ticks=AxisColor, col.axis=AxisLabelColor, cex.axis=AxisFontSize )
                }
            }
        }
        # draw 3D plots --------------------------------------------------------
        if( PlotDimension == 3 )
        {
            suppressPackageStartupMessages(library("scatterplot3d"))
            # 3D SCATTER-PLOT AND SAVE AS FILE ---------------------------------
            if( PlotSave )
            {
                png( PlotName, width=Width, height=Height, units="in", res=150 )
                if( PointPch == 21 )
                {
                    scatterplot3d( x=UmapResult$X, y=UmapResult$Y, z=UmapResult$Z, xlab="D-1", ylab="D-2", zlab="D-3", label.tick.marks=FALSE,
                        pch=PointPch, cex.symbols=PointSize, col.axis=GridColor, col.grid=GridColor, y.axis.offset=0.7,
                        type=Plot3D.Type, angle=Plot3D.Angle, lty.axis=1, lty.grid=2, lty.hplot=3, color=GridColor, bg=UmapResult$colors, lwd=2
                    )    
                }else{
                    scatterplot3d( x=UmapResult$X, y=UmapResult$Y, z=UmapResult$Z, xlab="D-1", ylab="D-2", zlab="D-3", label.tick.marks=FALSE,
                        pch=PointPch, cex.symbols=PointSize, col.axis=GridColor, col.grid=GridColor, y.axis.offset=0.7,
                        type=Plot3D.Type, angle=Plot3D.Angle, lty.axis=1, lty.grid=2, lty.hplot=3, color=UmapResult$colors
                    )
                }
                dev.off()
            }
            # Display on Screen ------------------------------------------------
            if( DisplayOnScreen )
            {
                x11( width=Width, height=Height, bg="white", family="Helvetica")
                if( PointPch == 21 )
                {
                    scatterplot3d( x=UmapResult$X, y=UmapResult$Y, z=UmapResult$Z, xlab="D-1", ylab="D-2", zlab="D-3", label.tick.marks=FALSE,
                        pch=PointPch, cex.symbols=PointSize, col.axis=GridColor, col.grid=GridColor, y.axis.offset=0.7,
                        type=Plot3D.Type, angle=Plot3D.Angle, lty.axis=1, lty.grid=2, lty.hplot=3, color=GridColor, bg=UmapResult$colors, lwd=2
                    )    
                }else{
                    scatterplot3d( x=UmapResult$X, y=UmapResult$Y, z=UmapResult$Z, xlab="D-1", ylab="D-2", zlab="D-3", label.tick.marks=FALSE,
                        pch=PointPch, cex.symbols=PointSize, col.axis=GridColor, col.grid=GridColor, y.axis.offset=0.7,
                        type=Plot3D.Type, angle=Plot3D.Angle, lty.axis=1, lty.grid=2, lty.hplot=3, color=UmapResult$colors
                    )
                }            
            }
        }
        #-----------------------------------------------------------------------
        return(UmapResult)
    }
#--------------------------------------------------------------------------------------------------#
#---/ UMAP result scatter plot (use ggplot2) /-----------------------------------------------------#
#' @description draw scatter-plot of umap clustering result 
#' @param UmapResult        result table of umap analysis
#' @param DisplaySampleGroup sample group ids for display on plot. if NULL, use all groups. default = NULL
#' @param PlotDimension     2D-plot or 3D-plot. default = 3 (3D) 
#' @param PlotName          plot filename. if NULL, default name will be used.
#' @param PlotTitle         plot title. if NULL, no title on plot.
#' @param Width             plot width. unit = 'inch'. 
#' @param Height            plot height unit = 'inch'. 
#' @param PlotSave          save plot as file or not. default = TRUE
#' @param DisplayOnScreen   display plot on screen with X11. default = FALSE
#' @param AxisColor         axis-line color for 2D-plot.
#' @param GridColor         axis & grid color of 3D-plot. box-border color in 2D-plot.
#' @param AxisLabelColor    axis-label font color. if NULL axis-color will be used as axis-label font color.
#' @param AxisFontSize      axis-label font size. default = 1.5
#' @param AxisTitleFontSize axis-title font size. default = 2
#' @param PointSize         scatter-plot point size. default = 3
#' @param PointColor        scatter-plot point colors. if NULL or less than the number of group, pre-defined colors will be used.
#' @param PointPch          scatter-plot point shape. either '19' or '21'. default = 19
#' @export 
    WTS_plot.UMAP.Clustering.Scatter2 <- function( UmapResult, DisplaySampleGroup=NULL,
        PlotDimension=3, PlotName=NULL, PlotTitle=NULL, Width=8, Height=8, PlotSave=TRUE, DisplayOnScreen=FALSE, 
        AxisColor="#2d2d2d", GridColor="#d6d6d6", AxisLabelColor=NULL, AxisFontSize=10, AxisTitleFontSize=12, PointSize=7, PointColor=NULL, PointPch=21
    )
    {
        suppressPackageStartupMessages(library("ggplot2"))
        suppressPackageStartupMessages(library("Cairo"))
        grDevices::X11.options(type='cairo')
        options(device='x11')
        options(bitmapType='cairo')
        #-----------------------------------------------------------------------
        UMAP_COLORS <- c(
            "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#80B1D3","#A65628","#F781BF","#999999","#66C2A5","#8DA0CB",
            "#A6D854","#FFD92F","#E5C494","#FFFF33","#FDB462","#B3DE69","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F","#FCCDE5"
        )
        #-----------------------------------------------------------------------
        if( all(is.na(UmapResult$Z)) ){ CompN = 2 }else{ CompN = 3 }
        #-----------------------------------------------------------------------
        if( CompN < PlotDimension )
        {
            message("|---> Analysis Dimension is less than Plot-Dimension. Analysis Dimension will be used as Plot-Dimension.")
            message(sprintf("|---> Analysis Dimension = %s", CompN))
            PlotDimension <- CompN
        }
        #-----------------------------------------------------------------------
        if( any(is.null(UmapResult$group) | is.na(UmapResult$group)) ){ UmapResult$group <- UmapResult$seq_id }
        #-----------------------------------------------------------------------
        if( is.null(PointColor) )
        {
            GroupColors <- UMAP_COLORS
        }else{
            if( length(PointColor) < length(unique(UmapResult$group)) )
            {
                message("|---> provided colors are not enough to show groups. pre-defined default colors will be used.")
                GroupColors <- UMAP_COLORS
            }else{
                GroupColors <- PointColor[1:length(unique(UmapResult$group))]
            }
        }
        names(GroupColors) <- unique(UmapResult$group)
        UmapResult$colors  <- GroupColors[ UmapResult$group ]
        # plot filename --------------------------------------------------------
        if( is.null(PlotName) )
        { 
            message("|---> no plot filename. default filename will be used.")
            if( PlotDimension == 2 ){ PlotName <- "./UMAP.clustering.2D.scatterplot.png" }else{ PlotName <- "./UMAP.clustering.3D.scatterplot.png" }
        }
        # plot title and margin-------------------------------------------------
        TopMargin <- 5
        if( is.null(PlotTitle) )
        {
            message("|---> no plot title. ")
            TopMargin <- 2
            PlotTitle <- NULL
        }
        # axis label font color ------------------------------------------------
        if( is.null(AxisLabelColor) ){ AxisLabelColor <- AxisColor }
        # group filtering 
        if( is.null(DisplaySampleGroup) )
        { PlotSampleGroup <- UmapResult$group }else{ PlotSampleGroup <- unlist(strsplit(DisplaySampleGroup, ","))
        }
        UmapPlotInput <- UmapResult %>% filter( group %in% PlotSampleGroup )
        # group colors ---------------------------------------------------------
        GROUP_COLORS        <- UmapPlotInput$colors      
        names(GROUP_COLORS) <- UmapPlotInput$group
        # draw 2D plots --------------------------------------------------------
        if( PlotDimension == 2 )
        {
            if( CompN == 2 )
            {   
                # draw default plot --------------------------------------------
                if( PointPch == 19 )
                {
                    UmapPlot <- ggplot( UmapPlotInput, aes( x=X, y=Y, color=group )) + geom_point( size=PointSize, pch=PointPch, alpha=0.9 ) +
                        labs( x="D-1", y="D-2", color="Sample Group") +
                        scale_color_manual(values=GROUP_COLORS)
                }else{
                    UmapPlot <- ggplot( UmapPlotInput, aes( x=X, y=Y, fill=group )) + geom_point( size=PointSize, pch=PointPch, color="white", alpha=0.9  ) +
                        labs( x="D-1", y="D-2", fill="Sample Group") +
                        scale_fill_manual(values=GROUP_COLORS)
                }
                # add plot title -----------------------------------------------
                if( !is.null(PlotTitle) ){ UmapPlot <- UmapPlot + ggtitle( PlotTitle ) }
                # add plot theme -----------------------------------------------
                UmapPlot <- UmapPlot + theme(
                    panel.background = element_rect(fill='white'), 
                    panel.grid.major = element_line(color = GridColor, linetype = 'dotted'),
                    axis.text        = element_text(family="Helvetica", colour=AxisLabelColor, size=AxisFontSize),
                    axis.title       = element_text(family="Helvetica", colour=AxisLabelColor, size=AxisTitleFontSize),
                    axis.line        = element_line(colour=AxisColor),
                    legend.title     = element_text(family="Helvetica", colour=AxisLabelColor, size=AxisTitleFontSize, face="bold"),
                    legend.text      = element_text(family="Helvetica", colour=AxisLabelColor, size=AxisTitleFontSize),
                    plot.margin      = unit(c(Height*0.05, Width*0.05, Height*0.05, Width*0.03), "in"),
                    plot.title       = element_text(family="Helvetica", colour=AxisLabelColor, size=AxisTitleFontSize*1.1),
                )
                # 2D SCATTER-PLOT AND SAVE AS FILE -----------------------------
                if( PlotSave )
                { suppressWarnings(ggsave(UmapPlot, file=PlotName, width=Width, height=Height, unit='in', dpi=150, type='cairo')) }
                # Display on Screen --------------------------------------------
                if( DisplayOnScreen )
                {
                    x11( width=Width, height=Height, bg="white", family="Helvetica")
                    print(UmapPlot)
                }
                #---------------------------------------------------------------
            }else{
                # draw default plot --------------------------------------------
                if( PointPch == 19 )
                {
                    UmapPlotList <- list(
                        D1_D2 <- ggplot( UmapPlotInput, aes( x=X, y=Y, color=group )) + geom_point( size=PointSize, pch=PointPch, alpha=0.9 ) +
                            labs( x="D-1", y="D-2", color="Sample Group") +
                            scale_color_manual(values=GROUP_COLORS),
                        D1_D3 <- ggplot( UmapPlotInput, aes( x=X, y=Z, color=group )) + geom_point( size=PointSize, pch=PointPch, alpha=0.9 ) +
                            labs( x="D-1", y="D-3", color="Sample Group") +
                            scale_color_manual(values=GROUP_COLORS),
                        D2_D3 <- ggplot( UmapPlotInput, aes( x=Y, y=Z, color=group )) + geom_point( size=PointSize, pch=PointPch, alpha=0.9 ) +
                            labs( x="D-2", y="D-3", color="Sample Group") +
                            scale_color_manual(values=GROUP_COLORS)
                    )
                }else{
                    UmapPlotList <- list(
                        D1_D2 <- ggplot( UmapPlotInput, aes( x=X, y=Y, fill=group )) + geom_point( size=PointSize, pch=PointPch, color="white", alpha=0.9  ) +
                            labs( x="D-1", y="D-2", fill="Sample Group") +
                            scale_fill_manual(values=GROUP_COLORS),
                        D1_D2 <- ggplot( UmapPlotInput, aes( x=X, y=Z, fill=group )) + geom_point( size=PointSize, pch=PointPch, color="white", alpha=0.9  ) +
                            labs( x="D-1", y="D-3", fill="Sample Group") +
                            scale_fill_manual(values=GROUP_COLORS),
                        D1_D2 <- ggplot( UmapPlotInput, aes( x=Y, y=Z, fill=group )) + geom_point( size=PointSize, pch=PointPch, color="white", alpha=0.9  ) +
                            labs( x="D-2", y="D-3", fill="Sample Group") +
                            scale_fill_manual(values=GROUP_COLORS)
                    )
                }
                # add plot title -----------------------------------------------
                if( !is.null(PlotTitle) ){ UmapPlotList <- lapply(UmapPlotList, function(pu) return( pu + ggtitle( PlotTitle ) ) ) }
                # add plot theme -----------------------------------------------
                UmapPlotList <- lapply( UmapPlotList, function(pu)
                {
                    pu <- pu + theme(
                        panel.background = element_rect(fill='white'), 
                        panel.grid.major = element_line(color = GridColor, linetype = 'dotted'),
                        axis.text        = element_text(family="Helvetica", colour=AxisLabelColor, size=AxisFontSize),
                        axis.title       = element_text(family="Helvetica", colour=AxisLabelColor, size=AxisTitleFontSize),
                        axis.line        = element_line(colour=AxisColor),
                        legend.title     = element_text(family="Helvetica", colour=AxisLabelColor, size=AxisTitleFontSize, face="bold"),
                        legend.text      = element_text(family="Helvetica", colour=AxisLabelColor, size=AxisTitleFontSize),
                        plot.margin      = unit(c(Height*0.05, Width*0.05, Height*0.05, Width*0.03), "in"),
                        plot.title       = element_text(family="Helvetica", colour=AxisLabelColor, size=AxisTitleFontSize*1.1),
                    )
                    return(pu)
                })
                UmapPlot <- suppressWarnings(cowplot::plot_grid( plotlist=UmapPlotList, ncol=3, rel_widths=c(1,1,1)))
                # 2D SCATTER-PLOT GROUPS AND SAVE AS FILE ----------------------
                if( PlotSave )
                { suppressWarnings(ggsave(UmapPlot, file=PlotName, width=Width*3, height=Height, unit='in', dpi=150, type='cairo')) }
                # Display on Screen --------------------------------------------
                if( DisplayOnScreen )
                {
                    x11( width=Width*3, height=Height, bg="white", family="Helvetica")
                    print(UmapPlot)
                }
            }
        }
        # draw 3D plots --------------------------------------------------------
        if( PlotDimension == 3 )
        {
            suppressPackageStartupMessages(library("gg3D"))
            suppressPackageStartupMessages(library("cowplot"))
            # 3D SCATTER-PLOT AND SAVE AS FILE ---------------------------------
            suppressWarnings(
                if( PointPch == 19 )
                {
                    preplot <- ggplot( UmapPlotInput, aes(x=X, y=Y, z=Z, color=group) ) + geom_point( pch=PointPch, size=PointSize, alpha=0.9 ) +
                        scale_color_manual(values=GROUP_COLORS) + labs( color="Sample Group" ) + 
                        theme( 
                            panel.background = element_rect(fill='white'), 
                            panel.grid.major = element_line(color = GridColor, linetype = 'dotted'),
                            legend.title     = element_text(family="Helvetica", colour=AxisLabelColor, size=AxisTitleFontSize, face="bold"),
                            legend.text      = element_text(family="Helvetica", colour=AxisLabelColor, size=AxisTitleFontSize)
                        )
                    UmapPlotLegend <- get_legend(preplot)
                    
                    UmapNoLgdPlot <- ggplot( UmapPlotInput, aes(x=X, y=Y, z=Z, color=group) ) +
                        scale_color_manual(values=GROUP_COLORS) + 
                        axes_3D( theta=135, phi=45 ) + 
                        stat_3D( theta=135, phi=45, size=PointSize, alpha=0.9, pch=PointPch ) 
                }else{
                    preplot <- ggplot( UmapPlotInput, aes(x=X, y=Y, z=Z, fill=group) ) + geom_point( pch=PointPch, size=PointSize, color="white", alpha=0.9 ) +
                        scale_fill_manual(values=GROUP_COLORS) + labs( fill="Sample Group" ) + 
                        theme( 
                            panel.background = element_rect(fill='white'), 
                            panel.grid.major = element_line(color = GridColor, linetype = 'dotted'),
                            legend.title     = element_text(family="Helvetica", colour=AxisLabelColor, size=AxisTitleFontSize, face="bold"),
                            legend.text      = element_text(family="Helvetica", colour=AxisLabelColor, size=AxisTitleFontSize)
                        )
                    UmapPlotLegend <- get_legend(preplot)

                    UmapNoLgdPlot <- ggplot( UmapPlotInput, aes(x=X, y=Y, z=Z, fill=group) ) +
                        scale_fill_manual(values=GROUP_COLORS) + 
                        axes_3D( theta=135, phi=45 ) + 
                        stat_3D( theta=135, phi=45, size=PointSize, alpha=0.9, pch=PointPch, color="white", geom="point" ) 
                }
            )
            # add plot title ---------------------------------------------------
            if( !is.null(PlotTitle) ){ UmapNoLgdPlot <- UmapNoLgdPlot + ggtitle( PlotTitle ) }
            # complete main plot -----------------------------------------------
            suppressWarnings(
                UmapNoLgdPlot <- UmapNoLgdPlot +
                    axes_3D( theta=135, phi=45, color=AxisColor ) + 
                    stat_3D( theta=135, phi=45, size=PointSize, alpha=0.9, pch=PointPch ) +
                    labs_3D( labs=c("D-1","D-2","D-3"), size=c(8,8,8) ) +
                    theme(
                        panel.background  = element_rect(fill='white'), 
                        panel.grid.major  = element_line(color = GridColor, linetype = 'dotted'),
                        axis.ticks.length = unit(0, "cm"),
                        axis.text         = element_blank(),
                        axis.title        = element_blank(),
                        legend.position   = "none",
                        plot.margin       = unit(c(Height*0.05, Width*0.05, Height*0.05, Width*0.03), "in"),
                        plot.title        = element_text(family="Helvetica", colour=AxisLabelColor, size=AxisTitleFontSize*1.1)
                    )
            )
            # add legend to main plot ------------------------------------------
            UmapPlot <- cowplot::plot_grid(UmapNoLgdPlot, UmapPlotLegend, ncol=2, rel_widths=c(0.8,0.2))
            # save plot --------------------------------------------------------
            if( PlotSave )
            { suppressWarnings(ggsave(UmapPlot, file=PlotName, width=Width, height=Height, unit='in', dpi=150, type='cairo')) }
            # Display on Screen --------------------------------------------
            if( DisplayOnScreen )
            {
                x11( width=Width*3, height=Height, bg="white", family="Helvetica")
                suppressWarnings(print(UmapPlot))
            } 
        } 
    }
#--------------------------------------------------------------------------------------------------#
#---/ PCA result biplot /--------------------------------------------------------------------------#
#' @description draw PCA result biplot
#' @param pca_sample_stats   REQUIRED. PCA analysis sample stats data. 
#' @param pca_feature_stats  PCA analysis feature stats data.
#' @param PlotType           plot type. either 'biplot' or 'sampleOnly'. default = 'biplot' ( = sample + feature ) 
#' @param PlotName           plot filename. if NULL, default name will be used.
#' @param PlotTitle          plot title. if NULL, no title on plot.
#' @param Width              plot width. unit = 'inch'. 
#' @param Height             plot height unit = 'inch'. 
#' @param PlotSave           save plot as file or not. default = TRUE
#' @param DisplayOnScreen    display plot on screen with X11. default = FALSE
#' @param AxisColor          axis-line color for 2D-plot.
#' @param GridColor          axis & grid color of 3D-plot. box-border color in 2D-plot.
#' @param AxisLabelColor     axis-label font color. if NULL axis-color will be used as axis-label font color.
#' @param AxisFontSize       axis-label font size. default = 1.5
#' @param AxisTitleFontSize  axis-title font size. default = 2
#' @param PointSize          scatter-plot point size. default = 3
#' @param PointPch           scatter-plot point shape. either '19' or '21'. default = 19
#' @param LabelsRepel        apply ggrepel to label or not. default = TRUE  
#' @export 
    WTS_plot.PCA.Clustering.Scatter <- function( pca_sample_stats=NULL, pca_feature_stats=NULL, PlotType="biplot",
        PlotName=NULL, PlotTitle=NULL, Width=8, Height=8, PlotSave=TRUE, DisplayOnScreen=FALSE, 
        AxisColor="#2d2d2d", GridColor="#d6d6d6", AxisLabelColor=NULL, AxisFontSize=10, AxisTitleFontSize=12, 
        PointSize=7, PointPch=21, LabelsRepel=TRUE
    )
    {
        suppressPackageStartupMessages(library("ggplot2"))
        suppressPackageStartupMessages(library("ggrepel"))
        suppressPackageStartupMessages(library("Cairo"))
        grDevices::X11.options(type='cairo')
        options(device='x11')
        options(bitmapType='cairo')
        # pre-defined colors ---------------------------------------------------
        PCA_COLORS <- c(
            "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#80B1D3","#A65628","#F781BF","#999999","#66C2A5","#8DA0CB",
            "#A6D854","#FFD92F","#E5C494","#FFFF33","#FDB462","#B3DE69","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F","#FCCDE5"
        )
        # plot type and data check ---------------------------------------------
        if( PlotType == "biplot" )
        {  
            if( any(is.null(pca_sample_stats) | nrow(pca_sample_stats) == 0 | is.null(pca_feature_stats) | nrow(pca_feature_stats) == 0 ) )
            {
                if( any(is.null(pca_sample_stats)  | nrow(pca_sample_stats) == 0)  ){ sample_stat_check  = 0 }else{ sample_stat_check  = 1 }
                if( any(is.null(pca_feature_stats) | nrow(pca_feature_stats) == 0) ){ feature_stat_check = 0 }else{ feature_stat_check = 1 }
                if( all(sample_stat_check == 0 & sample_stat_check == 0) )
                { stop("|---!!! No sample_stats and feature_stats are found. BOTH REQUIRED for biplot type. please check again. STOPPED.") 
                }else if( all(sample_stat_check == 1 & sample_stat_check == 0) )
                { stop("|---!!! No feature_stats is found. REQUIRED for biplot type. please check again. STOPPED.") 
                }else if( all(sample_stat_check == 0 & sample_stat_check == 1) ) 
                { stop("|---!!! No sample_stats is found. REQUIRED for biplot type. please check again. STOPPED.") }
            }
        }else if( PlotType == "sampleOnly" )
        {
            if( any(is.null(pca_sample_stats) | nrow(pca_sample_stats) == 0) )
            {  stop("|---!!! No sample_stats is found. REQUIRED for sampleOnly type. please check again. STOPPED.") }
        }
        # sample group colors --------------------------------------------------
        pca_sample_group_colors        <- PCA_COLORS[1:length(unique(pca_sample_stats$group))]
        names(pca_sample_group_colors) <- unique(pca_sample_stats$group)
        # plot x,y range -------------------------------------------------------
        plot_range_x <- range(c(range(pca_sample_stats$X), range(pca_feature_stats$X)))
        plot_range_y <- range(c(range(pca_sample_stats$Y), range(pca_feature_stats$Y)))
        # feature stats --------------------------------------------------------
        pca_feature_stats2     <- pca_feature_stats
        feature_scaling_factor <- min(c( as.integer(min(abs(plot_range_x)) / abs(max(pca_feature_stats2$X))) * 0.8, as.integer(min(abs(plot_range_y)) / abs(max(pca_feature_stats2$Y))) * 0.8 ))
        pca_feature_stats2$X   <- pca_feature_stats2$X * feature_scaling_factor
        pca_feature_stats2$Y   <- pca_feature_stats2$Y * feature_scaling_factor

        # draw main plots ------------------------------------------------------
        if( PlotType == "biplot" )
        {
            pca_plot <- ggplot() +
                geom_segment( data=pca_feature_stats2, aes(xend=X, yend=Y), x=0, y=0, color="#c1c1c1", alpha=1, linewidth=1, arrow=arrow( type="closed", angle=30, length=unit(0.2, "cm")) ) + 
                geom_text_repel( data=pca_feature_stats2, aes(x=X, y=Y, label=gene_symbol), point.padding = 0.5 ) +
                geom_vline( xintercept=0, color = AxisColor) +
                geom_hline( yintercept=0, color = AxisColor) +
                coord_cartesian( xlim = plot_range_x, ylim = plot_range_y ) +
                geom_point( data=pca_sample_stats, aes(x=X, y=Y, fill=group), size = PointSize, alpha=0.9, color="white", pch=PointPch )
        }else{
            pca_plot <- ggplot() +
                geom_vline( xintercept=0, color = AxisColor) +
                geom_hline( yintercept=0, color = AxisColor) +
                coord_cartesian( xlim = plot_range_x, ylim = plot_range_y ) +
                geom_point( data=pca_sample_stats, aes(x=X, y=Y, fill=group), size = PointSize, alpha=0.9, color="white", pch=PointPch )
        }
        # label repel ----------------------------------------------------------
        if( LabelsRepel )
        {
            if( length(which(colnames(pca_sample_stats) == "sample_name")) == 1 )
            {
                pca_plot <- pca_plot + 
                    geom_label_repel(data=pca_sample_stats, aes(x=X, y=Y, label=sample_name), box.padding=0.5, point.padding=0.5, color=pca_sample_stats$colors)
            }else{
                pca_plot <- pca_plot + 
                    geom_label_repel(data=pca_sample_stats, aes(x=X, y=Y, label=seq_id), box.padding=0.5, point.padding=0.5, color=pca_sample_stats$colors)             
            }
        }
        # plot title -----------------------------------------------------------
        if( !is.null( PlotTitle ) )
        {
            pca_plot <- pca_plot + ggtitle( PlotTitle )
        }
        # finalize plot -------------------------------------------------------- 
        pca_plot <- pca_plot +
            scale_fill_manual( values = pca_sample_group_colors ) +
            labs( x="PC1", y="PC2", fill="Sample Group" ) + 
            theme( 
                panel.background  = element_rect(fill='white', color=GridColor), 
                panel.grid.major  = element_line(color = GridColor, linetype = 'dotted'),
                axis.ticks.length = unit(0, "cm"),
                axis.text         = element_text(family="Helvetica", colour=AxisLabelColor, size=AxisFontSize),
                axis.title        = element_text(family="Helvetica", colour=AxisLabelColor, size=AxisTitleFontSize),
                plot.margin       = unit(c(Height*0.05, Width*0.05, Height*0.05, Width*0.03), "in"),
                legend.text       = element_text(family="Helvetica", colour=AxisLabelColor, size=AxisFontSize),
                legend.title      = element_text(family="Helvetica", colour=AxisLabelColor, size=AxisTitleFontSize),
                plot.title        = element_text(family="Helvetica", colour=AxisLabelColor, size=AxisTitleFontSize*1.1)
            )
        #-----------------------------------------------------------------------

        # save plot ------------------------------------------------------------
        if( PlotSave )
        { 
            # plot filename ----------------------------------------------------
            if( is.null(PlotName) )
            { 
                message("|---> no plot filename. default filename will be used.")
                if( PlotType == "biplot" )
                { PlotName <- "./PCA.clustering.result.biplot.png" }else{ PlotName <- "./PCA.clustering.result.smaple.only.scatter.plot.png" }
            }
            suppressWarnings(ggsave(pca_plot, file=PlotName, width=Width, height=Height, unit='in', dpi=150, type='cairo')) 
        }
        # display on screen ----------------------------------------------------
        if( DisplayOnScreen )
        {
            x11( width=Width*3, height=Height, bg="white", family="Helvetica")
            suppressWarnings(print(pca_plot))
        } 
        #-----------------------------------------------------------------------
        return(pca_plot)
    }
#--------------------------------------------------------------------------------------------------#
#---/ Relative sample distance heatmap /-----------------------------------------------------------#
#' @description draw sample distance heatmap
#' @param DistanceMatrix  REQUIRED. imput sample distance calculation result matrix
#' @param FontSize        label font size. default = 12
#' @param Clustering      run row/column clustgering. default = TRUE
#' @param CellSize        heatmap cell size. unit ="mm", default = 6
#' @export 
    WTS_plot.Sample.Distance.Heatmap <- function( DistanceMatrix=NULL, 
        FontSize=12, Clustering=TRUE, CellSize = 6
    )
    {
        suppressPackageStartupMessages(library("ComplexHeatmap"))
        suppressPackageStartupMessages(library("circlize"))
        #-----------------------------------------------------------------------
        DIST_COLORS <- c("#002c46", "#7fb6d7")
        #-----------------------------------------------------------------------
        if( is.null(DistanceMatrix) ){ stop("|---!!! No input distance matrix found. REQUIRED. please check again. STOPPED.") }
        #-----------------------------------------------------------------------
        HeatmapColorFun <- colorRamp2( 
            breaks       = c( 0, 0.05, 0.1, 0.15, 0.5 ), 
            colors       = c( "#002c46", "#004d7a", "#328bbf", "#7fb6d7", "#FFFFFF" ), 
            transparency = 0
        )
        #-----------------------------------------------------------------------
        DistHeatmap <- Heatmap( DistanceMatrix, col = HeatmapColorFun, na_col = "#ffffff", 
            # borders
            border_gp = gpar(col = "#eaeaea", lwd = 0.5), rect_gp = gpar(col = "#eaeaea", lwd = 0.25),
            # clustering
            cluster_rows = Clustering, cluster_columns = Clustering, show_row_dend = Clustering, show_column_dend = Clustering, 
            # row labels              
            row_names_side = "right", row_names_gp = gpar(fontsize = FontSize, col = "#2d2d2d" ), row_names_centered = FALSE,
            # column labels 
            column_names_side = "bottom", column_names_rot = 90, column_names_gp = gpar(fontsize = FontSize, col = "#2d2d2d"), column_names_centered = FALSE,  
            # legend
            heatmap_legend_param = list(
                title         = "Distance", 
                title_gp      = gpar(fontsize=FontSize*1.2),
                labels_gp     = gpar(fontsize=FontSize), 
                legend_width  = unit(20,"mm"),
                legend_height = unit(30,"mm")
            ),
            # heatmap size 
            width  = unit(CellSize*ncol(DistanceMatrix),"mm"), height = unit(CellSize*nrow(DistanceMatrix),"mm")
        )
        #-----------------------------------------------------------------------
        return(DistHeatmap)
    }
#--------------------------------------------------------------------------------------------------#
#====================================================================================================================================================#

#===| DEG RESULT PLOTS |=============================================================================================================================#
#---/ Volcano plot /-------------------------------------------------------------------------------#
#' @description draw volcano-plot og DEG analysis result
#' @param DegAnalysisResult deg analysis result table
#' @param ColorUp           up-regulation gene color.   default = NULL ( = use pre-defined color )
#' @param ColorDown         down-regulation gene color. default = NULL ( = use pre-defined color )
#' @param ColorNormal       nono-DEG gene color.        default = NULL ( = use pre-defined color )
#' @param AxisColor         axis-line color
#' @param GridColor         background grid line color
#' @param AxisLabelColor    axis label font color. if NULL, 'AxisColor' will be used.
#' @param FoldChangeGuide   foldchange guideline. default = NULL. if NOT NULL, DEG-tag will be re-defined.
#' @param ScoreValue        score metric selection
#' @param ScoreGuide        score guideline. default = NULL. if NOT NULL, DEG-tag will be re-defined.
#' @param PointSize         scatter-plot point size. default = 1
#' @param AxisFontSize      axis-label font size. default = 10 
#' @param AxisTitleFontSize axis-title font size. default = 12
#' @param Width             plot width. unit = 'inch'.
#' @param Height            plot height unit = 'inch'. 
#' @param PointPch          scatter-plot point shape. either '19' or '21'. default = 19
#' @param PlotSave          save plot as file or not. default = TRUE
#' @param DisplayOnScreen   display plot on screen with X11. default = FALSE
#' @param PlotName          plot filename. if NULL, default name will be used.
#' @param PlotTitle         plot title. if NULL, no title on plot.
#' @export 
    WTS_plot.DEG.Analysis.Result.Volcano <- function( DegAnalysisResult, 
        ColorUp=NULL, ColorDown=NULL, ColorNormal=NULL, AxisColor="#2d2d2d", GridColor="#d6d6d6", AxisLabelColor=NULL,
        FoldChangeGuide=NULL, ScoreValue=c("pvalue","fdr","score"), ScoreGuide=NULL,
        PointSize=3, AxisFontSize=10, AxisTitleFontSize=12, Width=8, Height=8, PointPch=19,
        PlotSave=TRUE, DisplayOnScreen=FALSE, PlotName=NULL, PlotTitle=NULL
    )
    {
        suppressPackageStartupMessages(library("ggplot2"))
        suppressPackageStartupMessages(library("dplyr"))
        # pre-define colors ----------------------------------------------------
        UP_GENE_COLOR   <- "#9f3122"
        DOWN_GENE_COLOR <- "#347019"
        NON_DEG_COLOR   <- "#d6d6d6"
        # score metric and y-limit ---------------------------------------------
        if( ScoreValue == "pvalue" )
        {
            VolcanoInput <- DegAnalysisResult[which(!is.na(DegAnalysisResult$pvalue) & DegAnalysisResult$keep), ]
            if( length(which(VolcanoInput$pvalue == 0)) > 0 )
            {
                max_y <- max(-log10(VolcanoInput[which(VolcanoInput$pvalue != 0), "pvalue"]))
                max_y <- max_y * 1.3
                VolcanoInput$y_value <- sapply( VolcanoInput$pvalue, function(p) if( p == 0 ){ max_y }else{ -log10(p) } )
            }else{
                VolcanoInput$y_value <- sapply( VolcanoInput$pvalue, function(p) if( p == 0 ){ max_y }else{ -log10(p) } )
                max_y <- max(VolcanoInput$y_value)*1.3
            }            
        }else if( ScoreValue == "fdr" )
        {
            VolcanoInput <- DegAnalysisResult[which(!is.na(DegAnalysisResult$pvalue) & DegAnalysisResult$keep), ]
            if( length(which(VolcanoInput$padj == 0)) > 0 )
            {
                max_y <- max(-log10(VolcanoInput[which(VolcanoInput$padj != 0), "padj"]))
                max_y <- max_y * 1.3
                VolcanoInput$y_value <- sapply( VolcanoInput$padj, function(p) if( p == 0 ){ max_y }else{ -log10(p) } )
            }else{
                VolcanoInput$y_value <- sapply( VolcanoInput$padj, function(p) if( p == 0 ){ max_y }else{ -log10(p) } )
                max_y <- max(VolcanoInput$y_value)*1.3
            }
        }else{
            VolcanoInput         <- DegAnalysisResult[which(!is.na(DegAnalysisResult$pvalue) & DegAnalysisResult$keep), ]
            VolcanoInput$y_value <- abs(VolcanoInput$score)
            max_y                <- max(VolcanoInput$y_value) * 1.3
        }
        # x-limit --------------------------------------------------------------
        max_x <- as.integer(max(abs(VolcanoInput$log2_foldchange))) * 1.2 
        # guideline or re-define DEGs ------------------------------------------
        if( is.null(FoldChangeGuide) )
        {
            FoldChangeGuide <- VolcanoInput %>% filter( keep ) %>% filter( deg != "" ) %>% .$log2_foldchange %>% abs() %>% min()
        }else{
            VolcanoInput[which(abs(VolcanoInput$log2_foldchange) < FoldChangeGuide), "deg" ] <- ""
        }
        ##
        if( is.null(ScoreGuide) )
        {
            ScoreGuide <- VolcanoInput %>% filter( keep ) %>% filter( deg != "" ) %>% .$y_value %>% min()
        }else{
            VolcanoInput[which(abs(VolcanoInput$score) < ScoreGuide), "deg" ] <- ""
        }
        # gene colors ----------------------------------------------------------
        if( is.null(ColorUp)     ){ UpGeneColor   <- UP_GENE_COLOR   }else{ UpGeneColor   <- ColorUp     }
        if( is.null(ColorDown)   ){ DownGeneColor <- DOWN_GENE_COLOR }else{ DownGeneColor <- ColorDown   }
        if( is.null(ColorNormal) ){ NonDegColor   <- NON_DEG_COLOR   }else{ NonDegColor   <- ColorNormal }
        #-----------------------------------------------------------------------
        VolcanoInput$deg_type <- sapply( VolcanoInput$deg, function(g)
        {
            if( g == "UP" ){ gt <- "UP-Regulation" }else if( g == "DOWN" ){ gt <- "DOWN-Regulation" }else{ gt <- "Non-DEGs" }
            return(gt)
        })
        GENE_COLORS <- c( "UP-Regulation"=UpGeneColor, "DOWN-Regulation"=DownGeneColor, "Non-DEGs"=NonDegColor )
        #-----------------------------------------------------------------------
        if( is.null(AxisLabelColor) ){ AxisLabelColor <- AxisColor }
        #-----------------------------------------------------------------------
        if( PointPch == 21 )
        {
            VPLOT <- ggplot( VolcanoInput, aes(x=log2_foldchange, y=y_value, fill=deg_type) ) + 
                geom_point( size = PointSize, pch=PointPch, color="#f4f4f4", alpha=0.9 ) +
                labs( x="FoldChange", y="Score Metric", fill="Genes" ) +
                scale_fill_manual( values = GENE_COLORS )
        }else{
            VPLOT <- ggplot( VolcanoInput, aes(x=log2_foldchange, y=y_value, color=deg_type) ) + 
                geom_point( size = PointSize, pch=PointPch, alpha=0.6 ) +
                labs( x="FoldChange", y="Score Metric", color="Genes" ) +
                scale_color_manual( values = GENE_COLORS )
        }
        #-----------------------------------------------------------------------
        if( is.null(PlotTitle) )
        {
            message("|---> no plot tile.")
        }else{
            VPLOT <- VPLOT + ggtitle( PlotTitle )
        }
        #-----------------------------------------------------------------------
        VPLOT <- VPLOT +
            coord_cartesian( xlim = c(-max_x, max_x), ylim = c(0, max_y) ) +
            geom_hline(yintercept=ScoreGuide, linetype=2, color='#1F78B4') +
            geom_vline(xintercept=c( -FoldChangeGuide, FoldChangeGuide ), linetype=2, color='#1F78B4') +
            theme(
                panel.background = element_rect(fill='white'), panel.grid.major = element_line(color = GridColor, linetype = 'dotted'),
                axis.text        = element_text(family="Helvetica", colour=AxisLabelColor, size=AxisFontSize),
                axis.title       = element_text(family="Helvetica", colour=AxisLabelColor, size=AxisTitleFontSize),
                axis.line        = element_line(colour=AxisColor),
                legend.title     = element_text(family="Helvetica", colour=AxisLabelColor, size=AxisTitleFontSize),
                legend.text      = element_text(family="Helvetica", colour=AxisLabelColor, size=AxisTitleFontSize),
                plot.margin      = unit(c(Height*0.05, Width*0.05, Height*0.05, Width*0.03), "in"),
                plot.title       = element_text(family="Helvetica", colour=AxisLabelColor, size=AxisTitleFontSize*1.1)
            )
        #-----------------------------------------------------------------------
        if( PlotSave )
        {
            if( is.null(PlotName) )
            {
                message("|---> no plot filename. default name will be used.")
                PlotName <- "./DEG.analysis.result.volcano.png"
            }
            message("|---> save plot as png file.")
            suppressWarnings(ggsave( VPLOT, file=PlotName, width=Width, height=Height, unit='in', dpi=150, type='cairo'))
        }
        #-----------------------------------------------------------------------
        if( DisplayOnScreen )
        {
            x11( width=Width, height=Height, bg="white")
            print(VPLOT)
        }
        #-----------------------------------------------------------------------
        return(VPLOT)
    }
#--------------------------------------------------------------------------------------------------#
#---/ CPM correlation pair-plots /-----------------------------------------------------------------#
#' @description draw correlation and scatter pair plots with CPM.
#' @param CpmMatrix                 CPM value matrix. REQUIRD.
#' @param SampleList                sample list. if NULL, all samples of input CPM matrix will be used.
#' @param CorrFontSize              correlation label font size. default = 2
#' @param CorrFonColors             correlation label font color. default = 'black'
#' @param ScatterPointSize          scatter plot point size. default = 0.6
#' @param ScatterPointColors        scatter plot point color. if NULL, pre-defined color will be used.
#' @param ScatterPointTransparency  scatter plot point transparency. default = 0.7 
#' @param SampleLabelFontSize       sample name font size. default = 2 
#' @param PlotTitle                 plot title. if NULL, no title.
#' @param PlotSave                  plot save as png file. default = TRUE 
#' @param PlotName                  plot file name. if NULL, pre-defined formatted name will be used.
#' @param DisplayOnScreen           display plot on screen. default = FALSE
#' @param PlotWidth                 plot png file width. unit = 'inch', default = 12
#' @param PlotHeight                plot png file height. unit = 'inch', default = 12
#' @export 
    WTS_plot.CPM.Correlation.Pairs <- function( CpmMatrix=NULL, SampleList=NULL, 
        CorrFontSize=2, CorrFonColors=NULL, ScatterPointSize=0.6, ScatterPointColors=NULL, ScatterPointTransparency=0.7, SampleLabelFontSize=1, 
        PlotTitle=NULL, PlotSave=TRUE, PlotName=NULL, DisplayOnScreen=FALSE, PlotWidth=11, PlotHeight=10  
    )
    {
        # cpm matrix input check -----------------------------------------------
        if( is.null(CpmMatrix) ){ stop("|---!!! no input CPM matrix. REQUIRED. please check again. STOPPED.") }
        # sample list check ----------------------------------------------------
        if( is.null(SampleList) )
        { 
            message("|--->> no sample list. all samples in input CPM matrix will be used.")
            SampleList <- colnames(CpmMatrix)    
        }
        # pre-define colors ----------------------------------------------------
        scatter_point_colors    <- "#cccccc"
        sample_label_font_color <- "#000000"
        corr_label_font_color   <- "#000000"
        # sample check ---------------------------------------------------------
        if( !all(SampleList %in% colnames(CpmMatrix)) ){ stop("|---!!! not all samples in cpm matrix. please check again. STOPPED.") }
        # colors ---------------------------------------------------------------
        if( is.null(CorrFonColors)         ){ CORR_FONT_COLOR <- corr_label_font_color   }else{ CORR_FONT_COLOR <- CorrFonColors         }   
        if( is.null(ScatterPointColors)    ){ POINT_COLOR     <- scatter_point_colors    }else{ POINT_COLOR     <- ScatterPointColors    }   
        # transparency ---------------------------------------------------------
        if( is.null(ScatterPointTransparency) ){ ALPHA <- 0.7 }else{ ALPHA <- ScatterPointTransparency }
        # point size -----------------------------------------------------------
        if( is.null(ScatterPointSize) ){ POINT_SZIE <- 0.7 }else{ POINT_SZIE <- ScatterPointSize }
        # functions ------------------------------------------------------------
        panel.corr <- function(x, y)
        {
            par(usr = c(0, 1, 0, 1))
            r   <- round(cor(x, y), digits=3)
            txt <- paste0("Corr:\n", r)
            text( 0.5, 0.5, txt, cex = CorrFontSize, col=CORR_FONT_COLOR ) 
        }
        panel.scat <- function(x, y)
        {
            points(x, y, pch = 19, cex = POINT_SZIE, col = alpha(POINT_COLOR, ALPHA))
        }
        # input data -----------------------------------------------------------
        InputData <- as.data.frame(CpmMatrix[ , SampleList ])
        lable_names         <- names(SampleList)
        names(lable_names)  <- SampleList
        colnames(InputData) <- lable_names[colnames(InputData)]
        maxDataValue = round(max(InputData), 0)
        # plot save ------------------------------------------------------------
        if( PlotSave )
        {
            if( is.null(PlotName) ){ PlotName <- "./CPM.correlation.scatter.pair.plot.png" }
            # draw plot --------------------------------------------------------
            png( PlotName, width=PlotWidth, height=PlotHeight, units="in",res=150 )
            pairs( 
                InputData,
                lower.panel = panel.scat,
                upper.panel = panel.corr,
                cex.labels  = SampleLabelFontSize,
                gap         = 0.5,
                main        = PlotTitle#,
                #xlim        = c(0, maxDataValue),
                #ylim        = c(0, maxDataValue)
            )
            dev.off()
        }
        # plot display ---------------------------------------------------------
        if( DisplayOnScreen )
        {
            x11( width=PlotWidth, height=PlotHeight, bg="white")
            pairs( 
                InputData,
                lower.panel = panel.scat,
                upper.panel = panel.corr,
                cex.labels  = SampleLabelFontSize,
                gap         = 0.5,
                main        = PlotTitle
            )
        }
        #----------------------------------------------------------------------#
        return("...pair plot created.")
    }
#--------------------------------------------------------------------------------------------------#
#====================================================================================================================================================#

#===| PATHWAY ANALYSIS PLOTS |=======================================================================================================================#
#---/ Top-ranked result dotplot or barplot /-------------------------------------------------------#
#' @description draw deg-enrich-pathway-analysis top-ranked result as dotplot or barplot
#' @param PathwayAnalysisResult pathway analysis result table
#' @param GeneType              analysis gene type. either 'UP','DOWN' or 'ALL_DEGs'. default = 'UP'.
#' @param PlotType              plot type. either 'dot' or 'bar'. default ='dot'
#' @param TopReults_N           number of top ranked pathways to draw plot
#' @param ColorUp               UP-DEGs analysis result colors. if NULL, pre-defined color will be used.  
#' @param ColorDown             DOWN-DEGs analysis result colors. if NULL, pre-defined color will be used.
#' @param ColorMerged           UP.DOWN merged DEGs analysis result colors. if NULL, pre-defined color will be used.
#' @param AxisColor             axis-line color
#' @param GridColor             background grid line color
#' @param AxisLabelColor        axis label font color. if NULL, 'AxisColor' will be used. 
#' @param ScoreMetric           score metric for ranking. either 'score','pvalue' or 'enrich_factor'. default = 'score'
#' @param AxisFontSize          axis-label font size. default = 10 
#' @param AxisTitleFontSize     axis-title font size. default = 12
#' @param Width                 plot width. unit = 'inch'.
#' @param Height                plot height unit = 'inch'. 
#' @param PathwayNameLength     length limit of displaying pathway name strings
#' @param PlotSave              save plot as file or not. default = TRUE
#' @param DisplayOnScreen       display plot on screen with X11. default = FALSE
#' @param PlotName              plot filename. if NULL, default name will be used.
#' @param PlotTitle             plot title. if NULL, no title on plot.
#' @export 
    WTS_plot.Pathway.Analysis.TopRank.Result <- function( PathwayAnalysisResult, GeneType="UP", PlotType="dot", TopReults_N=15, ScoreMetric="score", 
        ColorUp=NULL, ColorDown=NULL, ColorMerged=NULL, AxisColor="#2d2d2d", GridColor="#d6d6d6", AxisLabelColor=NULL,         
        AxisFontSize=10, AxisTitleFontSize=12, Width=8, Height=8, PathwayNameLength=50,
        PlotSave=TRUE, DisplayOnScreen=FALSE, PlotName=NULL, PlotTitle=NULL    
    )
    {
        suppressPackageStartupMessages(library("ggplot2"))
        suppressPackageStartupMessages(library("stringr"))
        suppressPackageStartupMessages(library("dplyr"))
        if( GeneType %in% PathwayAnalysisResult$input_geneset )
        {
            # scale down function --------------------------------------------------
            ScaleDown <- function( x, y )
            {
                cv <- (x-min(x))*((max(y)-min(y))/(max(x)-min(x))) + min(y)
                cv <- round(cv, 2)
                return(cv)
            }
            # pre-defined color ----------------------------------------------------
            UP_GENE_PATH_COLOR     <- c("#9f3122","#cf9890")
            DOWN_GENE_PATH_COLOR   <- c("#347019","#99b78c")
            MERGED_GENE_PATH_COLOR <- c("#1a5387","#8ca9c3")
            #-----------------------------------------------------------------------
            InputDataRaw <- PathwayAnalysisResult %>% filter( input_geneset == GeneType )
            # score metric and plot-data -------------------------------------------
            if( ScoreMetric == "pvalue" )
            {
                InputData <- InputDataRaw %>% arrange( pvalue, dplyr::desc(score) ) %>% dplyr::slice(1:TopReults_N)
                if( length(which(InputData$pvalue == 0)) > 0 ){
                    max_pvale_score <- max(-log10(InputData[which(InputData$pvalue != 0), "pvalue"])) * 1.2
                }
                InputData$RANK_VALUE <- sapply(InputData$pvalue, function(pv) ifelse( pv == 0, max_pvale_score, -log10(pv) ))
            }else{
                InputData <- InputDataRaw %>% arrange( dplyr::desc(score), pvalue ) %>% dplyr::slice(1:TopReults_N)
                InputData$RANK_VALUE <- InputData$score
            }
            # dotplot only params --------------------------------------------------
            if( PlotType == "dot" )
            {
                InputData$COLOR_FACTOR <- ScaleDown(x=InputData$enrich_factor, y=c(2,10))
                InputData$SIZE_FACTOR  <- InputData$enrich_genes
            }
            #----------------------------------------------------------------------#
            InputData$pathway <- gsub("_", " ", InputData$pathway)
            InputData$pathway <- factor(InputData$pathway, levels=InputData$pathway[TopReults_N:1])
            # colors I -------------------------------------------------------------
            if( is.null(ColorUp)     ){ COLOR_UP     <- UP_GENE_PATH_COLOR[1]     }else{ COLOR_UP     <- ColorUp     }  
            if( is.null(ColorDown)   ){ COLOR_DOWN   <- DOWN_GENE_PATH_COLOR[1]   }else{ COLOR_DOWN   <- ColorDown   }  
            if( is.null(ColorMerged) ){ COLOR_MERGED <- MERGED_GENE_PATH_COLOR[1] }else{ COLOR_MERGED <- ColorMerged }    
            # colors II ------------------------------------------------------------
            COLOR_UP_MID     <- colorRampPalette( c("#ffffff", COLOR_UP)     )(10)[5]
            COLOR_DOWN_MID   <- colorRampPalette( c("#ffffff", COLOR_DOWN)   )(10)[5]
            COLOR_MERGED_MID <- colorRampPalette( c("#ffffff", COLOR_MERGED) )(10)[5]
            # base plots -----------------------------------------------------------
            if( PlotType == "bar" )
            {
                # bar colors -------------------------------------------------------
                if( all(unique(InputData$input_geneset) == "UP") ){ RANK_COLOR <- COLOR_UP 
                }else if( all(unique(InputData$input_geneset) == "DOWN") ){ RANK_COLOR <- COLOR_DOWN 
                }else{ RANK_COLOR <- COLOR_MERGED }
                # bar plot ---------------------------------------------------------
                PlotObj <- ggplot(InputData, aes( x=RANK_VALUE, y=pathway) ) +
                    geom_bar(stat='identity', position='dodge', width=0.9, fill=RANK_COLOR ) +
                    labs( x="SCORE", y="" ) +
                    theme( legend.position = "none" )    
            }else{
                # dot plot ---------------------------------------------------------
                PlotObj <- ggplot(InputData) +
                    geom_point(aes( x=RANK_VALUE, y=pathway, color=COLOR_FACTOR, size=SIZE_FACTOR) ) +
                    labs( x="SCORE", y="", color="Enrichment", size="Genes") +
                    scale_size( range=c(4,12) )
                # dot colors -------------------------------------------------------
                if( all(unique(InputData$input_geneset) == "UP") ){ 
                    PlotObj <- PlotObj + scale_color_gradient2('Enrichment', low='#f4f4f4', mid=COLOR_UP_MID,     high=COLOR_UP,     midpoint=4 ) 
                }else if( all(unique(InputData$input_geneset) == "DOWN") ){ 
                    PlotObj <- PlotObj + scale_color_gradient2('Enrichment', low='#f4f4f4', mid=COLOR_DOWN_MID,   high=COLOR_DOWN,   midpoint=4 ) 
                }else{ 
                    PlotObj <- PlotObj + scale_color_gradient2('Enrichment', low='#f4f4f4', mid=COLOR_MERGED_MID, high=COLOR_MERGED, midpoint=4 ) 
                }
            }
            # plot title -----------------------------------------------------------
            if( !is.null( PlotTitle ) ){ PlotObj <- PlotObj + ggtitle( PlotTitle ) }
            # axis-label color -----------------------------------------------------
            if( is.null(AxisLabelColor) ){ AxisLabelColor <- AxisColor }
            # final plot object ----------------------------------------------------
            PlotObj <- PlotObj + scale_y_discrete(labels=function(lbl) str_wrap(lbl, width = PathwayNameLength)) +
                theme(
                    panel.background = element_rect(fill='white'), 
                    panel.grid.major = element_line(color = GridColor, linetype = 'dotted'),
                    axis.text.y      = element_text(family="Nunito", colour=AxisLabelColor, size=AxisFontSize, face="bold"),
                    axis.text.x      = element_text(family="Nunito", colour=AxisLabelColor, size=AxisFontSize),
                    axis.title       = element_text(family="Nunito", colour=AxisLabelColor, size=AxisTitleFontSize, , face="bold"),
                    axis.line        = element_line(colour=AxisColor),
                    legend.title     = element_text(family="Nunito", colour=AxisLabelColor, size=AxisTitleFontSize*0.9),
                    legend.text      = element_text(family="Nunito", colour=AxisLabelColor, size=AxisTitleFontSize*0.8),
                    plot.margin      = unit(c(Height*0.02, Width*0.04, Height*0.04, Width*0.02), "in"),
                    plot.title       = element_text(family="Nunito", colour=AxisLabelColor, size=AxisTitleFontSize*1.1),
                )
        }else{
            if( GeneType == "UP" ){ LABELS = "No Pathway Analysis Reusults of 'UP'-regulated Genes" 
            }else if( GeneType == "DOWN" ){ LABELS = "No Pathway Analysis Reusults of 'DOWN'-regulated Genes" 
            }else if( GeneType == "ALL_DEGs" ){ LABELS = "No Pathway Analysis Reusults of DEGs" }
            PlotObj <- ggplot() + annotate("text", x=0,y=0,label=LABELS) + theme_void()
        }
        #-----------------------------------------------------------------------
        if( PlotSave )
        {
            if( is.null(PlotName) )
            {
                message("|---> No input of plot filename. Default filename was used and it cound be overwritten on pre-existed file.")
                if( GeneType == "ALL_DEGs" ){ GeneTypePrefix <- "All.DEG" }else{ GeneTypePrefix <- toupper(GeneType) }
                PlotName <- sprintf("./%s_pathway.analysis.%s.genes.top%s.%splot.png", unique(InputData$deg_analysis_id), GeneTypePrefix, TopReults_N, PlotType )
            }
            suppressWarnings(ggsave( PlotObj, file=PlotName, width=Width, height=Height, unit='in', dpi=150, type='cairo'))
        }
        #-----------------------------------------------------------------------
        if( DisplayOnScreen )
        {
            x11( width=Width, height=Height, bg="white")
            print(PlotObj)
        }
        #-----------------------------------------------------------------------
        return(PlotObj)
    }
#--------------------------------------------------------------------------------------------------#
#---/ Top-ranked pathways and DEGs network /-------------------------------------------------------#
#' @description create network of top-ranked pathways and their DEGs 
#' @param TopRankNetData         input network node and edge data list
#' @param NetworkLayout          network layout. optimal=graphopt, large_net=lgl gem=gem, mid_net=dh, tree1=tree, tree2=layout_with_sugiyama. default = "optimal"
#' @param UseVisNetwork          create network using visNetwork instead igraph. Only work in browser.
#' @param VisNetHighlightDegree  highlight degree length when mouse-hovering. Only work in visNetwork.
#' @export 
    WTS_plot.Pathway.Analysis.TopRank.Network <- function( TopRankNetData, NetworkLayout="optimal", UseVisNetwork=FALSE, VisNetHighlightDegree=NULL
    )
    {
        suppressPackageStartupMessages(library("stringr"))
        suppressPackageStartupMessages(library("igraph"))
        suppressPackageStartupMessages(library("Hmisc"))
        # pre-defined colors ---------------------------------------------------
        Color_UpGene       <- "#dba9a7"
        Color_DownGene     <- "#b5c6df"
        ColorPathway       <- "#ebcd7f"
        ColorAllDegPathway <- "#fb8632"
        # create network -------------------------------------------------------
        if( is.null(TopRankNetData) )
        {
            stop("|---!!! no pathway analysis top-rank result network data list object found. REQUIRED. please check again. STOPPED.")
        }else{
            NET <- graph.data.frame(TopRankNetData$edge[,c("from","to")], directed=FALSE)
        }
        # extract nodes from network -------------------------------------------
        NET_node <- get.data.frame(NET, 'vertices')
        # node properties ------------------------------------------------------
            # node type
            NET_node$type <- TopRankNetData$node[match(NET_node$name, TopRankNetData$node$id), "type"]
            # node color
            NET_node$color <- sapply( NET_node$type, function(nt) 
            {
                if( nt == "pathway" ){ ColorPathway }else if( nt == "up_gene" ){ Color_UpGene }else{ Color_DownGene }
            })
            # node frame color
            NET_node$frame.color <- "#FFFFFF"
            up_down_both_pathway <- TopRankNetData$edge %>% dplyr::rename(pathway=1, id=2) %>% left_join( ., TopRankNetData$node, by="id") %>% 
                group_by(pathway) %>% reframe( deg_type = length(unique(type))) %>% 
                filter( deg_type > 1 )
            if( nrow(up_down_both_pathway) > 0 ){ NET_node[which(NET_node$name %in% up_down_both_pathway$pathway), "frame.color"] <- ColorAllDegPathway }
            # node shape 
            NET_node$shape <- sapply( NET_node$type, function(y) ifelse( y == "pathway", "square", "circle") )
            # node size 
            NET_node$size <- sapply( NET_node$type, function(y) ifelse( y == "pathway", 4, 3) )
            # node label
            NET_node$label        <- NET_node$name
            NET_node$label.family <- "sans"
            NET_node$label.color  <- "#0f0f0f"
            NET_node$label.font   <- sapply( NET_node$type, function(y) ifelse( y == "pathway", 2, 1) )
            NET_node$label.cex    <- 0.6
        # apply to network nodes -----------------------------------------------
        V(NET)$label        <- stringr::str_wrap(NET_node[match(V(NET)$name, NET_node$name), "label"], 15)
        V(NET)$color        <- NET_node[match(V(NET)$name, NET_node$name), "color"       ]
        V(NET)$frame.color  <- NET_node[match(V(NET)$name, NET_node$name), "frame.color" ]
        V(NET)$shape        <- NET_node[match(V(NET)$name, NET_node$name), "shape"       ]
        V(NET)$label.family <- NET_node[match(V(NET)$name, NET_node$name), "label.family"]
        V(NET)$size         <- NET_node[match(V(NET)$name, NET_node$name), "size"        ]
        V(NET)$label.cex    <- NET_node[match(V(NET)$name, NET_node$name), "label.cex"   ]
        V(NET)$label.color  <- NET_node[match(V(NET)$name, NET_node$name), "label.color" ]
        V(NET)$label.font   <- NET_node[match(V(NET)$name, NET_node$name), "label.font"  ]
        # edge properties ------------------------------------------------------
        E(NET)$lty         <- 1
        E(NET)$color       <- "#cccccc"
        # network layout -------------------------------------------------------
        if(       NetworkLayout == "optimal"   ){ 
            LAYOUT       <- layout_with_graphopt(NET, niter=1500, charge=0.01)
            VisNetLayout <- "layout_with_graphopt" 
        }else if( NetworkLayout == "large_net" ){ 
            LAYOUT       <- layout_with_lgl(NET, repulserad=vcount(NET)^2 )
            VisNetLayout <- "layout_with_lgl" 
        }else if( NetworkLayout == "gem"       ){ 
            LAYOUT       <- layout_with_gem(NET, maxiter = 100*vcount(NET)^2, temp.min=1/20 ) 
            VisNetLayout <- "layout_with_gem"
        }else if( NetworkLayout == "mid_net"   ){ 
            LAYOUT <- layout_with_dh(NET, maxiter=20, fineiter=20)            
            VisNetLayout <- "layout_with_dh"
        }else if( NetworkLayout == "tree1"     ){ 
            LAYOUT       <- layout_as_tree(NET, mode="all")
            V(NET)$label <- stringr::str_wrap(NET_node[match(V(NET)$name, NET_node$name), "label"], 10)
            VisNetLayout <- "layout_as_tree"
        }else if( NetworkLayout == "tree2"     ){ 
            LAYOUT       <- layout_with_sugiyama(NET, maxiter=500, vgap=4, hgap=10)
            V(NET)$label <- stringr::str_wrap(NET_node[match(V(NET)$name, NET_node$name), "label"], 10)
            VisNetLayout <- "layout_with_sugiyama"
        }     
        # make network object --------------------------------------------------
        if( UseVisNetwork )
        {
            suppressPackageStartupMessages(library("visNetwork"))
            # visnetwork nodes -------------------------------------------------
            VISNET_node <- get.data.frame(NET, 'vertices')
            VISNET_node <- VISNET_node %>% dplyr::rename( id=name ) %>% mutate( label = gsub("\\\n", " ", label) )
            VISNET_node$shape <- sapply( NET_node$type, function(y) ifelse( y == "pathway", "dot", "ellipse") )
            VISNET_node$value <- sapply( NET_node$type, function(y) ifelse( y == "pathway", 4, 3) )
            # visnetwork edges -------------------------------------------------
            VISNET_edge <- get.data.frame(NET, 'edges')
            VISNET_edge$shadow <- FALSE 
            # draw network -----------------------------------------------------
            NETWORK <- visNetwork( VISNET_node, VISNET_edge ) %>% visIgraphLayout( layout = VisNetLayout )
            # hightlight edges option ------------------------------------------
            if( all(!is.null(VisNetHighlightDegree) & is.numeric(VisNetHighlightDegree)) )
            { NETWORK <- NETWORK %>% visOptions(highlightNearest = list(enabled = T, degree = VisNetHighlightDegree, hover = T, labelOnly=F)) }
        }else{
            # draw network -----------------------------------------------------
            NETWORK <- plot( NET , vertex.label.degree = pi/2, layout = LAYOUT, asp=0)
        }
        return(NETWORK)

    }
#--------------------------------------------------------------------------------------------------#
#====================================================================================================================================================#

#===| GASEA PLOTS |==================================================================================================================================#
#---/ GSEA plots /---------------------------------------------------------------------------------#
#' @description draw GSEA plot of selected pathway
#' @param AnalysisID      DEG analysis id 
#' @param InputGeneRank   gene expression profile rank. named vector. names should be entrez id.
#' @param PathwayID       pathway-id. 
#' @param PathwayName     pathway-name. This will be used also as PlotTitle.
#' @param FontSize        label font size. default = 10
#' @param Width           plot width. unit = 'inch'.
#' @param Height          plot height unit = 'inch'. 
#' @param ColorUp         UP-DEGs analysis result colors. if NULL, pre-defined color will be used.  
#' @param ColorDown       DOWN-DEGs analysis result colors. if NULL, pre-defined color will be used.
#' @param ColorLines      GSEA NES line color. if NULL, pre-defined color will be used.
#' @param AxisColor       axis-line color  
#' @param GridColor       background grid line color
#' @param AxisLabelColor  axis label font color. if NULL, 'AxisColor' will be used. 
#' @param PlotSave        save plot as file or not. default = TRUE
#' @param DisplayOnScreen display plot on screen with X11. default = FALSE
#' @param PlotName        plot filename. if NULL, default name will be used.
#' @export 
    WTS_plot.Top.Rank.Profile.GSEA.plot <- function( AnalysisID, InputGeneRank, 
        PathwayID="pw01", PathwayName="cpdb", 
        FontSize=12, Width=6, Height=7, 
        ColorUp=NULL, ColorDown=NULL, ColorLines=NULL, AxisColor="#2d2d2d", GridColor="#d6d6d6", AxisLabelColor=NULL, 
        PlotSave=TRUE, DisplayOnScreen=FALSE, PlotName=NULL
    ) 
    {       
        suppressPackageStartupMessages(library("ggplot2"))
        suppressPackageStartupMessages(library("fgsea"))
        suppressPackageStartupMessages(library("RMySQL"))
        suppressPackageStartupMessages(library("stringr")) 
        suppressPackageStartupMessages(library("dplyr"))
        # pre-defined color ----------------------------------------------------
        UP_GENE_COLOR   <- "#9f3122"
        DOWN_GENE_COLOR <- "#14426c"
        LINE_COLOR      <- "#4DAF4A"
        # database connection --------------------------------------------------
        source("/data/wts/params/ruo_wts_db.R")
        dbCon <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db='public_database' )
        #on.exit(dbDisconnect(dbCon))
        # pathway set id -------------------------------------------------------
        if( all(is.null(PathwayID) & is.null(PathwayName)) )
        { 
            dbDisconnect(dbCon)
            stop("|---!!! STOPPED. no both pathway set id and pathway set name. at least one of the pathway set identifier shoud be exist.") 
        }else{
            if( is.null(PathwayID) )
            {
                PathwayID <- dbGetQuery(dbCon, sprintf("SELECT path_id FROM pathwaysets WHERE pathway = '%s'", PathwayName))$path_id
                if( length( PathwayID) != 1 ){ 
                    dbDisconnect(dbCon)
                    stop("|---!!! STOPPED. no pathway set has such name or multiple-pathway sets found." ) 
                }
            }
        }
        # pathway set ----------------------------------------------------------
        PathwaySet     <- dbGetQuery(dbCon, sprintf("SELECT path_id,pathway,entrez FROM pathwaysets WHERE path_id = '%s'", PathwayID))
        PATHWAYSET     <- unlist(strsplit(PathwaySet$entrez, ";"))
        PathwayName    <- PathwaySet$pathway
        dbDisconnect(dbCon)

        PathwayName <- gsub("GOBP_", "", PathwayName)
        PathwayName <- gsub("_", " ", PathwayName)

        # colors ---------------------------------------------------------------        
        if( is.null(ColorUp)    ){ COLOR_UP   <- UP_GENE_COLOR   }else{ COLOR_UP <- ColorUp     }
        if( is.null(ColorDown)  ){ COLOR_DOWN <- DOWN_GENE_COLOR }else{ COLOR_DOWN <- ColorDown }
        if( is.null(ColorLines) ){ COLOR_LINE <- LINE_COLOR      }else{ COLOR_LINE <- ColorLines }
        #-----------------------------------------------------------------------
        rnk <- rank(-InputGeneRank) ; ord = order(rnk)
        #-----------------------------------------------------------------------
        gseaParam <- 1
        lfcAdj    <- InputGeneRank[ord]
        LFC.adj   <- lfcAdj / max(abs(lfcAdj))
        #-----------------------------------------------------------------------
        pathway.stat <- unname(as.vector(na.omit(match(PATHWAYSET, names(LFC.adj)))))
        pathway.stat <- sort(pathway.stat)
        #-----------------------------------------------------------------------
        gseaRes <- calcGseaStat(stats = LFC.adj, selectedStats = pathway.stat, gseaParam = 1, returnAllExtremes = TRUE)
        bottoms <- gseaRes$bottoms 
        tops    <- gseaRes$tops
        #-----------------------------------------------------------------------
        n      <- length(LFC.adj)
        xs     <- as.vector(rbind(pathway.stat - 1, pathway.stat))
        ys     <- as.vector(rbind(bottoms, tops))
        toPlot <- data.frame(x=c(0, xs, n + 1), y=c(0, ys, 0))
        diff   <- (max(tops) - min(bottoms)) / 8
        #-----------------------------------------------------------------------
        gseaBreakcs <- round(seq(as.integer(range(toPlot$y)[1]*10)/10, as.integer(range(toPlot$y)[2]*10)/10, 0.1), 1)
        #-----------------------------------------------------------------------
        fc2plot       <- data.frame(x = rnk, lfc = InputGeneRank)
        fc2plot$col_y <- sapply( fc2plot$lfc, function(q) ifelse( q >= 0,  COLOR_UP, COLOR_DOWN))
        fc2plot$y     <- min(bottoms)
        fc2plot$lfc2  <- sapply( fc2plot$lfc, function(q) {
            if( q >= 1.5 ){ sv = 1.5 }else if( q < 1.5 & q > -1.5){ sv = q }else if( q <= -1.5 ){ sv = -1.5 }
            return(sv)
        })
        if( max(abs(fc2plot$lfc)) == max(fc2plot$lfc) )
        {
            nl           <- (min(fc2plot$lfc)/max(fc2plot$lfc))
            fc2plot$lfc3 <- scales::rescale(fc2plot$lfc, to=c(nl, 1) )
        }else{
            pl           <- -(max(fc2plot$lfc)/min(fc2plot$lfc))
            fc2plot$lfc3 <- scales::rescale(fc2plot$lfc, to=c(-1, pl) )
        }
        #-----------------------------------------------------------------------
        x = y = NULL
        if( is.null(AxisLabelColor) ){ AxisLabelColor <- AxisColor }
        # plot A ---------------------------------------------------------------
        g1 = ggplot(toPlot, aes( x = x, y = y) ) +
            geom_point( color = COLOR_LINE, size=0.1) +
            geom_line(  color = COLOR_LINE, size =2 ) +
            geom_hline( yintercept = max(tops),    colour = COLOR_UP,   size = 0.9, linetype="dashed") +
            geom_hline( yintercept = min(bottoms), colour = COLOR_DOWN, size = 0.9, linetype="dashed") +
            geom_hline( yintercept = 0, colour = AxisColor ) + 
            theme_bw() +
            ggtitle( str_wrap(PathwayName, width = 40)  ) +
            scale_y_continuous(breaks = gseaBreakcs ) +
            geom_segment(data=data.frame(x=pathway.stat), mapping=aes(x=x, y=min(bottoms)-0.2, xend=x, yend=min(bottoms)), size=0.8 )  +
            theme(
                axis.text.y  = element_text(family="Helvetica", colour=AxisLabelColor, size=FontSize-2, hjust=0.25),
                axis.text.x  = element_text(family="Helvetica", colour=AxisLabelColor, size=FontSize-3, vjust = 1, hjust = 0.5),
                axis.title.y = element_text(family="Helvetica", face='bold', colour=AxisLabelColor, size=FontSize, vjust=1),
                axis.title.x = element_text(family="Helvetica", face='bold', colour=AxisLabelColor, size=FontSize, vjust=1),
                axis.line.y  = element_line(colour = AxisColor),
                plot.title   = element_text(family="Helvetica", face='bold', colour=AxisLabelColor, size=FontSize+3),
                panel.border = element_blank(),
                panel.grid.minor = element_blank()
            ) +
            labs(x="Rank in Ordered Dataset", y="Enrichment Score") 

        g2 = g1 + geom_segment(data=fc2plot, mapping=aes(x=x, y=y-0.2 , xend=x, yend=y-0.13, color = lfc2), size = 0.3) +
            scale_color_gradient2('Rank Stats', low=COLOR_DOWN, mid='#ffffff', high=COLOR_UP, midpoint=0, breaks = c(-1e10, 1e10))      
        # plot B ---------------------------------------------------------------
        g3 = ggplot() + 
            geom_segment( data=fc2plot, mapping=aes(x=x, y=lfc3 , xend=x, yend=0 ), size = 0.8, color ='#d6d6d6') + 
            geom_hline( yintercept = 0, colour = AxisColor )+
            coord_cartesian( ylim = c(-1, 1) ) +
            theme(
                axis.text.y  = element_text(family="Helvetica", colour=AxisLabelColor, size=FontSize-2, hjust=0.25),
                axis.text.x  = element_text(family="Helvetica", colour=AxisLabelColor, size=FontSize-3, vjust = 1, hjust = 0.5),
                axis.title.y = element_text(family="Helvetica", face='bold', colour=AxisLabelColor, size=FontSize, vjust=1),
                panel.border = element_blank(),
                panel.background = element_rect(fill = '#ffffff'),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_line(color = GridColor, linetype = 'dotted'),
                axis.line.y  = element_line(colour = AxisColor)
            ) +
            labs( x = '', y = 'Ranked Metric (Scaled)')
        # merge plot -----------------------------------------------------------
        gseaPlot <- cowplot::plot_grid(g2, g3, ncol=1, rel_heights = c(1,0.4))
        #-----------------------------------------------------------------------
        if( PlotSave )
        {
            if( is.null(PlotName) )
            {
                message("|--->> no plot filename. default name will be used.")
                PlotName <- sprintf("./%s_GSEA.plot_%s.%s.png", AnalysisID, PathwayID, PathwayName)
            }
            message("|--->> GSEA plot saved as png file.")
            suppressWarnings(ggsave( gseaPlot, file=PlotName, width=Width, height=Height, unit='in', dpi=150, type='cairo'))
        }
        #-----------------------------------------------------------------------
        if( DisplayOnScreen )
        {
            x11( width=Width, height=Height, bg="white")
            print(gseaPlot)
        }
        #-----------------------------------------------------------------------
        return(gseaPlot)
    }
#--------------------------------------------------------------------------------------------------#
#====================================================================================================================================================#

#===| HEATMAP |======================================================================================================================================#
#' @description draw gene expression profiles (matrix) heatmap
#' @param ExpressionValueMatrix  input data matrix. if numeric data.frame, automatically converted into numeric matrix
#' @param Log2                   log2 transformation of input data. default = FALSE                       
#' @param Scale                  scaling of values. 'none'=no scaling. 'row'=scaling row direction. 'column'=scaling column direction. default = 'none'.
#' @param RowClustering          perform rwo-wise clustering. default = TRUE ( no dendrogem will be shown )
#' @param ColumnClustering       perform column-wise clustering. default = FALSE ( no dendrogem will be shown )
#' @param ColorUp                UP-regulated gene color. if NULL, pre-defined color will be used.
#' @param ColorDown              DOWN-regulated gene color. if NULL, pre-defined color will be used.
#' @param ColorMid               Middle Value color (generally 0). if NULL, pre-defined color will be used.
#' @param ColorNA                missing value color. if NULL, pre-defined color will be used.
#' @param ColorBorder            heatmap border and heatmap cell border color. if NULL, pre-defined color will be used.
#' @param HeatmapBorderLineWidth heatmap border line thickness. default = 0.5
#' @param HeatmapCellLineWidth   heatmap cell border line thickness. default = 0.3
#' @param RowTitle               row sample group title. defulat = NULL (none)
#' @param RowLabels              row labels. default = NULL ( = rownames of input data )  
#' @param RowLabelsSide          row labels location. either 'left' or 'right'. default = 'left' 
#' @param RowLabelsRotation      row labels roation angle. default = 0
#' @param RowLabelsSize          row labels font size. default = 10
#' @param RowLabelsColors        row labels font colors. default = NULL ( = black )
#' @param ColumnTitle            column sample group title. defulat = NULL (none) 
#' @param ColumnLabels           column labels. default = NULL ( = colnames of input data )    
#' @param ColumnLabelsSide       column labels location. either 'top' or 'bottom'. default = 'top'
#' @param ColumnLabelsRotation   column labels roation angle. default = 90
#' @param ColumnLabelsSize       column labels font size. default = 10
#' @param ColumnNamesCentered    column labels aligning to center of column. default = 'TRRUE'
#' @param CellWidth              heatmap cell width. unit = 'mm'. default = 20
#' @param CellHeight             heatmap cell height. unit = "mm", default = NA ( = use all height of canvas) 
#' @param TopAnnot               annotation located at top of heatmap
#' @param BottomAnnot            annotation located at bottom of heatmap
#' @param LeftAnnot              annotation located at left of heatmap
#' @param RightAnnot             annotation located at right of heatmap
#' @param ShowLegend             display legend or not. default = TRUE
#' @param LegendTitle            title of legend. if NULL, "VALUES" will be used.
#' @param LegendFontSize         legend font size 
#' @param LegendWidth            lengend width. unit ='mm', default = 25 
#' @param LegendHeight           lengend height. unit ='mm', default = 50 
#' @export 
    WTS_plot.Gene.Expression.Profiles.Heatmap <- function( ExpressionValueMatrix, 
        Log2=FALSE, Scales="none", RowClustering=TRUE, ColumnClustering=FALSE, 
        ColorUp=NULL, ColorDown=NULL, ColorMid=NULL, ColorNA=NULL, ColorBorder=NULL,
        HeatmapBorderLineWidth=0.5, HeatmapCellLineWidth=0.3,
        RowTitle=NULL, RowLabels=NULL, RowLabelsSide="left", RowLabelsRotation=0, RowLabelsSize=10, RowLabelsColors=NULL,
        ColumnTitle=NULL, ColumnLabels=NULL, ColumnLabelsSide="top", ColumnLabelsRotation=90, ColumnLabelsSize=10, ColumnNamesCentered=TRUE,
        CellWidth=15, CellHeight=NA,
        TopAnnot=NULL, BottomAnnot=NULL, LeftAnnot=NULL, RightAnnot=NULL,
        ShowLegend=TRUE, LegendTitle=NULL, LegendFontSize=10, LegendWidth=NA, LegendHeight=NA
    )
    {
        suppressPackageStartupMessages(library("ComplexHeatmap"))
        suppressPackageStartupMessages(library("circlize"))
        suppressPackageStartupMessages(library("dplyr"))
        # pre-defined color ----------------------------------------------------
        UP_GENE_COLOR       <- "#ff2020"
        DOWN_GENE_COLOR     <- "#208020"
        MIDDLE_COLOR        <- "#ffffff"   
        MISSING_VALUE_COLOR <- "#eaeaea"  
        BORDER_COLOR        <- "#727272"   
        # data format check ----------------------------------------------------
        if( !is.matrix(ExpressionValueMatrix) )
        {
            if( is.numeric(as.matrix(ExpressionValueMatrix)) ){ MAT <- as.matrix(ExpressionValueMatrix)
            }else{
                stop("!---!!! Input data is not 'matirx' and cannot be converted into numeric-matrix. Please check input data. STOPPED.")
            }
        }
        # log transformation ---------------------------------------------------
        if( Log2 )
        { MAT <- log2(ExpressionValueMatrix + 1) }else{ MAT <- ExpressionValueMatrix }
        # scaling --------------------------------------------------------------
        if( Scales == "row" ){
            MAT <- t(apply(MAT, 1, scale))
            colnames(MAT) <- colnames(ExpressionValueMatrix)
        }else if( Scales == "column" ){
            MAT <- apply(MAT, 2, scale)
            rownames(MAT) <- rownames(ExpressionValueMatrix)
        }else{
            MAT <- MAT
        }
        # data digits adjust ---------------------------------------------------
        MAT <- round(MAT, 4)
        # heatmap color --------------------------------------------------------
        if( is.null(ColorUp)    ){ COLOR_UP     <- UP_GENE_COLOR       }else{ COLOR_UP     <- ColorUp     }
        if( is.null(ColorDown)  ){ COLOR_DOWN   <- DOWN_GENE_COLOR     }else{ COLOR_DOWN   <- ColorDown   }
        if( is.null(ColorMid)   ){ COLOR_MID    <- MIDDLE_COLOR        }else{ COLOR_MID    <- ColorMid    }
        if( is.null(ColorNA)    ){ NA_COLOR     <- MISSING_VALUE_COLOR }else{ NA_COLOR     <- ColorNA     }
        if( is.null(ColorBorder)){ COLOR_BORDER <- BORDER_COLOR        }else if(is.na(ColorBorder)){ COLOR_BORDER <- NA }else{ COLOR_BORDER <- ColorBorder }
        # heatmap color function -----------------------------------------------
        data_value_range <- round(range(MAT),4)
        data_value_limit <- max(abs(data_value_range))
        HEATMAP_COLOR_FUN <- colorRamp2( 
            breaks       = c( -data_value_limit, 0, data_value_limit ), 
            colors       = c( COLOR_DOWN, COLOR_MID, COLOR_UP), 
            transparency = 0
        )
        # heatmap size ---------------------------------------------------------
        if(is.na(CellWidth) ){ WIDTH  <- NULL }else{ WIDTH  <- unit(CellWidth*ncol(MAT),"mm")  }
        if(is.na(CellHeight)){ HEIGHT <- NULL }else{ HEIGHT <- unit(CellHeight*nrow(MAT),"mm") }
        # row/column title -----------------------------------------------------
        if( is.null(RowTitle)    ){ ROW_TITLE    <- character(0) }else{ ROW_TITLE    <- RowTitle    }
        if( is.null(ColumnTitle) ){ COLUMN_TITLE <- character(0) }else{ COLUMN_TITLE <- ColumnTitle }
        # row/column label -----------------------------------------------------
        if( is.null(RowLabels)    ){ ROW_LABEL    <- rownames(MAT) }else{ ROW_LABEL    <- RowLabels    }
        if( is.null(ColumnLabels) ){ COLUMN_LABEL <- colnames(MAT) }else{ COLUMN_LABEL <- ColumnLabels }
        # row label colors -----------------------------------------------------
        if( is.null(RowLabelsColors) ){ ROW_LABEL_COLOR <- "#000000" }else{ ROW_LABEL_COLOR <- RowLabelsColors }
        # legend ---------------------------------------------------------------
        if( is.null(LegendTitle)    ){ LEGEND_TITLE     <- "VALUES" }else{ LEGEND_TITLE     <- LegendTitle    }
        if( is.null(LegendFontSize) ){ LEGEND_FONT_SIZE <- 10       }else{ LEGEND_FONT_SIZE <- LegendFontSize }
        if( is.na(LegendWidth)      ){ LEGEND_WIDTH     <- 20       }else{ LEGEND_WIDTH     <- LegendWidth    }
        if( is.na(LegendHeight)     ){ LEGEND_HEIGHT    <- 30       }else{ LEGEND_HEIGHT    <- LegendHeight   }
        LEGEND_TITLE_FONT_SIZE <- LEGEND_FONT_SIZE * 1.2
        # draw heatmap ---------------------------------------------------------
        HM <- Heatmap( MAT,
            # colors
            col = HEATMAP_COLOR_FUN, na_col = NA_COLOR,
            # borders
            border_gp = gpar(col = COLOR_BORDER, lwd = HeatmapBorderLineWidth), rect_gp = gpar(col = COLOR_BORDER, lwd = HeatmapCellLineWidth),
            # clustering
            cluster_rows = RowClustering, cluster_columns = ColumnClustering, show_row_dend = FALSE, show_column_dend = FALSE, 
            # row/column title
            row_title = ROW_TITLE, column_title = COLUMN_TITLE,
            # row labels              
            row_labels = ROW_LABEL, row_names_side = RowLabelsSide, row_names_rot = RowLabelsRotation, 
            row_names_gp = gpar(fontsize = RowLabelsSize, col = ROW_LABEL_COLOR ), row_names_centered = TRUE,
            # column labels 
            column_labels = COLUMN_LABEL, column_names_side = ColumnLabelsSide, column_names_rot = ColumnLabelsRotation,
            column_names_gp = gpar(fontsize = ColumnLabelsSize), column_names_centered = ColumnNamesCentered,  
            # add annotation
            top_annotation = TopAnnot, bottom_annotation = BottomAnnot, left_annotation = LeftAnnot, right_annotation = RightAnnot,
            # legend
            heatmap_legend_param = list(
                title         = LEGEND_TITLE, 
                title_gp      = gpar(fontsize=LEGEND_TITLE_FONT_SIZE),
                labels_gp     = gpar(fontsize=LEGEND_FONT_SIZE), 
                legend_width  = unit(LEGEND_WIDTH,"mm"),
                legend_height = unit(LEGEND_HEIGHT,"mm")
            ),
            # heatmap size 
            width  = WIDTH, height = HEIGHT

        )
        #-----------------------------------------------------------------------
        return(HM)
    }
#====================================================================================================================================================#

#===| GO ENRICHMENT ANALYSIS RESULT PLOTS |==========================================================================================================#
#---/ Top-ranked result dotplot or barplot /-------------------------------------------------------#
#
#--->> "use top-ranked result dotplot,barplot module in pathway plots"
#
#--------------------------------------------------------------------------------------------------#
#---/ create GO enrich top-ranked ontologies as network with ancestors /---------------------------#
#' @description top-ranked GO term ancestor tree network with expresssion profile summary
#' @param GoNetworDataList       REQUIRED. network node, edge data list. 
#' @param GoEnrichResult         REQUIRED. GO enrich analysis result table.
#' @param GoTerms                REQUIRED. GO terms info table.
#' @param IsAllDegsNetwork       use GO enrich result of 'ALL_DEGs'. default = FALSE
#' @param ScoreCutOff            score cut-off for node colors. if NULL, automatically calucalted from GO enrich result. default = NULL.
#' @param TopRank_N              top-rank #N of GO enrich analysis result to show on network.
#' @param NodeShape              node shape. default = 'circle'.
#' @param MarkRootNode           make unique colors for GO root node. default = TRUE
#' @param MarkTopRankNodes       add node stroke to GO enrich top-ranked nodes. default = TRUE
#' @param NetworkLayout          network layout. optimal=graphopt, large_net=lgl gem=gem, mid_net=dh, tree1=tree, tree2=layout_with_sugiyama. default = "optimal"
#' @param UseVisNetwork          draw network using 'visNetwork' package. default = FALSE
#' @param VisNetHighlightDegree  enable edge and node highlight option. work ONLY UseVisNetwork = TRUE. default = NULL
#' @export
    WTS_plot.GO.Top.Enrich.As.Network <- function( GoNetworDataList=NULL, GoEnrichResult=NULL, GoTerms=NULL,
        IsAllDegsNetwork=FALSE, ScoreCutOff=NULL, TopRank_N=10, 
        NodeShape="circle", MarkRootNode=TRUE, MarkTopRankNodes=TRUE, NetworkLayout="optimal",
        UseVisNetwork=FALSE, VisNetHighlightDegree=NULL
    )
    {
        suppressPackageStartupMessages(library("stringr"))
        suppressPackageStartupMessages(library("igraph"))
        suppressPackageStartupMessages(library("Hmisc"))
        suppressPackageStartupMessages(library("dplyr"))
        # pre-defined colors ---------------------------------------------------
        Color_UpOnly           <- "#dba9a7"
        Color_UpDown           <- "#fdc299"
        Color_DownUp           <- "#c0d9b2"
        Color_DownOnly         <- "#b5c6df"
        Color_AllDegs          <- "#ebcd7f"
        Color_DefaultNode      <- "#eaeaea"
        ColorFrame_UpOnly      <- "#B85450"
        ColorFrame_UpDown      <- "#fb8632"
        ColorFrame_DownUp      <- "#82B366"
        ColorFrame_DownOnly    <- "#6C8EBF"
        # scale down function --------------------------------------------------
        ScaleDown <- function( x, y )
        {
            cv <- (x-min(x))*((max(y)-min(y))/(max(x)-min(x))) + min(y)
            cv <- round(cv, 2)
            return(cv)
        }
        # create network -------------------------------------------------------
        if( is.null(GoNetworDataList) )
        {
            stop("|---!!! no GO network data list object found. REQUIRED. please check again. STOPPED.")
        }else{
            NET <- graph.data.frame(GoNetworDataList$edge[,c("subject","object")], directed=T)
        }
        # extract nodes from network -------------------------------------------
        NET_node <- get.data.frame(NET, 'vertices')
        # go terms check -------------------------------------------------------
        if( is.null(GoTerms) )
        {
            stop("|---!!! no GO Terms data found. REQUIRED. please check again. STOPPED.")
        }
        # set score cut off ----------------------------------------------------
        if( is.null(GoEnrichResult) )
        { 
            stop("|---!!! no GO enrich analysis result found. REQUIRED. please check again. STOPPED.")   
        }else{
            if( is.null(ScoreCutOff) )
            {
                if( IsAllDegsNetwork )
                {
                    score_cut_off <- GoEnrichResult %>% filter( input_geneset == "ALL_DEGs" ) %>% arrange(dplyr::desc(score)) %>% dplyr::slice(1:TopRank_N) %>% .$score %>% min()
                }else{
                    score_cut_off <- min(
                        GoEnrichResult %>% filter( input_geneset == "UP"   ) %>% arrange(dplyr::desc(score)) %>% dplyr::slice(1:TopRank_N) %>% .$score %>% min(),
                        GoEnrichResult %>% filter( input_geneset == "DOWN" ) %>% arrange(dplyr::desc(score)) %>% dplyr::slice(1:TopRank_N) %>% .$score %>% min()
                    )
                }
                score_cut_off <- score_cut_off/2
            }else{
                score_cut_off <- ScoreCutOff
            }
            # add enrich data (add color) ------------------------------------------
            if( IsAllDegsNetwork )
            {
                GoRes       <- GoEnrichResult %>% filter( input_geneset == "ALL_DEGs" ) %>% filter( ext_id %in% NET_node$name ) %>% dplyr::select(c("ext_id","score","enrich_genes")) 
                GoRes$color <- sapply(GoRes$score, function(sc) if( sc >=  score_cut_off ){ Color_AllDegs }else{ Color_DefaultNode })
                GoRes[which( GoRes$ext_id %in% (GoEnrichResult %>% filter( input_geneset == "ALL_DEGs" ) %>% arrange(dplyr::desc(score)) %>% dplyr::slice(1:TopRank_N) %>% .$ext_id)), "toprank" ] <- 1
                GoRes[which( is.na(GoRes$toprank)), "toprank"] <- 0
            }else{
                GoResUp   <- GoEnrichResult %>% filter( input_geneset == "UP"   ) %>% filter( ext_id %in% NET_node$name )
                GoResDown <- GoEnrichResult %>% filter( input_geneset == "DOWN" ) %>% filter( ext_id %in% NET_node$name )
                GoTermIDS <- union(GoResUp$ext_id, GoResDown$ext_id)
                GoRes <- data.frame(
                    ext_id      = GoTermIDS,
                    score_up    = GoResUp[match(GoTermIDS, GoResUp$ext_id), "score"],
                    score_down  = GoResDown[match(GoTermIDS, GoResDown$ext_id), "score"],
                    enrich_up   = GoResUp[match(GoTermIDS, GoResUp$ext_id), "enrich_genes"],
                    enrich_down = GoResDown[match(GoTermIDS, GoResDown$ext_id), "enrich_genes"]
                )
                toprank_profiles <- union(
                    GoEnrichResult %>% filter( input_geneset == "UP"   ) %>% arrange(dplyr::desc(score)) %>% dplyr::slice(1:TopRank_N) %>% .$ext_id,
                    GoEnrichResult %>% filter( input_geneset == "DOWN" ) %>% arrange(dplyr::desc(score)) %>% dplyr::slice(1:TopRank_N) %>% .$ext_id
                )
                GoRes[which(GoRes$ext_id %in% toprank_profiles), "toprank" ] <- 1
                GoRes[which( is.na(GoRes$toprank)),    "toprank"    ] <- 0
                GoRes[which(is.na(GoRes$score_up)),    "score_up"   ] <- 0
                GoRes[which(is.na(GoRes$score_down)),  "score_down" ] <- 0
                GoRes[which(is.na(GoRes$enrich_up)),   "enrich_up"  ] <- 0
                GoRes[which(is.na(GoRes$enrich_down)), "enrich_down"] <- 0
                # node color -------------------------------------------------------
                GoRes$color <- apply( GoRes[,c("score_up","score_down")], 1, function(sc)
                {
                    if( all(sc >= score_cut_off) )
                    { 
                        if( sc[1] >= sc[2] ){ clr <- Color_UpDown }else{ clr <- Color_DownUp }
                    }else{
                        if(       sc[1] >= score_cut_off & sc[2] <  score_cut_off ){ clr <- Color_UpOnly
                        }else if( sc[1]  < score_cut_off & sc[2] >= score_cut_off ){ clr <- Color_DownOnly
                        }else{ clr <- Color_DefaultNode }
                    }
                    return(clr)
                })
                # node frame color -------------------------------------------------
                GoRes$frame.color <- apply( GoRes[,c("score_up","score_down","toprank")], 1, function(sc)
                {
                    if( all(sc[1] >= score_cut_off & sc[2] >= score_cut_off) )
                    { 
                        if( sc[1] >= sc[2] ){ clr <- ColorFrame_UpDown }else{ clr <- ColorFrame_DownUp }
                    }else{
                        if(       sc[1] >= score_cut_off & sc[2] <  score_cut_off & sc[3] == 1 ){ clr <- ColorFrame_UpOnly
                        }else if( sc[1]  < score_cut_off & sc[2] >= score_cut_off & sc[3] == 1 ){ clr <- ColorFrame_DownOnly
                        }else{ clr <- "#FFFFFF" }
                    }
                    return(clr)
                })
            }
            GoRes$frame.width <- sapply( GoRes$toprank, function(tr) ifelse( tr == 1, 1.5, 0.5) )
        }
        # node properties-------------------------------------------------------
            NET_node$status       <- GoRes[match(NET_node$name, GoRes$ext_id), "color"]
            # node colors ------------------------------------------------------
            NET_node$color                           <- NET_node$status 
            NET_node[is.na(NET_node$color), "color"] <- Color_DefaultNode
            # node frame colors ------------------------------------------------
            NET_node$frame.color  <- GoRes[match(NET_node$name, GoRes$ext_id), "frame.color"]
            NET_node[is.na(NET_node$frame.color), "frame.color"] <- "#FFFFFF"
            # node frame thickness ---------------------------------------------
            NET_node$frame.width  <- GoRes[match(NET_node$name, GoRes$ext_id), "frame.width"]
            NET_node[is.na(NET_node$frame.width), "frame.width"] <- 1
            # node shape -------------------------------------------------------
            NET_node$shape <- NodeShape
            # node size --------------------------------------------------------
            NET_node$size <- sapply( NET_node$color, function(y) ifelse( y == Color_DefaultNode, 1, 2) )
            # node labels ----------------------------------------------------------
            NET_node$label        <- capitalize(GoTerms[match(NET_node$name, GoTerms$GO_TermID), "GO_Name"])
            NET_node$label.family <- "sans"
            NET_node$label.font   <- sapply( NET_node$color, function(y) ifelse( y == Color_DefaultNode, 1, 2) )
            NET_node$label.color  <- sapply( NET_node$color, function(y) ifelse( y == Color_DefaultNode, "#3d3d3d", "#0f0f0f") )
            NET_node$label.cex    <- 0.5
            # root nodes -------------------------------------------------------
            NET_node[which(NET_node$name %in% c("GO:0008150","GO:0005575","GO:0003674")), "size"]        <- 2.
            NET_node[which(NET_node$name %in% c("GO:0008150","GO:0005575","GO:0003674")), "color"]       <- "#a7c993"
            NET_node[which(NET_node$name %in% c("GO:0008150","GO:0005575","GO:0003674")), "frame.color"] <- "#FFFFFF"
            NET_node[which(NET_node$name %in% c("GO:0008150","GO:0005575","GO:0003674")), "label.color"] <- "#286a00"
            NET_node[which(NET_node$name %in% c("GO:0008150","GO:0005575","GO:0003674")), "label.font"]  <- 1
        # apply to network nodes -----------------------------------------------
        V(NET)$label        <- stringr::str_wrap(NET_node[match(V(NET)$name, NET_node$name), "label"], 20)
        V(NET)$color        <- NET_node[match(V(NET)$name, NET_node$name), "color"       ]
        V(NET)$frame.color  <- NET_node[match(V(NET)$name, NET_node$name), "frame.color" ]
        V(NET)$frame.width  <- NET_node[match(V(NET)$name, NET_node$name), "frame.width" ]
        V(NET)$shape        <- NET_node[match(V(NET)$name, NET_node$name), "shape"       ]
        V(NET)$label.family <- NET_node[match(V(NET)$name, NET_node$name), "label.family"]
        V(NET)$label.font   <- NET_node[match(V(NET)$name, NET_node$name), "label.font"  ]
        V(NET)$label.color  <- NET_node[match(V(NET)$name, NET_node$name), "label.color" ]
        V(NET)$size         <- NET_node[match(V(NET)$name, NET_node$name), "size"        ]
        V(NET)$label.cex    <- NET_node[match(V(NET)$name, NET_node$name), "label.cex"   ]
        # edge properties ------------------------------------------------------
        E(NET)$arrow.size  <- 0.5
        E(NET)$arrow.width <- 0.5
        E(NET)$lty         <- 1
        E(NET)$color       <- "#cccccc"
        # network layout -------------------------------------------------------
        if(       NetworkLayout == "optimal"   ){ 
            LAYOUT       <- layout_with_graphopt(NET, niter=1500, charge=0.01)
            VisNetLayout <- "layout_with_graphopt" 
        }else if( NetworkLayout == "large_net" ){ 
            LAYOUT       <- layout_with_lgl(NET, repulserad=vcount(NET)^2 )
            VisNetLayout <- "layout_with_lgl" 
        }else if( NetworkLayout == "gem"       ){ 
            LAYOUT       <- layout_with_gem(NET, maxiter = 100*vcount(NET)^2, temp.min=1/20 ) 
            VisNetLayout <- "layout_with_gem"
        }else if( NetworkLayout == "mid_net"   ){ 
            LAYOUT <- layout_with_dh(NET, maxiter=50, fineiter=100)            
            VisNetLayout <- "layout_with_dh"
        }else if( NetworkLayout == "tree1"     ){ 
            LAYOUT       <- layout_as_tree(NET, mode="all")
            V(NET)$label <- stringr::str_wrap(NET_node[match(V(NET)$name, NET_node$name), "label"], 10)
            VisNetLayout <- "layout_as_tree"
        }else if( NetworkLayout == "tree2"     ){ 
            LAYOUT       <- layout_with_sugiyama(NET, maxiter=500, vgap=4, hgap=10)
            V(NET)$label <- stringr::str_wrap(NET_node[match(V(NET)$name, NET_node$name), "label"], 10)
            VisNetLayout <- "layout_with_sugiyama"
        }     
        # make network object --------------------------------------------------
        if( UseVisNetwork )
        {
            suppressPackageStartupMessages(library("visNetwork"))
            # visnetwork nodes -------------------------------------------------
            VISNET_node <- get.data.frame(NET, 'vertices')
            VISNET_node <- VISNET_node %>% dplyr::rename( id=name, value=size ) %>% mutate( label = gsub("\\\n", " ", label) )
            VISNET_node$shape <- "dot"
            # visnetwork edges -------------------------------------------------
            VISNET_edge <- get.data.frame(NET, 'edges')
            VISNET_edge$arrows <- "to"
            VISNET_edge$shadow <- FALSE 
            # draw network -----------------------------------------------------
            NETWORK <- visNetwork( VISNET_node, VISNET_edge ) %>% visIgraphLayout( layout = VisNetLayout )
            # hightlight edges option ------------------------------------------
            if( all(!is.null(VisNetHighlightDegree) & is.numeric(VisNetHighlightDegree)) )
            { NETWORK <- NETWORK %>% visOptions(highlightNearest = list(enabled = T, degree = VisNetHighlightDegree, hover = T, labelOnly=F)) }
        }else{
            # draw network -----------------------------------------------------
            NETWORK <- plot( NET , vertex.label.degree = pi/2, layout = LAYOUT, asp=0)
        }
        return(NETWORK)
    }
#--------------------------------------------------------------------------------------------------#
#====================================================================================================================================================#
   

