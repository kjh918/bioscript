setwd("/storage/home/kangsm/analysis/cnv_seqc2/analysis_hg38")

load("concordance_validation_by_params.Rdata")

SEQC2_COLOR = c("FD1_T"="#D5E8D4", "FD2_T"="#82B366", "FD3_T"="#2D7600", "NV1_T"="#DAE8FC", "NV2_T"="#6C8EBF", "NV3_T"="#1BA1E2")

DATA = pcn_Tonly_valid
TAG  = "PureCN_Tonly"

CN_plot <- ggplot(DATA, aes( x = id, y = concordance_CN, fill = id)) + 
    geom_bar( width = 0.75, stat = "identity", position = "dodge") +
    geom_hline( yintercept = 0, color = "#6f6d6d" ) + 
    coord_cartesian( ylim = c(0,100) ) + 
    labs( x = "", y = "Concordance (%)", fill = "SEQC2 ID" ) +
    theme(
        panel.background = element_rect(fill='white'), 
        panel.grid.major = element_line(color = '#dcdbdb', linetype = 'dotted', linewidth = 0.5),
        plot.margin  = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        plot.title   = element_text(face = "bold", size = 12),
        axis.line.y  = element_line(colour = "#6f6d6d"),
        axis.text.y  = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 12, margin = margin(t = 5)),
        axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 5)),
        axis.ticks.x = element_blank(), 
        legend.key.size = unit(3, 'mm')
    ) + 
    scale_fill_manual(values=SEQC2_COLOR) +
    ggtitle( "CopyNumber" )

ABS_CNA_plot <- ggplot(DATA, aes( x = id, y = concordance_AbsoluteCNA, fill = id)) + 
    geom_bar( width = 0.75, stat = "identity", position = "dodge") +
    geom_hline( yintercept = 0, color = "#6f6d6d" ) + 
    coord_cartesian( ylim = c(0,100) ) + 
    labs( x = "", y = "Concordance (%)", fill = "SEQC2 ID" ) +
    theme(
        panel.background = element_rect(fill='white'), 
        panel.grid.major = element_line(color = '#dcdbdb', linetype = 'dotted', linewidth = 0.5),
        plot.margin  = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        plot.title   = element_text(face = "bold", size = 12),
        axis.line.y  = element_line(colour = "#6f6d6d"),
        axis.text.y  = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 12, margin = margin(t = 5)),
        axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 5)),
        axis.ticks.x = element_blank(), 
        legend.key.size = unit(3, 'mm')
    ) + 
    scale_fill_manual(values=SEQC2_COLOR) +
    ggtitle( "AbsoluteCN CNA" )
PLD_CNA_plot <- ggplot(DATA, aes( x = id, y = concordance_PloidyCNA, fill = id)) + 
    geom_bar( width = 0.75, stat = "identity", position = "dodge") +
    geom_hline( yintercept = 0, color = "#6f6d6d" ) + 
    coord_cartesian( ylim = c(0,100) ) + 
    labs( x = "", y = "Concordance (%)", fill = "SEQC2 ID" ) +
    theme(
        panel.background = element_rect(fill='white'), 
        panel.grid.major = element_line(color = '#dcdbdb', linetype = 'dotted', linewidth = 0.5),
        plot.margin  = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        plot.title   = element_text(face = "bold", size = 12),
        axis.line.y  = element_line(colour = "#6f6d6d"),
        axis.text.y  = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 12, margin = margin(t = 5)),
        axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 5)),
        axis.ticks.x = element_blank(), 
        legend.key.size = unit(3, 'mm')
    ) + 
    scale_fill_manual(values=SEQC2_COLOR) +
    ggtitle( "PloidyCN CNA" )

conPlot <- cowplot::plot_grid( CN_plot, ABS_CNA_plot, PLD_CNA_plot, ncol=1)
ggsave(conPlot, width=5, height=12, unit='in', dpi=150, type='cairo', file = sprintf("%s.params.concordance.png", TAG )  )
    


