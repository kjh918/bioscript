

setwd("/storage/home/kangsm/analysis/TMB")

D1 <- openxlsx::read.xlsx("TMB_Harmoniation_Project_PhaseII_Data.xlsx", sheet=2)
D2 <- openxlsx::read.xlsx("TMB_Harmoniation_Project_PhaseII_Data.xlsx", sheet=3)
#
FMI_fit <- lm( FMI ~ WES.TMB, data=D1)
TSO_fit <- lm( ILLUM ~ WES.TMB, data=D1)
#
gp_clinicl_all <- ggplot(D1) +
    geom_point(aes(x = WES.TMB, y = FMI), color = "#9f3122", size=4, alpha=0.5) +
    geom_point(aes(x = WES.TMB, y = ILLUM), color = "#1a5387", size=4, alpha=0.5) +
    coord_cartesian(xlim = c(0,280), ylim = c(0,280)) +
    geom_segment(x=0, xend=280, y=0, yend=280, linetype=2, linewidth=1) +
    geom_abline( intercept = TSO_fit$coefficients[1], slope = TSO_fit$coefficients[2], col = "#1a5387") +
    geom_abline( intercept = FMI_fit$coefficients[1], slope = FMI_fit$coefficients[2], col = "#9f3122") +
        theme(
        axis.line.x  = element_line(colour = "#6f6d6d"), 
        axis.line.y  = element_line(colour = "#6f6d6d"),
        panel.background = element_rect(fill='white'), 
        panel.grid.major = element_line(color = '#d1d1d1', linetype = 'dotted'), 
        legend.position  = "none",
        plot.margin  = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
    ) + labs(x="WES TMB", y="Panel TMB")

gp_clinicl_100 <- ggplot(D1) +
    geom_point(aes(x = WES.TMB, y = FMI), color = "#9f3122", size=4, alpha=0.5) +
    geom_point(aes(x = WES.TMB, y = ILLUM), color = "#1a5387", size=4, alpha=0.5) +
    coord_cartesian(xlim = c(0,100), ylim = c(0,100)) +
    geom_segment(x=0, xend=100, y=0, yend=100, linetype=2, linewidth=1) +
    geom_abline( intercept = TSO_fit$coefficients[1], slope = TSO_fit$coefficients[2], col = "#1a5387") +
    geom_abline( intercept = FMI_fit$coefficients[1], slope = FMI_fit$coefficients[2], col = "#9f3122") +
        theme(
        axis.line.x  = element_line(colour = "#6f6d6d"), 
        axis.line.y  = element_line(colour = "#6f6d6d"),
        panel.background = element_rect(fill='white'), 
        panel.grid.major = element_line(color = '#d1d1d1', linetype = 'dotted'), 
        legend.position  = "none",
        plot.margin  = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
    ) + labs(x="WES TMB", y="Panel TMB")
#
FMI_fit2 <- lm( FMI ~ WES.TMB, data=D2)
TSO_fit2 <- lm( ILLUM ~ WES.TMB, data=D2)
#
gp2_clinicl_all <- ggplot(D2) +
    geom_point(aes(x = WES.TMB, y = FMI), color = "#9f3122", size=4, alpha=0.5) +
    geom_point(aes(x = WES.TMB, y = ILLUM), color = "#1a5387", size=4, alpha=0.5) +
    coord_cartesian(xlim = c(0,280), ylim = c(0,280)) +
    geom_segment(x=0, xend=280, y=0, yend=280, linetype=2, linewidth=1) +
    geom_abline( intercept = TSO_fit2$coefficients[1], slope = TSO_fit2$coefficients[2], col = "#1a5387") +
    geom_abline( intercept = FMI_fit2$coefficients[1], slope = FMI_fit2$coefficients[2], col = "#9f3122") +
        theme(
        axis.line.x  = element_line(colour = "#6f6d6d"), 
        axis.line.y  = element_line(colour = "#6f6d6d"),
        panel.background = element_rect(fill='white'), 
        panel.grid.major = element_line(color = '#d1d1d1', linetype = 'dotted'), 
        legend.position  = "none",
        plot.margin  = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
    ) + labs(x="WES TMB", y="Panel TMB")

gp2_clinicl_100 <- ggplot(D2) +
    geom_point(aes(x = WES.TMB, y = FMI), color = "#9f3122", size=4, alpha=0.5) +
    geom_point(aes(x = WES.TMB, y = ILLUM), color = "#1a5387", size=4, alpha=0.5) +
    coord_cartesian(xlim = c(0,100), ylim = c(0,100)) +
    geom_segment(x=0, xend=100, y=0, yend=100, linetype=2, linewidth=1) +
    geom_abline( intercept = TSO_fit2$coefficients[1], slope = TSO_fit2$coefficients[2], col = "#1a5387") +
    geom_abline( intercept = FMI_fit2$coefficients[1], slope = FMI_fit2$coefficients[2], col = "#9f3122") +
        theme(
        axis.line.x  = element_line(colour = "#6f6d6d"), 
        axis.line.y  = element_line(colour = "#6f6d6d"),
        panel.background = element_rect(fill='white'), 
        panel.grid.major = element_line(color = '#d1d1d1', linetype = 'dotted'), 
        legend.position  = "none",
        plot.margin  = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
    ) + labs(x="WES TMB", y="Panel TMB")

#
D3 = rbind(
    data.frame( sampleid=D1$TMB.ID, WES.TMB=D1$WES.TMB, FMI=D1$FMI, TSO=D1$ILLUM ),
    data.frame( sampleid=D2$TMB.ID, WES.TMB=D2$WES.TMB, FMI=D2$FMI, TSO=D2$ILLUM )
)
#
FMI_fit3 <- lm( FMI ~ WES.TMB, data=D3)
TSO_fit3 <- lm( TSO ~ WES.TMB, data=D3)
#
gp3_clinicl_all <- ggplot(D3) +
    geom_point(aes(x = WES.TMB, y = FMI), color = "#9f3122", size=4, alpha=0.5) +
    geom_point(aes(x = WES.TMB, y = TSO), color = "#1a5387", size=4, alpha=0.5) +
    coord_cartesian(xlim = c(0,280), ylim = c(0,280)) +
    geom_segment(x=0, xend=280, y=0, yend=280, linetype=2, linewidth=1) +
    geom_abline( intercept = TSO_fit3$coefficients[1], slope = TSO_fit3$coefficients[2], col = "#1a5387") +
    geom_abline( intercept = FMI_fit3$coefficients[1], slope = FMI_fit3$coefficients[2], col = "#9f3122") +
        theme(
        axis.line.x  = element_line(colour = "#6f6d6d"), 
        axis.line.y  = element_line(colour = "#6f6d6d"),
        panel.background = element_rect(fill='white'), 
        panel.grid.major = element_line(color = '#d1d1d1', linetype = 'dotted'), 
        legend.position  = "none",
        plot.margin  = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
    ) + labs(x="WES TMB", y="Panel TMB")

gp3_clinicl_100 <- ggplot(D3) +
    geom_point(aes(x = WES.TMB, y = FMI), color = "#9f3122", size=4, alpha=0.5) +
    geom_point(aes(x = WES.TMB, y = TSO), color = "#1a5387", size=4, alpha=0.5) +
    coord_cartesian(xlim = c(0,100), ylim = c(0,100)) +
    geom_segment(x=0, xend=100, y=0, yend=100, linetype=2, linewidth=1) +
    geom_abline( intercept = TSO_fit3$coefficients[1], slope = TSO_fit3$coefficients[2], col = "#1a5387") +
    geom_abline( intercept = FMI_fit3$coefficients[1], slope = FMI_fit3$coefficients[2], col = "#9f3122") +
        theme(
        axis.line.x  = element_line(colour = "#6f6d6d"), 
        axis.line.y  = element_line(colour = "#6f6d6d"),
        panel.background = element_rect(fill='white'), 
        panel.grid.major = element_line(color = '#d1d1d1', linetype = 'dotted'), 
        legend.position  = "none",
        plot.margin  = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
    ) + labs(x="WES TMB", y="Panel TMB")

#
cowplot::plot_grid( gp_clinicl_all, gp2_clinicl_all, gp3_clinicl_all, gp_clinicl_100,  gp2_clinicl_100,  gp3_clinicl_100, nrow=2 )


Get.WES.TMB <- function( panel_tmb, fit_slope, fit_y ){ (panel_tmb - fit_y)/fit_slope }


Get.WES.TMB( 10, FMI_fit$coefficients[2],  FMI_fit$coefficients[1] )
Get.WES.TMB( 10, FMI_fit2$coefficients[2], FMI_fit2$coefficients[1] )
Get.WES.TMB( 10, FMI_fit3$coefficients[2], FMI_fit3$coefficients[1] )

Get.WES.TMB( 10, TSO_fit$coefficients[2],  TSO_fit$coefficients[1] )
Get.WES.TMB( 10, TSO_fit2$coefficients[2], TSO_fit2$coefficients[1] )
Get.WES.TMB( 10, TSO_fit3$coefficients[2], TSO_fit3$coefficients[1] )

10 * FMI_fit$coefficients[2] + FMI_fit$coefficients[1]
10 * FMI_fit2$coefficients[2] + FMI_fit2$coefficients[1]
10 * FMI_fit3$coefficients[2] + FMI_fit3$coefficients[1]

