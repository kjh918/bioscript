

RUN_SCRIPT=/storage/home/kangsm/runScripts
BASE_DIR=/data/wes

SEQ_FOLDER="23WES03"
SAMPLE_LIST=$(cat ${BASE_DIR}/${SEQ_FOLDER}/meta/sample.list)

for SEQ_ID in $SAMPLE_LIST
do
    BAM_DIR=${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/bam
    cp ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.realign.recal.bam ${BAM_DIR}/${SEQ_ID}.analysisReady.bam
    cp ${BAM_DIR}/${SEQ_ID}.merged.dup.marked.realign.recal.bai ${BAM_DIR}/${SEQ_ID}.analysisReady.bam.bai
    echo ${SEQ_ID}
done

for SEQ_ID in $SAMPLE_LIST
do
    LOG_DIR=${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/log
    TMB_DIR=${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/tmb
    rm -r ${TMB_DIR}
    rm -r ${LOG_DIR}
    mkdir -p ${LOG_DIR}   
     mkdir -p ${TMB_DIR} 
    echo ${SEQ_ID}
done

for SEQ_ID in $SAMPLE_LIST
do
    COV_FILE=${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/qcfiles/${SEQ_ID}.alfred.target.coverage.txt
    BED_FILE=/storage/references_and_index/hg19/bed/Twist_WES_2.0/GRCh37_twist.exome.2.0.target.bed
    MAF_FILE=${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/vcf/${SEQ_ID}.mutect2.filtered.vep.annotated.maf

    # Rscript ${RUN_SCRIPT}/WES_run.TMB_Calculation.R \
    #     --seqFolder ${SEQ_FOLDER} \
    #     --seqID ${SEQ_ID} \
    #     --baseDir /data/wes \
    #     --cov ${COV_FILE} \
    #     --bed ${BED_FILE} \
    #     --maf ${MAF_FILE} \
    #     --useReference TRUE \
    #     --useReferenceSize FALSE \
    #     --DP 10 \
    #     --variantCallMode "tonly"

    TMB_DIR=${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/tmb
    # cp ${TMB_DIR}/${SEQ_ID}.TMB.calculation.result.txt ${TMB_DIR}/${SEQ_ID}.TMB.GCX2.method.res.txt

    Rscript ${RUN_SCRIPT}/WES_run.TMB_Calculation.R \
        --seqFolder ${SEQ_FOLDER} \
        --seqID ${SEQ_ID} \
        --baseDir /data/wes \
        --cov ${COV_FILE} \
        --bed ${BED_FILE} \
        --maf ${MAF_FILE} \
        --useReference FALSE \
        --useReferenceSize FALSE \
        --DP 10 --tumorALT 2 --vafSNV 0.02 --vafINDEL 0.05 --popAF 0.01 \
        --excludeSilent FALSE \
        --excludeClinvar FALSE \
        --excludeAlphaMissense FALSE \
        --variantCallMode "tonly"

    cp ${TMB_DIR}/${SEQ_ID}.TMB.calculation.result.txt ${TMB_DIR}/${SEQ_ID}.TMB.GCX2.method.res.txt
    
    echo ${SEQ_ID}
done


TMB_REF = TMB_GCX = TMB_REF2 = TMB_GCX2 = data.frame()

SEQ_FOLDER = "23WES03"
SAMPLE_LIST = read.delim(sprintf("/data/wes/%s/meta/sample.list",SEQ_FOLDER), header=F)[,1]

for( SEQ_ID in SAMPLE_LIST )
{
    if( SEQ_ID != "23WES04_30" )
    {
        res_dir = sprintf("/data/wes/%s/%s/tmb", SEQ_FOLDER, SEQ_ID )
        #tmb.ref = read.delim(sprintf("%s/%s.TMB.REF.method.res.txt", res_dir, SEQ_ID))
        tmb.ref2 = read.delim(sprintf("%s/%s.TMB.REF2.method.res.txt", res_dir, SEQ_ID))
        #tmb.gcx = read.delim(sprintf("%s/%s.TMB.GCX.method.res.txt", res_dir, SEQ_ID))
        #tmb.gcx2 = read.delim(sprintf("%s/%s.TMB.GCX2.method.res.txt", res_dir, SEQ_ID))

        #TMB_REF = rbind(TMB_REF, tmb.ref )
        TMB_REF2 = rbind(TMB_REF2, tmb.ref2 )
        #TMB_GCX = rbind(TMB_GCX, tmb.gcx )
        #TMB_GCX2 = rbind(TMB_GCX2, tmb.gcx2 )
    }
}

TB = data.frame(
    TMB_REF[,1:2],
    REF1=TMB_REF$TMB_value,
    REF2=TMB_REF2$TMB_value,
    GCX1=TMB_GCX$TMB_value,
    GCX2=TMB_GCX2$TMB_value,
    REF1.var=TMB_REF$variants,
    REF2.var=TMB_REF2$variants,
    GCX1.var=TMB_GCX$variants,
    GCX2.var=TMB_GCX2$variants
)

lm( REF2 ~ REF1+0, data=TB )

lm( GCX ~ REF1+0, data=TB )

plot(x=TB$REF1, y=TB$REF1, xlim=c(0,170), ylim=c(0,170),
    xlab="TMB Harmonization Project TMB (Mut/Mb)",
    ylab="TMB (other methods, Mut/Mb)"

)
points(x=TB$REF1, y=TB$REF2, col="red")
points(x=TB$REF1, y=TB$GCX1, col="blue")
points(x=TB$REF1, y=TB$GCX2, col="green")

abline(lm( REF2 ~ REF1+0, data=TB ), col="red")
abline(lm( GCX1 ~ REF1+0, data=TB ), col="blue")


plot(x=TB$REF1, y=TB$REF1, xlim=c(0,170), ylim=c(0,170),
    xlab="TMB Harmonization Project TMB (Mut/Mb)",
    ylab="TMB (other methods, Mut/Mb)"

)
points(x=TB$REF1, y=TB$REF2, col="red")
points(x=TB$REF1, y=TB$GCX1, col="blue")
points(x=TB$REF1, y=TB$GCX2, col="green")

abline(lm( REF2 ~ REF1+0, data=TB ), col="red")
abline(lm( GCX1 ~ REF1+0, data=TB ), col="blue")

lm.REF2 = lm( REF2 ~ REF1, data =TB)
lm.GCX1 = lm( GCX1 ~ REF1, data =TB)
lm.GCX2 = lm( GCX2 ~ REF1, data =TB)



tmb1 = ggplot(TB) +
    geom_point(aes(x=REF1, y=REF2), col="#474f1b", size=4, alpha=0.5) +
    geom_point(aes(x=REF1, y=GCX1), col="#1a5387", size=4, alpha=0.5) +
    geom_point(aes(x=REF1, y=GCX2), col="#9f3122", size=4, alpha=0.5) +
    coord_cartesian(xlim = c(0, 170), ylim = c(0,170)) +
    geom_segment(x=0, xend=170, y=0, yend=170, linetype=2, linewidth=1) +
    theme(
        axis.line.x  = element_line(colour = "#6f6d6d"), 
        axis.line.y  = element_line(colour = "#6f6d6d"),
        panel.background = element_rect(fill='white'), 
        panel.grid.major = element_line(color = '#d1d1d1', linetype = 'dotted'), 
        legend.position  = "none",
        plot.margin  = unit(c(0.2, 0.2, 0.2, 0.2), "cm")

    ) +
    labs( x = "TMB Harmonization Project TMB (Mut/Mb)", y = "TMB (other methods, Mut/Mb)") +
    geom_abline( intercept = lm.GCX1$coefficients[1], slope = lm.GCX1$coefficients[2], col = "#1a5387") +
    geom_abline( intercept = lm.GCX2$coefficients[1], slope = lm.GCX2$coefficients[2], col = "#9f3122") +
    geom_abline( intercept = lm.REF2$coefficients[1], slope = lm.REF2$coefficients[2], col = "#474f1b")     

tmb2 = ggplot(TB) +
    geom_point(aes(x=REF1, y=REF2), col="#474f1b", size=4, alpha=0.5) +
    geom_point(aes(x=REF1, y=GCX1), col="#1a5387", size=4, alpha=0.5) +
    geom_point(aes(x=REF1, y=GCX2), col="#9f3122", size=4, alpha=0.5) +
    coord_cartesian(xlim = c(0, 45), ylim = c(0,45)) +
    geom_segment(x=0, xend=45, y=0, yend=45, linetype=2, linewidth=1) +
    theme(
        axis.line.x  = element_line(colour = "#6f6d6d"), 
        axis.line.y  = element_line(colour = "#6f6d6d"),
        panel.background = element_rect(fill='white'), 
        panel.grid.major = element_line(color = '#d1d1d1', linetype = 'dotted'), 
        legend.position  = "none",
        plot.margin  = unit(c(0.2, 0.2, 0.2, 0.2), "cm")

    ) +
    labs( x = "TMB Harmonization Project TMB (Mut/Mb)", y = "TMB (other methods, Mut/Mb)") +
    geom_abline( intercept = lm.GCX1$coefficients[1], slope = lm.GCX1$coefficients[2], col = "#1a5387") +
    geom_abline( intercept = lm.GCX2$coefficients[1], slope = lm.GCX2$coefficients[2], col = "#9f3122") +
    geom_abline( intercept = lm.REF2$coefficients[1], slope = lm.REF2$coefficients[2], col = "#474f1b")     

cowplot::plot_grid(tmb1, tmb2, nrow=1)


CHK= rbind(
    data.frame(TB[,c(1,2)], group="REF_1", variants=TB$REF1.var),
    data.frame(TB[,c(1,2)], group="REF_2", variants=TB$REF2.var),
    data.frame(TB[,c(1,2)], group="GCX_1", variants=TB$GCX1.var),
    data.frame(TB[,c(1,2)], group="GCX_2", variants=TB$GCX2.var)
)

TB <- TB %>% left_join(., sinfo[,c("seq_id","sample_group","sample_tissue","sample_info")], by="seq_id")

TB2 = TB[which(TB$sample_tissue %in% c("liver","lung","colorectal","prostate")),]

TB2 <- TB2 %>% dplyr::rename(Cancer=sample_tissue, TMB=REF1)

ggplot(TB2, aes(x=sample_tissue, y=REF1, color=sample_tissue))


gp_tmb = ggplot(TB2, aes(x=Cancer, y=TMB, fill=Cancer)) + 
    geom_boxplot(outlier.shape = NA, alpha=0.7) +
    geom_jitter(width=0.2, height=0.35, size=2, color = "#2c2b2b") +  
    geom_hline(yintercept=10, linetype="dashed", color="#0a80db") +
    labs(x = "Cancer Type", y = "Tumor Mutational Burden (TMB)") + 
    coord_cartesian(ylim = c(0,50) ) +
    theme(
        axis.line.x  = element_line(colour = "#6f6d6d"), 
        axis.line.y  = element_line(colour = "#6f6d6d"),
        axis.text.x  = element_text(face = "bold", size = 15),
        axis.text.y  = element_text(face = "bold", size = 15),
        axis.title.x = element_text(face = "bold", size = 20, margin = margin(t = 7)),
        axis.title.y = element_text(face = "bold", size = 20, margin = margin(r = 7)),
        panel.background = element_rect(fill='#f8f7f7'), 
        panel.grid.major = element_line(color = "#dcdbdb", linetype = 'dotted'),
        panel.grid.minor = element_line(color = "#dcdbdb", linetype = 'dotted'), 
        legend.position  = "none",
        plot.title=element_text(face="bold", size=14),
        plot.margin  = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    ) 



PS <- TB2 %>% filter( sample_group == "VOS_CRC_23_0016") 
PS$seq_id = c("TISSUE","OG-p2","OG-p4","OG-p8","OG-p10.1","OG-p10.2","OG-p10.3")
PS$seq_id = factor(PS$seq_id, levels=c("TISSUE","OG-p2","OG-p4","OG-p8","OG-p10.1","OG-p10.2","OG-p10.3"))
PS$group = c("TISSUE","OG-P2","OG-P4","OG-P8","OG-P10","OG-P10","OG-P10")
PS$group = factor(PS$group, levels=c("TISSUE","OG-P2","OG-P4","OG-P8","OG-P10"))

gp_tmb2 = ggplot(PS, aes(x=group, y=TMB, group=sample_group)) + 
    geom_hline(yintercept=10.03, linetype="dashed",color="#0a80db", alpha=0.5) +
    geom_line(linewidth=0.5, color="#bfbfbf") + 
    geom_point(aes(col=sample_origin), size=7, alpha=0.9) + 
    coord_cartesian(ylim = c(0,15) ) +
    labs(x = "Organoid Passages", y = "Tumor Mutational Burden (TMB)") + 
    theme(
        axis.line.x  = element_line(colour = "#6f6d6d"), 
        axis.line.y  = element_line(colour = "#6f6d6d"),
        axis.text.x  = element_text(face = "bold", size = 15),
        axis.text.y  = element_text(face = "bold", size = 15),
        axis.title.x = element_text(face = "bold", size = 17, margin = margin(t = 7)),
        axis.title.y = element_text(face = "bold", size = 20, margin = margin(r = 7)),
        panel.background = element_rect(fill='#f8f7f7'), 
        panel.grid.major = element_line(color = "#dcdbdb", linetype = 'dotted'),
        panel.grid.minor = element_line(color = "#dcdbdb", linetype = 'dotted'), 
        legend.position  = "none",
        plot.title=element_text(face="bold", size=14),
        plot.margin  = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    ) +
    scale_color_manual(values=c("TS"="#2055DA", "ORG"="#a72e2e")) 

ggsave(gp_tmb2, file = "TMB_plot_example2.png", width=7, height=6, unit='in', dpi=150, type='cairo')





