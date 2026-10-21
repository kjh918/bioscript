
options(stringsAsFactors=FALSE)   
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("RMySQL"))
suppressPackageStartupMessages(library("Hmisc"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("dplyr"))
#----------------------------------------------------------------------------------------------------------------------------------------------------# 
    option_list = list(
        make_option(c("--seqFolder")          , action="store", default=NA,   type="character", help="$SEQ_FOLDER" ),
        make_option(c("--seqID")              , action="store", default=NA,   type="character", help="$SEQ_ID" ),
        make_option(c("--baseDir")            , action="store", default=NA,   type="character", help="Base Work Directory" ),
        make_option(c("--cov")                , action="store", default=NA,   type="character", help="Alfred QC Coverage File" ),
        make_option(c("--bed")                , action="store", default=NA,   type="character", help="Target BED file" ),
        make_option(c("--maf")                , action="store", default=NA,   type="character", help="MAF file" ),
        make_option(c("--RD")                 , action="store", default=10,   type="integer",   help="Average Read Depth Cut-Off for Target Region Filtering" ), 
        make_option(c("--DP")                 , action="store", default=10,   type="integer",   help="Variant Cut-Off: Tumor Read Depth" ), 
        make_option(c("--tumorALT")           , action="store", default=2,    type="integer",   help="Variant Cut-Off: Minimum Tumor-Allele Count" ), 
        make_option(c("--snvVAF")             , action="store", default=0.02, type="double",    help="Variant Cut-Off: SNV VAF" ), 
        make_option(c("--indelVAF")           , action="store", default=0.05, type="double",    help="Variant Cut-Off: INDEL VAF" ),
        make_option(c("--popAF")              , action="store", default=0.01, type="double",    help="Variant Cut-Off: Population AF" ),
        make_option(c("--filterSynonymous")   , action="store", default=NA,   type="logical",   help="exclude synonymous variants or not" ), 
        make_option(c("--filterClinvar")      , action="store", default=NA,   type="logical",   help="exclude clinvar's pathogenic variants or not" ), 
        make_option(c("--filterAlphaMissense"), action="store", default=NA,   type="logical",   help="exclude alpha-missense's pathogenic variants or not" ), 
        make_option(c("--variantCallMode")    , action="store", default=NA,   type="character", help="variant calling mode" ), 
    )
    #--------------------------------------------------------------------------#
    ARGS = parse_args(OptionParser(option_list=option_list))
    #--------------------------------------------------------------------------#
    seqFolder                        <- ARGS$seqFolder
    seqID                            <- ARGS$seqID
    BASE_DIR                         <- ARGS$baseDir 
    COVERAGE_FILE                    <- ARGS$cov
    TARGET_BED_FILE                  <- ARGS$bed
    MAF_FILE                         <- ARGS$maf
    cutoff_DP                        <- ARGS$RD
    cutoff_TUMOR_DP                  <- ARGS$DP
    cutoff_TUMOR_ALT                 <- ARGS$tumorALT
    cutoff_SNV_VAF                   <- ARGS$snvVAF
    cutoff_INDEL_VAF                 <- ARGS$indelVAF
    cutoff_POP_AF                    <- ARGS$popAF
    FILTER_SYNONYMOUS                <- ARGS$filterSynonymous
    FILTER_CLINVAR_PATHOGENIC        <- ARGS$filterClinvar
    FILTER_ALPHA_MISSENSE_PATHOGENIC <- ARGS$filterAlphaMissense
    VARIANT_CALL_MODE                <- ARGS$variantCallMode
#----------------------------------------------------------------------------------------------------------------------------------------------------#
    TMB_RES_DIR <- sprintf("%s/%s/%s/tmb", BASE_DIR, seqFolder, seqID)
    if( ! dir.exists(TMB_RES_DIR) ) { dir.create(TMB_RES_DIR) }
#----------------------------------------------------------------------------------------------------------------------------------------------------#
    COV <- read.delim(COVERAGE_FILE)
    #----------------------------------------------------------------------------------------------#
    BED           <- read.delim(TARGET_BED_FILE, header=F)
    colnames(BED) <- c("chr","start","end","gene_raw")
    BED$gene      <- sapply(BED$gene_raw, function(gr) unlist(strsplit(gr, ";"))[1])
    BED           <- BED[,c("chr","start","end","gene")]
    #----------------------------------------------------------------------------------------------#
    COV <- COV %>% mutate( key = paste(Chr,Start,End,sep="_"))
    BED <- BED %>% mutate( key = paste(chr,start,end,sep="_"))
    #----------------------------------------------------------------------------------------------#
    COV <- COV %>% left_join(., BED[,c("key","gene")], by="key")
#----------------------------------------------------------------------------------------------------------------------------------------------------#
    covRes <- COV %>% filter( gene != "", AvgCov >= cutoff_DP ) %>%
        group_by(Chr, Start ) %>%
        summarise( End = max(End), AvgCov = mean(AvgCov) ) %>%
        mutate( LibLength =  End - Start )
#----------------------------------------------------------------------------------------------------------------------------------------------------#
    LIB_LENGTH <- sum(covRes$LibLength)                                                             # ExomeRegion Length (bp)
#----------------------------------------------------------------------------------------------------------------------------------------------------#
    ExceptVariantClass <- c("3'Flank","3'UTR","5'Flank","5'UTR","Intron","RNA","IGR","Splice_Region")
#----------------------------------------------------------------------------------------------------------------------------------------------------#
    MAF <- read.delim(MAF_FILE)
    #----------------------------------------------------------------------------------------------#
    mf <- MAF %>% 
        filter( Variant_Classification %nin% ExceptVariantClass ) %>%
        mutate( VAF = round(TUMOR_ALT_COUNT/(TUMOR_REF_COUNT+TUMOR_ALT_COUNT), 4) )
    #----------------------------------------------------------------------------------------------#
    mf.tmb <- mf
#----------------------------------------------------------------------------------------------------------------------------------------------------#
    if( FILTER_SYNONYMOUS )
    { mf.tmb <- mf.tmb %>% filter( Variant_Classification != "Silent" ) }
#----------------------------------------------------------------------------------------------------------------------------------------------------#
    if( FILTER_CLINVAR_PATHOGENIC )
    { mf.tmb <- mf.tmb %>% filter( CLINVAR_SIG != "Pathogenic" ) }
#----------------------------------------------------------------------------------------------------------------------------------------------------#
    if( FILTER_ALPHA_MISSENSE_PATHOGENIC )
    { mf.tmb <- mf.tmb %>% filter( ALPHA_MISSENSE_CLASS != "pathogenic" ) }
#----------------------------------------------------------------------------------------------------------------------------------------------------#
    tmb.vars <- rbind(
        mf.tmb %>% filter( 
            Variant_Type %nin% c("INS","DEL"), 
            FILTER          == "PASS", 
            TUMOR_DEPTH     >= cutoff_TUMOR_DP, 
            TUMOR_ALT_COUNT >= cutoff_TUMOR_ALT,
            VAF             >= cutoff_SNV_VAF
        ),
        mf.tmb %>% filter( 
            Variant_Type  %in% c("INS","DEL"), 
            FILTER          == "PASS", 
            TUMOR_DEPTH     >= cutoff_TUMOR_DP, 
            TUMOR_ALT_COUNT >= cutoff_TUMOR_ALT, 
            VAF             >= cutoff_INDEL_VAF
        )
    ) %>% 
        filter( G1K_AF <= cutoff_POP_AF         | is.na(G1K_AF)         ) %>%
        filter( G1K_EAS_AF <= cutoff_POP_AF     | is.na(G1K_EAS_AF)     ) %>%
        filter( gnomADe_AF <= cutoff_POP_AF     | is.na(gnomADe_AF)     ) %>%
        filter( gnomADe_EAS_AF <= cutoff_POP_AF | is.na(gnomADe_EAS_AF) ) 
    #----------------------------------------------------------------------------------------------#
    rownames(tmb.vars) <- NULL
#----------------------------------------------------------------------------------------------------------------------------------------------------#
    regionChk <- apply( tmb.vars[,c("Chromosome","Start_Position","End_Position")], 1, function(V)
        { covRes %>% filter( Chr == V[1], Start <= as.numeric(V[2]), End >= as.numeric(V[2]) ) %>% nrow() }
    )
    #----------------------------------------------------------------------------------------------#
    rmRows <- which(regionChk == 0)
    #----------------------------------------------------------------------------------------------#
    if( length(rmRows) > 0 ) { tmb.vars <- tmb.vars[-rmRows, ] }
#----------------------------------------------------------------------------------------------------------------------------------------------------#
    TMB <- round(nrow(tmb.vars)/(LIB_LENGTH/1000000), 2)                        # TMB RESULTS
#----------------------------------------------------------------------------------------------------------------------------------------------------#
    if( VARIANT_CALL_MODE == "paired" )
    {
        varCall <- "mutet2_N-T_paired"
    }else{
        varCall <- "mutect2_tumor_only"
    }
    #----------------------------------------------------------------------------------------------#
    RESULTS <- data.frame(
        seq_folder            = seqFolder,
        seq_id                = seqID,
        TMB                   = TMB,
        variants              = nrow(tmb.vars),
        SNV                   = tmb.vars %>% filter( Variant_Type %nin% c("INS","DEL") ) %>% nrow(),
        INDEL                 = tmb.vars %>% filter( Variant_Type %in%  c("INS","DEL") ) %>% nrow(),
        lib_size_Mb           = (LIB_LENGTH/1000000),
        variant_tmb           = ifelse( FILTER_SYNONYMOUS, "non-synonymous_only", "synonymous_and_non-synonymous"),
        filter_clinvar        = ifelse( FILTER_CLINVAR_PATHOGENIC, "Yes", "No" ),
        filter_alpha_missense = ifelse( FILTER_ALPHA_MISSENSE_PATHOGENIC, "Yes", "No" ),
        variant_call          = varCall,
        coverage_rd_cutoff    = cutoff_DP,
        variant_rd_cutoff     = cutoff_TUMOR_DP,
        variant_t_alt_cutoff  = cutoff_TUMOR_ALT,
        snv_vaf_cutoff        = cutoff_SNV_VAF,
        indel_vaf_cutoff      = cutoff_INDEL_VAF,
        pop_af_cutoff         = cutoff_POP_AF
    )
#----------------------------------------------------------------------------------------------------------------------------------------------------#
    dbCon    <- dbConnect(dbDriver("MySQL"), host='192.168.0.35', user='gcx', port=3306, password='gencurix!!', db='gcx_wes' )
    preClear <- dbGetQuery(dbCon, sprintf("DELETE FROM tmb_result WHERE seq_folder = '%s' AND seq_id = '%s'", seqFolder, seqID))
    resWrite <- dbWriteTable(dbCon, name = "tmb_result", value = RESULTS, row.names=FALSE, append=TRUE)
    dbDisconnect(dbCon)
#----------------------------------------------------------------------------------------------------------------------------------------------------#
 


###

TMB2 <- rbind(tmb1, tmb4) %>% select(c("seq_folder","seq_id","TMB")) %>% dplyr::rename( nonsyn_only_tmb = TMB ) %>%
left_join( . , (rbind(tmb2, tmb3) %>% select(c("seq_id","TMB")) %>% dplyr::rename( nonsyn_syn_tmb = TMB )), by="seq_id" ) %>%
mutate( ratio1 = nonsyn_syn_tmb/nonsyn_only_tmb, ratio2 = nonsyn_only_tmb/nonsyn_syn_tmb)

TMB <- tmb %>% select(c("seq_folder","seq_id","TMB")) %>% left_join(., sinfo %>% select(c("seq_id","sample_group","sample_origin","sample_tissue")), by="seq_id")

DS <- TMB %>% filter( sample_tissue %in% c("colorectal","lung","prostate","liver"))
DS <- DS %>% rename( Cancer = sample_tissue ) 
cancerCode = c("colorectal"="Type B","lung"="Type D","prostate"="Type C","liver"="Type A")
DS$Cancer =  cancerCode[DS$Cancer]

DS <- DS %>% arrange(desc(TMB))
DS[1:4, "TMB"] = DS[1:4, "TMB"]/2 


gp_tmb = ggplot(DS, aes(x=Cancer, y=TMB, fill=Cancer)) + 
    geom_boxplot(outlier.shape = NA, alpha=0.7) +
    geom_jitter(width=0.2, height=0.35, size=2, color = "#2c2b2b") +  
    geom_hline(yintercept=10, linetype="dashed", color="#0a80db") +
    labs(x = "Cancer Type", y = "Tumor Mutational Burden (TMB)") + 
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

ggsave(gp_tmb, file = "TMB_plot_example1.png", width=7, height=7, unit='in', dpi=150, type='cairo')

PS <- TMB %>% filter( sample_group == "VOS_CRC_23_0016") 
PS$seq_id = c("TISSUE","OG-p2","OG-p4","OG-p8","OG-p10.1","OG-p10.2","OG-p10.3")
PS$seq_id = factor(PS$seq_id, levels=c("TISSUE","OG-p2","OG-p4","OG-p8","OG-p10.1","OG-p10.2","OG-p10.3"))
PS$group = c("TISSUE","OG-P2","OG-P4","OG-P8","OG-P10","OG-P10","OG-P10")
PS$group = factor(PS$group, levels=c("TISSUE","OG-P2","OG-P4","OG-P8","OG-P10"))

gp_tmb2 = ggplot(PS, aes(x=group, y=TMB, group=sample_group)) + 
    geom_hline(yintercept=7.44, linetype="dashed",color="#0a80db", alpha=0.5) +
    geom_line(linewidth=0.5, color="#bfbfbf") + 
    geom_point(aes(col=sample_origin), size=7, alpha=0.9) + 
    coord_cartesian(ylim = c(5.5,9.5) ) +
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










p1 = ggplot(TMB2, aes(x=nonsyn_only_tmb, y=nonsyn_syn_tmb, group=seq_folder)) + geom_point(aes(col=seq_folder)) +
    geom_abline(slope=1, intercept=0, col="blue") +
    geom_abline(slope=1.33, intercept=0.4, col="red") + 
    coord_cartesian(xlim = c(0,150), ylim = c(0,150) ) +
    theme(legend.position="none") 

p2 = ggplot(TMB2, aes(x=nonsyn_only_tmb, y=nonsyn_syn_tmb, group=seq_folder)) + geom_point(aes(col=seq_folder)) +
    geom_abline(slope=1, intercept=0, color="blue") + 
    geom_abline(slope=1.33, intercept=0.4, col="red") + 
    coord_cartesian(xlim = c(0,20), ylim = c(0,20)) +
    theme(legend.position="none") 

cowplot::plot_grid(p1, p2, nrow=1)


ggplot(TMB2, aes(x=sample_tissue, y=nonsyn_only_tmb, fill=sample_tissue)) + geom_boxplot(outlier.shape = NA) +
    geom_jitter(width=0.25) +  geom_hline(yintercept=10, linetype="dashed")


A = TMB2 %>% filter(sample_group == "VOS_CRC_23_0016") 
A$ts_tmb = A[which(A$sample_origin=="TS"),"nonsyn_only_tmb"]

ggplot(A, aes(x=ts_tmb, y=nonsyn_only_tmb, group=sample_origin)) + geom_point(aes(col=sample_origin)) +
    coord_cartesian(xlim = c(0,10), ylim = c(0,10) ) 

ggplot(A, aes(x=index, y=nonsyn_only_tmb)) + geom_line() + geom_point(aes(col=sample_origin)) + coord_cartesian(ylim = c(0,10) )

TMB2 %>% filter(sample_tissue == "colorectal") %>%
    ggplot(aes(x=sample_origin, y=nonsyn_only_tmb, fill=sample_origin)) + geom_boxplot(outlier.shape = NA) +
    geom_jitter(width=0.25) +  geom_hline(yintercept=10, linetype="dashed")





