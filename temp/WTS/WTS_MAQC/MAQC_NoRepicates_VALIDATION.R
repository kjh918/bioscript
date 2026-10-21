
library(seqc)
library(plyr)
library(dplyr)
library(Hmisc)
library(parallel)
library(data.table)
library(RMySQL)

#---------------------------------------------------------------------------------------------------
    source("/storage/home/kangsm/runScripts/WTS_Fun.AnalysisModules.R")
    dbCon     <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db="public_database" )
    GENE_INFO <- dbGetQuery(dbCon, "SELECT symbol,name,locus_group,locus_type,entrez,ens_geneid, hgnc_id FROM gene_info_hgnc") %>% unique()
    dbDisconnect(dbCon)

    gcx_maqc_seqid =c("MAQC_II_A_rep1","MAQC_II_A_rep2","MAQC_II_A_rep3","MAQC_II_B_rep1","MAQC_II_B_rep2","MAQC_II_B_rep3")

      <- Create.Expression.Values.Matrix(
        SeqIdList=gcx_maqc_seqid,
        ValueType="count",
        DataSource="rds"
    )
    rownames(GCX_MAQC_COUNT) = sapply( rownames(GCX_MAQC_COUNT), function(y) unlist(strsplit(y, "\\."))[1] )
    geneinfo   <- GENE_INFO[match(rownames(GCX_MAQC_COUNT), GENE_INFO$ens_geneid), c("symbol","entrez","ens_geneid","name","locus_type")] %>% dplyr::rename( gene=symbol, gene_type=locus_type )


    norep_seqid =c("MAQC_II_NoRep_A_S75_rep1","MAQC_II_NoRep_A_S75_rep2","MAQC_II_NoRep_B_S75_rep1","MAQC_II_NoRep_B_S75_rep2")
    
    QUANT <- read.delim(sprintf("/storage/home/kangsm/analysis/gcx.rnaseq.maqc/no_rep_subsampling_test/%s/quant/rsem/%s.genes.results", norep_seqid[1], norep_seqid[1]))
    QUANT <- QUANT[,c(1,5)]
    colnames(QUANT)[ncol(QUANT)] <- norep_seqid[1]
    for( k in 2:length(norep_seqid) )
    {
        quant <- read.delim(sprintf("/storage/home/kangsm/analysis/gcx.rnaseq.maqc/no_rep_subsampling_test/%s/quant/rsem/%s.genes.results", norep_seqid[k], norep_seqid[k]))
        quant <- quant[,c(1,5)]
        QUANT <- QUANT %>% left_join(., quant, by="gene_id")
        colnames(QUANT)[ncol(QUANT)] <- norep_seqid[k]
    }
    rownames(QUANT) = sapply(QUANT$gene_id, function(w) unlist(strsplit(w, "\\."))[1] )
    QUANT = QUANT[,-1]
    for( i in 1:4 ){ QUANT[,i] = round(QUANT[,i], 0) }

    deg_analysis_config = "/storage/home/kangsm/analysis/gcx.rnaseq.maqc/no_rep_subsampling_test/meta/WTS_MAQC_deg_analysis_set_pairs.yaml"
    ControlTreatPairList <- read_yaml(deg_analysis_config)
    ControlTreatPairList <- lapply(ControlTreatPairList, function(ctp) lapply(ctp, function(z) if(!is.null(z)){ unlist(strsplit(z, ","))} ))

    SampleList = c("MAQC_II_A","MAQC_II_B","MAQC_II_NoRep_A_S75","MAQC_II_NoRep_B_S75")
    names(SampleList) = SampleList

    ainfo <- data.frame()
    for( k in 1:length(ControlTreatPairList) )
    {
        ainfo <- rbind(ainfo, 
            rbind(
                data.frame(
                    seq_folder      = SeqFolderId,
                    analysis_id     = k,
                    analysis_group  = "CTRL",
                    sample_id       = SampleList[ names(SampleList) == ControlTreatPairList[[k]]$ctrl ],
                    sample_group    = ControlTreatPairList[[k]]$ctrl,
                    genome_assembly = GenomeAssembly
                ),
                data.frame(
                    seq_folder      = SeqFolderId,
                    analysis_id     = k,
                    analysis_group  = "TREAT",
                    sample_id       = SampleList[ names(SampleList) == ControlTreatPairList[[k]]$treat ],
                    sample_group    = ControlTreatPairList[[k]]$treat,
                    genome_assembly = GenomeAssembly
                ) 
            )
        )
    }
    ainfo$analysis_id <- paste0("AID", formatC(ainfo$analysis_id, digits=0, flag="0", width=3))

    ANALYSIS_SET_LIST <- list(
        rbind(
            data.frame(
                seq_folder      = "no_rep_subsampling_test",
                analysis_id     = "AID001",
                analysis_group  = "CTRL",
                sampel_id       = c("MAQC_II_A_rep1","MAQC_II_A_rep2","MAQC_II_A_rep3"),
                sampel_group    = "MAQC_II_A",
                genome_assembly = "hg19"
            ),
            data.frame(
                seq_folder      = "no_rep_subsampling_test",
                analysis_id     = "AID001",
                analysis_group  = "TREAT",
                sampel_id       = c("MAQC_II_B_rep1","MAQC_II_B_rep2","MAQC_II_B_rep3"),
                sampel_group    = "MAQC_II_B",
                genome_assembly = "hg19"
            )
        ),
        rbind(
            data.frame(
                seq_folder      = "no_rep_subsampling_test",
                analysis_id     = "AID002",
                analysis_group  = "CTRL",
                sampel_id       = c("MAQC_II_NoRep_A_S75_rep1","MAQC_II_NoRep_A_S75_rep2"),
                sampel_group    = "MAQC_II_NoRep_A_S75",
                genome_assembly = "hg19"
            ),
            data.frame(
                seq_folder      = "no_rep_subsampling_test",
                analysis_id     = "AID002",
                analysis_group  = "TREAT",
                sampel_id       = c("MAQC_II_NoRep_B_S75_rep1","MAQC_II_NoRep_B_S75_rep2"),
                sampel_group    = "MAQC_II_NoRep_B_S75",
                genome_assembly = "hg19"
            )
        )
    )

    ExprMatList = list( GCX_MAQC_COUNT, QUANT )


    deg_analysis_res_list = lapply( 1:length(ANALYSIS_SET_LIST), function(j) 
    {
        AnalysisDesign=ANALYSIS_SET_LIST[[j]] 
        ExprCountMatrix=ExprMatList[[j]]

        FilterDesign <- model.matrix( ~ AnalysisDesign$analysis_group )
        keepGenes    <- filterByExpr( ExprCountMatrix, design=FilterDesign, group=AnalysisDesign$analysis_group, min.count=10, min.total.count=15 )

        CountMat  <- ExprMatList[[j]][, AnalysisDesign$sampel_id]
        dds       <- suppressMessages(DESeqDataSetFromMatrix(countData=CountMat, colData=AnalysisDesign, design= ~ analysis_group))
        dds       <- suppressMessages(estimateSizeFactors(dds))
        dds       <- suppressMessages(estimateDispersions(dds))
        dds_res   <- suppressMessages(DESeq(dds, betaPrior=FALSE))
        deseq_res <- suppressMessages(results(dds_res, pAdjustMethod = "BH", independentFiltering=FALSE))
        DESEQ_RES <- data.frame( ens_geneid = rownames(deseq_res), deseq_res )
        DESEQ_RES$keep <- keepGenes[ DESEQ_RES$ens_geneid ]
        #DESEQ_RES[ which(rmGenes), "keep" ] <- FALSE
        DESEQ_RES <- Calculate.DESeq.Scores( DESeqAnalysisResult=DESEQ_RES, ResultSort=TRUE )
        GeneSpaceSize <- DESEQ_RES %>% filter( keep ) %>% nrow(.)

        DESEQ_RES$deg <- ""
        DESEQ_RES[which( 
            DESEQ_RES$log2FoldChange >= log2(1.3) &
            DESEQ_RES$pvalue         <= 0.05        &
            DESEQ_RES$padj           <= 0.05           
        ), "deg" ] <- "UP"
        DESEQ_RES[which( 
            DESEQ_RES$log2FoldChange <= -log2(1.3) &
            DESEQ_RES$pvalue         <= 0.05         &
            DESEQ_RES$padj           <= 0.05           
        ), "deg" ] <- "DOWN"

        if( length(which(DESEQ_RES$deg != "")) > GeneSpaceSize * 0.11 )
        {
            # adjust DEGs : too many DEGs --------------------------------------
            message(sprintf("|---> too many DEGs. DEG ratio = %s", round(length(which(DESEQ_RES$deg != ""))/GeneSpaceSize, 2) ))
            message("|---> reducing DEGs to 10% of gene-space-size...")
            deg_numbers      <- as.integer(GeneSpaceSize * 0.1)
            deg_score_cutoff <- rev(sort(abs(DESEQ_RES$score)))[deg_numbers]
            DESEQ_RES$deg    <- ""
            DESEQ_RES[which( abs(DESEQ_RES$score) >= deg_score_cutoff & DESEQ_RES$log2FoldChange > 0 ), "deg"] <- "UP"
            DESEQ_RES[which( abs(DESEQ_RES$score) >= deg_score_cutoff & DESEQ_RES$log2FoldChange < 0 ), "deg"] <- "DOWN"
        }
        return(DESEQ_RES)
    })

    A = lapply( deg_analysis_res_list, function(y) y %>% arrange(desc(score)) %>% mutate(rank =1:nrow(.)) )

    REP = A[[1]][,c(1,3,10,11)]
    NO_REP = A[[2]][,c(1,3,10,11)]

    RD = data.frame( REP=REP, NO_REP=NO_REP[match(REP$ens_geneid, NO_REP$ens_geneid), ])
    RD_FILTERED = RD[ apply(RD[,c("REP.deg","NO_REP.deg")], 1, function(z) ifelse( all(z == ""), FALSE, TRUE)  ), ]
    ggplot(RD_FILTERED) + geom_point(aes(x=REP.rank, y=NO_REP.rank))

    ggplot(RD_FILTERED) + geom_point(aes(x=REP.log2FoldChange, y=NO_REP.log2FoldChange))


    RD_FILTERED2 = RD_FILTERED[,c("REP.deg","NO_REP.deg")]
    RD_FILTERED2[RD_FILTERED2 == "UP"] = 1
    RD_FILTERED2[RD_FILTERED2 == "DOWN"] = 1
    RD_FILTERED2[RD_FILTERED2 == ""] = 0
    RD_FILTERED2[,1] = as.numeric(RD_FILTERED2[,1])
    RD_FILTERED2[,2] = as.numeric(RD_FILTERED2[,2])
    round((1-Mercator::binaryDistance(as.matrix(RD_FILTERED2), metric="sokalMichener")), 3)
    RD_FILTERED[names(which(apply(RD_FILTERED2, 1, function(z) ifelse(z[1] != z[2] , TRUE, FALSE)))), ]


    QUANT <- read.delim(sprintf("/data/wts/WTS_MAQC/%s/quant/rsem/%s.genes.results", gcx_maqc_seqid[1], gcx_maqc_seqid[1]))
    QUANT <- QUANT[,c(1,4)]
    colnames(QUANT)[ncol(QUANT)] <- paste0("effLen.", gcx_maqc_seqid[1])

    for( k in 2:length(gcx_maqc_seqid) )
    {
        quant <- read.delim(sprintf("/data/wts/WTS_MAQC/%s/quant/rsem/%s.genes.results", gcx_maqc_seqid[k], gcx_maqc_seqid[k]))
        quant <- quant[,c(1,4)]
        QUANT <- QUANT %>% left_join(., quant, by="gene_id")
        colnames(QUANT)[ncol(QUANT)] <- paste0("effLen.", gcx_maqc_seqid[k])
    }
    EFF_LEN = QUANT
    EFF_LEN$MAQC_II_A = apply(EFF_LEN[,c("effLen.MAQC_II_A_rep1","effLen.MAQC_II_A_rep2","effLen.MAQC_II_A_rep3")], 1, function(y) mean(y))
    EFF_LEN$MAQC_II_B = apply(EFF_LEN[,c("effLen.MAQC_II_B_rep1","effLen.MAQC_II_B_rep2","effLen.MAQC_II_B_rep3")], 1, function(y) mean(y))
    EFF_LEN$gene_id = sapply(EFF_LEN$gene_id, function(y) unlist(strsplit(y, "\\."))[1])
    EFF_LEN = EFF_LEN[match(GENE_SPACE, EFF_LEN$gene_id), ]   


    GCX_MAQC <- GCX_MAQC_COUNT[ GENE_SPACE, ]

    GCX_MAQC$MAQC_A_gcx_tpm_rep1 = countToTpm( counts=GCX_MAQC$MAQC_A_rep1, effLen=EFF_LEN$MAQC_A )
    GCX_MAQC$MAQC_A_gcx_tpm_rep2 = countToTpm( counts=GCX_MAQC$MAQC_A_rep2, effLen=EFF_LEN$MAQC_A )
    GCX_MAQC$MAQC_B_gcx_tpm_rep1 = countToTpm( counts=GCX_MAQC$MAQC_B_rep1, effLen=EFF_LEN$MAQC_B )
    GCX_MAQC$MAQC_B_gcx_tpm_rep2 = countToTpm( counts=GCX_MAQC$MAQC_B_rep2, effLen=EFF_LEN$MAQC_B )
    GCX_MAQC$MAQC_II_A_gcx_tpm_rep1 = tpm( counts=GCX_MAQC$MAQC_II_A_rep1, lengths=EFF_LEN$MAQC_II_A )
    GCX_MAQC$MAQC_II_A_gcx_tpm_rep2 = tpm( counts=GCX_MAQC$MAQC_II_A_rep2, lengths=EFF_LEN$MAQC_II_A )
    GCX_MAQC$MAQC_II_A_gcx_tpm_rep3 = tpm( counts=GCX_MAQC$MAQC_II_A_rep3, lengths=EFF_LEN$MAQC_II_A )
    GCX_MAQC$MAQC_II_B_gcx_tpm_rep1 = tpm( counts=GCX_MAQC$MAQC_II_B_rep1, lengths=EFF_LEN$MAQC_II_B )
    GCX_MAQC$MAQC_II_B_gcx_tpm_rep2 = tpm( counts=GCX_MAQC$MAQC_II_B_rep2, lengths=EFF_LEN$MAQC_II_B )
    GCX_MAQC$MAQC_II_B_gcx_tpm_rep3 = tpm( counts=GCX_MAQC$MAQC_II_B_rep3, lengths=EFF_LEN$MAQC_II_B )

    GCX_MAQC$MAQC_A = apply(GCX_MAQC[,c("MAQC_A_gcx_tpm_rep1","MAQC_A_gcx_tpm_rep2")], 1, function(y) mean(y))
    GCX_MAQC$MAQC_B = apply(GCX_MAQC[,c("MAQC_B_gcx_tpm_rep1","MAQC_B_gcx_tpm_rep2")], 1, function(y) mean(y))
    GCX_MAQC$MAQC_II_A = apply(GCX_MAQC[,c("MAQC_II_A_gcx_tpm_rep1","MAQC_II_A_gcx_tpm_rep2","MAQC_II_A_gcx_tpm_rep3")], 1, function(y) mean(y))
    GCX_MAQC$MAQC_II_B = apply(GCX_MAQC[,c("MAQC_II_B_gcx_tpm_rep1","MAQC_II_B_gcx_tpm_rep2","MAQC_II_B_gcx_tpm_rep3")], 1, function(y) mean(y))

    GCX_MAQC$MAQC_A = round(log2(GCX_MAQC$MAQC_A +1), 3)
    GCX_MAQC$MAQC_B = round(log2(GCX_MAQC$MAQC_B +1), 3)
    GCX_MAQC$MAQC_II_A = round(log2(GCX_MAQC$MAQC_II_A +1), 3)
    GCX_MAQC$MAQC_II_B = round(log2(GCX_MAQC$MAQC_II_B +1), 3)
    
    MAQC_GCX = data.frame(
        gene = rownames(GCX_MAQC),
        MAQC_A = GCX_MAQC$MAQC_A,
        MAQC_B = GCX_MAQC$MAQC_B,
        MAQC_II_A = GCX_MAQC$MAQC_II_A,
        MAQC_II_B = GCX_MAQC$MAQC_II_B
    )

    saveRDS(MAQC_GCX, file="MAQC_GCX_log2TPM_matrix.rds")

#---------------------------------------------------------------------------------------------------


    saveRDS(MAQC_RNASEQ, file="MAQC_RNAseq_log2TPM_matrix.rds")
    saveRDS(MAQC_RT_PCR, file="MAQC_RT-PCR_ct_values.rds")
    saveRDS(MAQC_ARRAY,  file="MAQC_Microarray_Affymetrix.rds")

    MAQC <- data.frame(
        ens_geneid = GENE_SPACE,
        MAQC_RT_PCR[match(GENE_SPACE, MAQC_RT_PCR$ens_geneid), c(3,4,5,6)]
    )
    colnames(MAQC) = c("ens_geneid","gene","entrez","RTPCR.MAQC_A","RTPCR.MAQC_B")

    MAQC <- data.frame(
        MAQC,
        MAQC_ARRAY[match(MAQC$gene, MAQC_ARRAY$Gene), c(2,3)]
    )
    colnames(MAQC)[6:7] = c("microarray.MAQC_A","microarray.MAQC_B")

    MAQC <- data.frame(
        MAQC,
        MAQC_RNASEQ[match(MAQC$ens_geneid, MAQC_RNASEQ$gene), -1]
    )

    MAQC <- data.frame(
        MAQC,
        MAQC_GCX[match(MAQC$ens_geneid, MAQC_GCX$gene), -1]
    )
    colnames(MAQC)[16:19] = c("GCX.SE.MAQC_A","GCX.SE.MAQC_B","GCX.PE.MAQC_A","GCX.PE.MAQC_B")

    saveRDS(MAQC, file="MAQC_data_matrix_full.rds")

#---------------------------------------------------------------------------------------------------

###
library(ggplot2)
library(Cairo)
grDevices::X11.options(type='cairo')
options(device='x11')
options(stringsAsFactors=FALSE) 
options(bitmapType='cairo')

MAQC$norm.RTPCR.MAQC_A = round(log2(2^(-(MAQC$RTPCR.MAQC_A - 35))+1), 3)
MAQC$norm.RTPCR.MAQC_B = round(log2(2^(-(MAQC$RTPCR.MAQC_B - 35))+1), 3)

MAQC_Rev1 = MAQC %>% filter( RTPCR.MAQC_A >= 11 & RTPCR.MAQC_A <= 32 ) 
MAQC_Rev2 = MAQC %>% filter( RTPCR.MAQC_B >= 11 & RTPCR.MAQC_B <= 32 ) 
MAQC_Rev3 = MAQC %>% filter( !is.na(RTPCR.MAQC_A) ) %>% filter( RTPCR.MAQC_A >= 11 & RTPCR.MAQC_A <= 32 ) 
MAQC_Rev4 = MAQC %>% filter( !is.na(RTPCR.MAQC_B) ) %>% filter( RTPCR.MAQC_B >= 11 & RTPCR.MAQC_B <= 32 ) 
MAQC_REV2 = MAQC_Rev3 %>% filter( !is.na(microarray.MAQC_A))
MAQC_Rev6 = MAQC_Rev4 %>% filter( !is.na(microarray.MAQC_B))

MAQC_REV <- MAQC %>% filter(
    RTPCR.MAQC_A >= 11 & RTPCR.MAQC_A <= 32, RTPCR.MAQC_B >= 11 & RTPCR.MAQC_B <= 32, 
    !is.na(RTPCR.MAQC_A), !is.na(RTPCR.MAQC_B)
)
MAQC_REV2 <- MAQC_REV %>% filter( !is.na(microarray.MAQC_A), !is.na(microarray.MAQC_B) )




MAQC_FC = MAQC_REV2 %>% mutate(
    fc.RTPCR      = round(RTPCR.MAQC_B/RTPCR.MAQC_A, 3),
    fc.microarray = round(microarray.MAQC_B/microarray.MAQC_A, 3),
    fc.cufflinks  = round(cufflinks.MAQC_B/cufflinks.MAQC_A, 3),
    fc.kallisto   = round(kallisto.MAQC_B/kallisto.MAQC_A, 3),
    fc.salmon     = round(salmon.MAQC_B/salmon.MAQC_A, 3),
    fc.htseq      = round(htseq.MAQC_B/htseq.MAQC_A, 3),
    fc.GCX.SE     = round(GCX.SE.MAQC_B/GCX.SE.MAQC_A, 3),
    fc.GCX.PE     = round(GCX.PE.MAQC_B/GCX.PE.MAQC_A, 3),
    fc.norm.RTPCR = round(norm.RTPCR.MAQC_B/norm.RTPCR.MAQC_A, 3)
) %>% dplyr::select(c("ens_geneid","gene","entrez","fc.RTPCR","fc.norm.RTPCR","fc.microarray","fc.cufflinks","fc.kallisto","fc.salmon","fc.htseq","fc.GCX.SE","fc.GCX.PE"))

MAQC_FC2 = MAQC_REV2 %>% mutate(
    fc.RTPCR      = round(log2(RTPCR.MAQC_B+1) - log2(RTPCR.MAQC_A+1), 3),
    fc.microarray = round(log2(microarray.MAQC_B+1) - log2(microarray.MAQC_A+1), 3),
    fc.cufflinks  = round(log2(cufflinks.MAQC_B+1) - log2(cufflinks.MAQC_A+1), 3),
    fc.kallisto   = round(log2(kallisto.MAQC_B+1) - log2(kallisto.MAQC_A+1), 3),
    fc.salmon     = round(log2(salmon.MAQC_B+1) - log2(salmon.MAQC_A+1), 3),
    fc.htseq      = round(log2(htseq.MAQC_B+1) - log2(htseq.MAQC_A+1), 3),
    fc.GCX.SE     = round(log2(GCX.SE.MAQC_B+1) - log2(GCX.SE.MAQC_A+1), 3),
    fc.GCX.PE     = round(log2(GCX.PE.MAQC_B+1) - log2(GCX.PE.MAQC_A+1), 3),
    fc.norm.RTPCR = round(log2(norm.RTPCR.MAQC_B+1) - log2(norm.RTPCR.MAQC_A+1), 3),
) %>% dplyr::select(c("ens_geneid","gene","entrez","fc.RTPCR","fc.norm.RTPCR","fc.microarray","fc.cufflinks","fc.kallisto","fc.salmon","fc.htseq","fc.GCX.SE","fc.GCX.PE"))


png( "MAQC_expression_foldchange_scatter.png", width=16, height=7.5, units="in",res=150 )
par(mfrow=c(3,7))
plot(x=MAQC_REV2$norm.RTPCR.MAQC_A, y=MAQC_REV2$microarray.MAQC_A, pch=19, col=alpha("#2d2d2d", 0.5), cex=0.8, xlab="MAQC-A : RT-PCR", ylab="MAQC-A : Microarray", main="Microarray", ylim=c(0,14))
plot(x=MAQC_REV2$norm.RTPCR.MAQC_A, y=MAQC_REV2$cufflinks.MAQC_A,  pch=19, col=alpha("#b85450", 0.5), cex=0.8, xlab="MAQC-A : RT-PCR", ylab="MAQC-A : RNAseq Cufflinks", main="Tophat-Cufflinks", ylim=c(0,14))
plot(x=MAQC_REV2$norm.RTPCR.MAQC_A, y=MAQC_REV2$kallisto.MAQC_A,   pch=19, col=alpha("#8000FF", 0.5), cex=0.8, xlab="MAQC-A : RT-PCR", ylab="MAQC-A : RNAseq Kallisto", main="Kallisto", ylim=c(0,14))
plot(x=MAQC_REV2$norm.RTPCR.MAQC_A, y=MAQC_REV2$salmon.MAQC_A,     pch=19, col=alpha("#FF7F0E", 0.5), cex=0.8, xlab="MAQC-A : RT-PCR", ylab="MAQC-A : RNAseq Salmon", main="Salmon", ylim=c(0,14))
plot(x=MAQC_REV2$norm.RTPCR.MAQC_A, y=MAQC_REV2$htseq.MAQC_A,      pch=19, col=alpha("#006EAF", 0.5), cex=0.8, xlab="MAQC-A : RT-PCR", ylab="MAQC-A : RNAseq HTseq", main="TopHat-HTseq", ylim=c(0,14))
plot(x=MAQC_REV2$norm.RTPCR.MAQC_A, y=MAQC_REV2$GCX.SE.MAQC_A,     pch=19, col=alpha("#2D7600", 0.5), cex=0.8, xlab="MAQC-A : RT-PCR", ylab="MAQC-A : RNAseq GCX-SE", main="STAR-RSEM (SE)", ylim=c(0,14))
plot(x=MAQC_REV2$norm.RTPCR.MAQC_A, y=MAQC_REV2$GCX.PE.MAQC_A,     pch=19, col=alpha("#00b22c", 0.5), cex=0.8, xlab="MAQC-A : RT-PCR", ylab="MAQC-A : RNAseq GCX-PE", main="STAR-RSEM (PE)", ylim=c(0,14))

plot(x=MAQC_REV2$norm.RTPCR.MAQC_B, y=MAQC_REV2$microarray.MAQC_B, pch=19, col=alpha("#2d2d2d", 0.5), cex=0.8, xlab="MAQC-B : RT-PCR", ylab="MAQC-B : Microarray", main="Microarray", ylim=c(0,14))
plot(x=MAQC_REV2$norm.RTPCR.MAQC_B, y=MAQC_REV2$cufflinks.MAQC_B,  pch=19, col=alpha("#b85450", 0.5), cex=0.8, xlab="MAQC-B : RT-PCR", ylab="MAQC-B : RNAseq Cufflinks", main="Tophat-Cufflinks", ylim=c(0,14))
plot(x=MAQC_REV2$norm.RTPCR.MAQC_B, y=MAQC_REV2$kallisto.MAQC_B,   pch=19, col=alpha("#8000FF", 0.5), cex=0.8, xlab="MAQC-B : RT-PCR", ylab="MAQC-B : RNAseq Kallisto", main="Kallisto", ylim=c(0,14))
plot(x=MAQC_REV2$norm.RTPCR.MAQC_B, y=MAQC_REV2$salmon.MAQC_B,     pch=19, col=alpha("#FF7F0E", 0.5), cex=0.8, xlab="MAQC-B : RT-PCR", ylab="MAQC-B : RNAseq Salmon", main="Salmon", ylim=c(0,14))
plot(x=MAQC_REV2$norm.RTPCR.MAQC_B, y=MAQC_REV2$htseq.MAQC_B,      pch=19, col=alpha("#006EAF", 0.5), cex=0.8, xlab="MAQC-B : RT-PCR", ylab="MAQC-B : RNAseq HTseq", main="TopHat-HTseq", ylim=c(0,14))
plot(x=MAQC_REV2$norm.RTPCR.MAQC_B, y=MAQC_REV2$GCX.SE.MAQC_B,     pch=19, col=alpha("#2D7600", 0.5), cex=0.8, xlab="MAQC-B : RT-PCR", ylab="MAQC-B : RNAseq GCX-SE", main="STAR-RSEM (SE)", ylim=c(0,14))
plot(x=MAQC_REV2$norm.RTPCR.MAQC_B, y=MAQC_REV2$GCX.PE.MAQC_B,     pch=19, col=alpha("#00b22c", 0.5), cex=0.8, xlab="MAQC-B : RT-PCR", ylab="MAQC-B : RNAseq GCX-PE", main="STAR-RSEM (PE)", ylim=c(0,14))

plot(x=log2(MAQC_FC$fc.norm.RTPCR), y=log2(MAQC_FC$fc.microarray), pch=19, col=alpha("#2d2d2d", 0.5), cex=0.8, xlab="FoldChange RT-PCR", ylab="FoldChange Microarray", main="Microarray", xlim=c(-4,4), ylim=c(-10,10))
plot(x=log2(MAQC_FC$fc.norm.RTPCR), y=log2(MAQC_FC$fc.cufflinks),  pch=19, col=alpha("#b85450", 0.5), cex=0.8, xlab="FoldChange RT-PCR", ylab="FoldChange RNAseq Cufflinks", main="Tophat-Cufflinks", xlim=c(-4,4), ylim=c(-10,10))
plot(x=log2(MAQC_FC$fc.norm.RTPCR), y=log2(MAQC_FC$fc.kallisto),   pch=19, col=alpha("#8000FF", 0.5), cex=0.8, xlab="FoldChange RT-PCR", ylab="FoldChange RNAseq Kallisto", main="Kallisto", xlim=c(-4,4), ylim=c(-10,10))
plot(x=log2(MAQC_FC$fc.norm.RTPCR), y=log2(MAQC_FC$fc.salmon),     pch=19, col=alpha("#FF7F0E", 0.5), cex=0.8, xlab="FoldChange RT-PCR", ylab="FoldChange RNAseq Salmon", main="Salmon", xlim=c(-4,4), ylim=c(-10,10))
plot(x=log2(MAQC_FC$fc.norm.RTPCR), y=log2(MAQC_FC$fc.htseq),      pch=19, col=alpha("#006EAF", 0.5), cex=0.8, xlab="FoldChange RT-PCR", ylab="FoldChange RNAseq HTseq", main="TopHat-HTseq", xlim=c(-4,4), ylim=c(-10,10))
plot(x=log2(MAQC_FC$fc.norm.RTPCR), y=log2(MAQC_FC$fc.GCX.SE),     pch=19, col=alpha("#2D7600", 0.5), cex=0.8, xlab="FoldChange RT-PCR", ylab="FoldChange RNAseq GCX-SE", main="STAR-RSEM (SE)", xlim=c(-4,4), ylim=c(-10,10))
plot(x=log2(MAQC_FC$fc.norm.RTPCR), y=log2(MAQC_FC$fc.GCX.PE),     pch=19, col=alpha("#00b22c", 0.5), cex=0.8, xlab="FoldChange RT-PCR", ylab="FoldChange RNAseq GCX-PE", main="STAR-RSEM (PE)", xlim=c(-4,4), ylim=c(-10,10))
dev.off()



ExprValueCorrMatrix <- data.frame(rbind(
    c(round(cor(x=MAQC_REV2$norm.RTPCR.MAQC_A, y=MAQC_REV2$microarray.MAQC_A, method="spearman"),3), round(cor(x=MAQC_REV2$norm.RTPCR.MAQC_B, y=MAQC_REV2$microarray.MAQC_B, method="spearman"),3)),
    c(round(cor(x=MAQC_REV2$norm.RTPCR.MAQC_A, y=MAQC_REV2$cufflinks.MAQC_A,  method="spearman"),3), round(cor(x=MAQC_REV2$norm.RTPCR.MAQC_B, y=MAQC_REV2$cufflinks.MAQC_B,  method="spearman"),3)),
    c(round(cor(x=MAQC_REV2$norm.RTPCR.MAQC_A, y=MAQC_REV2$kallisto.MAQC_A,   method="spearman"),3), round(cor(x=MAQC_REV2$norm.RTPCR.MAQC_B, y=MAQC_REV2$kallisto.MAQC_B,   method="spearman"),3)),
    c(round(cor(x=MAQC_REV2$norm.RTPCR.MAQC_A, y=MAQC_REV2$salmon.MAQC_A,     method="spearman"),3), round(cor(x=MAQC_REV2$norm.RTPCR.MAQC_B, y=MAQC_REV2$salmon.MAQC_B,     method="spearman"),3)),
    c(round(cor(x=MAQC_REV2$norm.RTPCR.MAQC_A, y=MAQC_REV2$htseq.MAQC_A,      method="spearman"),3), round(cor(x=MAQC_REV2$norm.RTPCR.MAQC_B, y=MAQC_REV2$htseq.MAQC_B,      method="spearman"),3)),
    c(round(cor(x=MAQC_REV2$norm.RTPCR.MAQC_A, y=MAQC_REV2$GCX.SE.MAQC_A,     method="spearman"),3), round(cor(x=MAQC_REV2$norm.RTPCR.MAQC_B, y=MAQC_REV2$GCX.SE.MAQC_B,     method="spearman"),3)),
    c(round(cor(x=MAQC_REV2$norm.RTPCR.MAQC_A, y=MAQC_REV2$GCX.PE.MAQC_A,     method="spearman"),3), round(cor(x=MAQC_REV2$norm.RTPCR.MAQC_B, y=MAQC_REV2$GCX.PE.MAQC_B,     method="spearman"),3))
))
colnames(ExprValueCorrMatrix) = c("MAQC_A","MAQC_B")
rownames(ExprValueCorrMatrix) = c("Microarray","Tophat_Cufflinks","Kallisto","Salmon","Tophat_HTseq","GCX_SE_STAR_RSEM","GCX_PE_STAR_RSEM")
 

FoldChangeCorrMatrix <- data.frame(rbind(
    c(round(cor(x=MAQC_FC2$fc.norm.RTPCR, y=MAQC_FC2$fc.microarray, method="spearman"),3)),
    c(round(cor(x=MAQC_FC2$fc.norm.RTPCR, y=MAQC_FC2$fc.cufflinks,  method="spearman"),3)),
    c(round(cor(x=MAQC_FC2$fc.norm.RTPCR, y=MAQC_FC2$fc.kallisto,   method="spearman"),3)),
    c(round(cor(x=MAQC_FC2$fc.norm.RTPCR, y=MAQC_FC2$fc.salmon,     method="spearman"),3)),
    c(round(cor(x=MAQC_FC2$fc.norm.RTPCR, y=MAQC_FC2$fc.htseq,      method="spearman"),3)),
    c(round(cor(x=MAQC_FC2$fc.norm.RTPCR, y=MAQC_FC2$fc.GCX.SE,     method="spearman"),3)),
    c(round(cor(x=MAQC_FC2$fc.norm.RTPCR, y=MAQC_FC2$fc.GCX.PE,     method="spearman"),3))
))
colnames(FoldChangeCorrMatrix) = c("FoldChange.MAQC_B/MAQC_A)")
rownames(FoldChangeCorrMatrix) = c("Microarray","Tophat_Cufflinks","Kallisto","Salmon","Tophat_HTseq","GCX_SE_STAR_RSEM","GCX_PE_STAR_RSEM")

CorMat <- cbind(ExprValueCorrMatrix, FoldChangeCorrMatrix)
openxlsx::write.xlsx(CorMat, "maqc_corr_result.xlsx")

CM = rbind(
    data.frame( analysis=rownames(CorMat), group="Expression Correlation MAQC-A", value=CorMat$MAQC_A ),
    data.frame( analysis=rownames(CorMat), group="Expression Correlation MAQC-B", value=CorMat$MAQC_B ),
    data.frame( analysis=rownames(CorMat), group="FoldChange Correlation MAQC-B / MAQC-A", value=CorMat[,3] ) 
)
CM$analysis = factor(CM$analysis, levels=c("Microarray","Tophat_Cufflinks","Kallisto","Salmon","Tophat_HTseq","GCX_SE_STAR_RSEM","GCX_PE_STAR_RSEM"))
CorrColors = c("Expression Correlation MAQC-A"="#6C8EBF", "Expression Correlation MAQC-B"="#82B366", "FoldChange Correlation MAQC-B / MAQC-A"="#B85450" )

    corr_barplot1 = ggplot(CM, aes(x=analysis, y=value, group=group, fill=group)) + geom_bar( stat='identity', position='dodge' ) +
       coord_cartesian(ylim=c(0,1)) +
       scale_fill_manual(values=CorrColors) +
       scale_y_continuous(expand=c(0,0)) +
       labs(x="", y="Correlation", fill="") +
       theme(
            panel.background = element_rect(fill='white'), 
            panel.grid.major = element_line(color = '#dcdbdb', linetype = 'dotted'),
            axis.text.x  = element_text(colour = "#000000", face="bold", size=15, family="Nunito"), 
            axis.text.y  = element_text(colour = "#000000", size=11),
            axis.title.y = element_text(colour = "#000000", face="bold", size=16, family="Nunito"), 
            axis.line.x  = element_line(colour = "#1e1e1e"), 
            axis.line.y  = element_line(colour = "#1e1e1e"),
            plot.margin  = unit(c(0.6, 0.6, 0.6, 0.6), "cm"),
            legend.text  = element_text(colour = "#000000", size=12, family="Nunito"),
            legend.position = "bottom"
       )

    corr_barplot2 = ggplot(CM, aes(x=analysis, y=value, group=group, fill=group)) + geom_bar( stat='identity', position='dodge' ) +
       coord_cartesian(ylim=c(0.6,1)) +
       scale_fill_manual(values=CorrColors) +
       scale_y_continuous(expand=c(0,0)) +
       labs(x="", y="Correlation", fill="") +
       theme(
            panel.background = element_rect(fill='white'), 
            panel.grid.major = element_line(color = '#dcdbdb', linetype = 'dotted'),
            axis.text.x  = element_text(colour = "#000000", face="bold", size=13, family="Nunito"), 
            axis.text.y  = element_text(colour = "#000000", size=11),
            axis.title.y = element_text(colour = "#000000", face="bold", size=16, family="Nunito"), 
            axis.line.x  = element_line(colour = "#1e1e1e"), 
            axis.line.y  = element_line(colour = "#1e1e1e"),
            plot.margin  = unit(c(0.6, 0.6, 0.6, 0.6), "cm"),
            legend.text  = element_text(colour = "#000000", size=12, family="Nunito"),
       )


    corr_barplot3 = ggplot(CM, aes(x=analysis, y=value, group=group, fill=group)) + geom_bar( stat='identity', position='dodge' ) +
       coord_cartesian(ylim=c(0.6,1)) +
       scale_fill_manual(values=CorrColors) +
       scale_y_continuous(expand=c(0,0)) +
       labs(x="", y="Correlation", fill="") +
       theme(
            panel.background = element_rect(fill='white'), 
            panel.grid.major = element_line(color = '#dcdbdb', linetype = 'dotted'),
            axis.text.x  = element_text(colour = "#000000", face="bold", size=13, family="Nunito"), 
            axis.text.y  = element_text(colour = "#000000", size=11),
            axis.title.y = element_text(colour = "#000000", face="bold", size=16, family="Nunito"), 
            axis.line.x  = element_line(colour = "#1e1e1e"), 
            axis.line.y  = element_line(colour = "#1e1e1e"),
            plot.margin  = unit(c(0.6, 0.6, 0.6, 0.6), "cm"),
            legend.text  = element_text(colour = "#000000", size=12, family="Nunito"),
            legend.position = "none"
       )


    ggsave(corr_barplot1, file = "maqc_corr_barplot_type1.png", width=25, height=8, unit='in', dpi=150)
    ggsave(corr_barplot2, file = "maqc_corr_barplot_type2.png", width=25, height=8, unit='in', dpi=150)
    ggsave(corr_barplot3, file = "maqc_corr_barplot_type3.png", width=25, height=8, unit='in', dpi=150)











