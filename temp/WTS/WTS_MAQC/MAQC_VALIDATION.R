
library(seqc)
library(plyr)
library(dplyr)
library(Hmisc)
library(parallel)
library(data.table)
library(RMySQL)

#---/ MAQC microarry (affymetrix) and RT-PCR : GSE564457 -------------------------------------------
maqc_array_data_dir <- "/data/wts/WTS_MAQC/performance/GSE56457"
maqc_array_data_file <- sprintf("%s/GSE56457_SEQC_USF_GEOSub_HuGene2.txt", maqc_array_data_dir)
maqc_rtpcr_data_file <- sprintf("%s/GSE56457_SEQC_qPCR_GEOSub.txt", maqc_array_data_dir)

# microarray
maqc_array <- read.delim(maqc_array_data_file)
MAQC_ARRAY <- maqc_rt_pcr[, c(1,2,3,4,5,6,7,8,9)]
colnames(MAQC_ARRAY) <- c("GENE","MAQC_A_microarray_rep1","MAQC_A_microarray_rep2","MAQC_A_microarray_rep3","MAQC_A_microarray_rep4",
    "MAQC_B_microarray_rep1","MAQC_B_microarray_rep2","MAQC_B_microarray_rep3","MAQC_B_microarray_rep4")
MAQC_ARRAY <- data.frame(
    Gene = MAQC_ARRAY$GENE,
    MAQC_A = apply( MAQC_ARRAY[,c(2:5)], 1, function(y) round(mean(y), 3) ),
    MAQC_B = apply( MAQC_ARRAY[,c(6:9)], 1, function(y) round(mean(y), 3) )
)

# RT-PCR
maqc_rtpcr <- read.delim(maqc_rtpcr_data_file)
MAQC_RT_PCR <- maqc_rtpcr[,c("ID_REF","ensembl_gene","official_symbol","entrez_gene_id","SEQC_RTPCR_A","SEQC_RTPCR_B")]
colnames(MAQC_RT_PCR) <- c("ID_REF","ens_geneid","Gene","entrez","MAQC_A_R_PCR", "MAQC_B_RTPCR")

#---------------------------------------------------------------------------------------------------

#---/ MAQC RNAseq : GSE83402 -----------------------------------------------------------------------
maqc_rnaseq_dir <- "/data/wts/WTS_MAQC/performance/GSE83402"
method_list     <- c("clout"="cufflinks", "htout"="htseq", "klout"="kallisto", "salout"="salmon")

# HTseq
    maqc_a_r1 <- fread(sprintf("%s/%s", maqc_rnaseq_dir, list.files(maqc_rnaseq_dir, pattern=sprintf("A-1_%s", "htout") )))
    maqc_a_r2 <- fread(sprintf("%s/%s", maqc_rnaseq_dir, list.files(maqc_rnaseq_dir, pattern=sprintf("A-2_%s", "htout") )))
    maqc_b_r1 <- fread(sprintf("%s/%s", maqc_rnaseq_dir, list.files(maqc_rnaseq_dir, pattern=sprintf("B-1_%s", "htout") )))
    maqc_b_r2 <- fread(sprintf("%s/%s", maqc_rnaseq_dir, list.files(maqc_rnaseq_dir, pattern=sprintf("B-2_%s", "htout") )))
    colnames(maqc_a_r1) = c("ID","MAQC_A_htseq_count_rep1")
    colnames(maqc_a_r2) = c("ID","MAQC_A_htseq_count_rep2")
    colnames(maqc_b_r1) = c("ID","MAQC_B_htseq_count_rep1")
    colnames(maqc_b_r2) = c("ID","MAQC_B_htseq_count_rep2")
    maqc <- maqc_a_r1 %>% left_join(., maqc_a_r2, by="ID" ) %>% left_join(., maqc_b_r1, by="ID" ) %>% left_join(., maqc_b_r2, by="ID" ) 
    MAQC_HTseq <- maqc[grep("^ENS", maqc$ID), ] %>% as.data.frame()
    
# Salmon
    maqc_a_r1 <- fread(sprintf("%s/%s", maqc_rnaseq_dir, list.files(maqc_rnaseq_dir, pattern=sprintf("A-1_%s", "salout") )))
    maqc_a_r2 <- fread(sprintf("%s/%s", maqc_rnaseq_dir, list.files(maqc_rnaseq_dir, pattern=sprintf("A-2_%s", "salout") )))
    maqc_b_r1 <- fread(sprintf("%s/%s", maqc_rnaseq_dir, list.files(maqc_rnaseq_dir, pattern=sprintf("B-1_%s", "salout") )))
    maqc_b_r2 <- fread(sprintf("%s/%s", maqc_rnaseq_dir, list.files(maqc_rnaseq_dir, pattern=sprintf("B-2_%s", "salout") )))
    colnames(maqc_a_r1) = c("tx_id","length","eff_length","MAQC_A_salon_tpm_rep1","MAQC_A_salon_count_rep1")
    colnames(maqc_a_r2) = c("tx_id","length","eff_length","MAQC_A_salon_tpm_rep2","MAQC_A_salon_count_rep2")
    colnames(maqc_b_r1) = c("tx_id","length","eff_length","MAQC_B_salon_tpm_rep1","MAQC_B_salon_count_rep1")
    colnames(maqc_b_r2) = c("tx_id","length","eff_length","MAQC_B_salon_tpm_rep2","MAQC_B_salon_count_rep2")
    maqc <- maqc_a_r1 %>% left_join(., maqc_a_r2[,c("tx_id","MAQC_A_salon_tpm_rep2","MAQC_A_salon_count_rep2")], by="tx_id") %>%
        left_join(., maqc_b_r1[,c("tx_id","MAQC_B_salon_tpm_rep1","MAQC_B_salon_count_rep1")], by="tx_id") %>% 
        left_join(., maqc_b_r2[,c("tx_id","MAQC_B_salon_tpm_rep2","MAQC_B_salon_count_rep2")], by="tx_id")
    MAQC_salmon <- maqc

# Kallisto
    maqc_a_r1 <- fread(sprintf("%s/%s", maqc_rnaseq_dir, list.files(maqc_rnaseq_dir, pattern=sprintf("A-1_%s", "klout") )))
    maqc_a_r2 <- fread(sprintf("%s/%s", maqc_rnaseq_dir, list.files(maqc_rnaseq_dir, pattern=sprintf("A-2_%s", "klout") )))
    maqc_b_r1 <- fread(sprintf("%s/%s", maqc_rnaseq_dir, list.files(maqc_rnaseq_dir, pattern=sprintf("B-1_%s", "klout") )))
    maqc_b_r2 <- fread(sprintf("%s/%s", maqc_rnaseq_dir, list.files(maqc_rnaseq_dir, pattern=sprintf("B-2_%s", "klout") )))
    colnames(maqc_a_r1) = c("tx_id","length","eff_length","MAQC_A_kallisto_count_rep1","MAQC_A_kallisto_tpm_rep1")
    colnames(maqc_a_r2) = c("tx_id","length","eff_length","MAQC_A_kallisto_count_rep2","MAQC_A_kallisto_tpm_rep2")
    colnames(maqc_b_r1) = c("tx_id","length","eff_length","MAQC_B_kallisto_count_rep1","MAQC_B_kallisto_tpm_rep1")
    colnames(maqc_b_r2) = c("tx_id","length","eff_length","MAQC_B_kallisto_count_rep2","MAQC_B_kallisto_tpm_rep2")
    maqc <- maqc_a_r1 %>% left_join(., maqc_a_r2[,c("tx_id","MAQC_A_kallisto_tpm_rep2","MAQC_A_kallisto_count_rep2")], by="tx_id") %>%
        left_join(., maqc_b_r1[,c("tx_id","MAQC_B_kallisto_tpm_rep1","MAQC_B_kallisto_count_rep1")], by="tx_id") %>% 
        left_join(., maqc_b_r2[,c("tx_id","MAQC_B_kallisto_tpm_rep2","MAQC_B_kallisto_count_rep2")], by="tx_id")
    MAQC_kallisto <- maqc

# Cufflinks
    maqc_a_r1 <- fread(sprintf("%s/%s", maqc_rnaseq_dir, list.files(maqc_rnaseq_dir, pattern=sprintf("A-1_%s", "clout") )))
    maqc_a_r2 <- fread(sprintf("%s/%s", maqc_rnaseq_dir, list.files(maqc_rnaseq_dir, pattern=sprintf("A-2_%s", "clout") )))
    maqc_b_r1 <- fread(sprintf("%s/%s", maqc_rnaseq_dir, list.files(maqc_rnaseq_dir, pattern=sprintf("B-1_%s", "clout") )))
    maqc_b_r2 <- fread(sprintf("%s/%s", maqc_rnaseq_dir, list.files(maqc_rnaseq_dir, pattern=sprintf("B-2_%s", "clout") )))
    maqc_a_r1 <- maqc_a_r1 %>% dplyr::select(c(4,5,10))
    maqc_a_r2 <- maqc_a_r2 %>% dplyr::select(c(4,10))
    maqc_b_r1 <- maqc_b_r1 %>% dplyr::select(c(4,10))
    maqc_b_r2 <- maqc_b_r2 %>% dplyr::select(c(4,10))
    colnames(maqc_a_r1) = c("gene_id","gene","MAQC_A_cufflinks_fpkm_rep1")
    colnames(maqc_a_r2) = c("gene_id","MAQC_A_cufflinks_fpkm_rep2")
    colnames(maqc_b_r1) = c("gene_id","MAQC_B_cufflinks_fpkm_rep1")
    colnames(maqc_b_r2) = c("gene_id","MAQC_B_cufflinks_fpkm_rep2")
    maqc <- maqc_a_r1 %>% left_join(., maqc_a_r2, by="gene_id") %>% left_join(., maqc_b_r1, by="gene_id") %>% left_join(., maqc_b_r2, by="gene_id")
    MAQC_cufflinks <- maqc

# hg19 ens-gene:ens-tx
    hg19_gtf     <- rtracklayer::readGFF("/storage/references_and_index/hg19/gtf/hg19.gtf")
    hg19_gene_tx <- hg19_gtf[,c("gene_id","transcript_id","gene_name","gene_type")] 
    hg19_gene_tx$ens_geneid <- sapply(hg19_gene_tx$gene_id, function(y) unlist(strsplit(y, "\\."))[1])
    hg19_gene_tx$ens_txid   <- sapply(hg19_gene_tx$transcript_id, function(y) unlist(strsplit(y, "\\."))[1])

# Salmon and Kallisto Gene-leve merge
    MAQC_salmon$ens_geneid <- hg19_gene_tx[match(MAQC_salmon$tx_id, hg19_gene_tx$ens_txid), "ens_geneid"]
    MAQC_salmon_gene <- MAQC_salmon %>% filter( !is.na(ens_geneid) ) %>% group_by(ens_geneid) %>%
        reframe( 
            MAQC_A_salon_count_rep1=sum(MAQC_A_salon_count_rep1),
            MAQC_A_salon_count_rep2=sum(MAQC_A_salon_count_rep2),
            MAQC_B_salon_count_rep1=sum(MAQC_B_salon_count_rep1),
            MAQC_B_salon_count_rep2=sum(MAQC_B_salon_count_rep2),
            EFF_LENGTH = mean(eff_length)
        )

    MAQC_kallisto$ens_geneid <- hg19_gene_tx[match(MAQC_kallisto$tx_id, hg19_gene_tx$ens_txid), "ens_geneid"]
    MAQC_kallisto_gene <- MAQC_kallisto %>% filter( !is.na(ens_geneid) ) %>% group_by(ens_geneid) %>%
        reframe( 
            MAQC_A_kallisto_count_rep1=sum(MAQC_A_kallisto_count_rep1),
            MAQC_A_kallisto_count_rep2=sum(MAQC_A_kallisto_count_rep2),
            MAQC_B_kallisto_count_rep1=sum(MAQC_B_kallisto_count_rep1),
            MAQC_B_kallisto_count_rep2=sum(MAQC_B_kallisto_count_rep2)
        )

# htseq to tpm
    MAQC_salmon_gene = as.data.frame(MAQC_salmon_gene)
    MAQC_HTseq$eff_length = MAQC_salmon_gene[match(MAQC_HTseq$ID, MAQC_salmon_gene$ens_geneid), "EFF_LENGTH"]

    countToTpm <- function(counts, effLen)
    {
        rate <- log(counts) - log(effLen)
        denom <- log(sum(exp(rate)))
        exp(rate - denom + log(1e6))
    }

    MAQC_HTseq <- MAQC_HTseq %>% filter( !is.na(eff_length) )
    MAQC_HTseq$MAQC_A_htseq_tpm_rep1 = countToTpm( counts=MAQC_HTseq$MAQC_A_htseq_count_rep1, effLen=MAQC_HTseq$eff_length )
    MAQC_HTseq$MAQC_A_htseq_tpm_rep2 = countToTpm( counts=MAQC_HTseq$MAQC_A_htseq_count_rep2, effLen=MAQC_HTseq$eff_length )
    MAQC_HTseq$MAQC_B_htseq_tpm_rep1 = countToTpm( counts=MAQC_HTseq$MAQC_B_htseq_count_rep1, effLen=MAQC_HTseq$eff_length )
    MAQC_HTseq$MAQC_B_htseq_tpm_rep2 = countToTpm( counts=MAQC_HTseq$MAQC_B_htseq_count_rep2, effLen=MAQC_HTseq$eff_length )

    GENE_SPACE <- MAQC_HTseq$ID 
    rm_id  = "ENSG00000215549"
    GENE_SPACE = GENE_SPACE[GENE_SPACE != rm_id]

# complete tpm with gene space : salmon 
    MAQC_salmon_gene_rev <- MAQC_salmon_gene[match(GENE_SPACE, MAQC_salmon_gene$ens_geneid), ]
    MAQC_salmon_gene_rev$MAQC_A_salmon_tpm_rep1 = countToTpm( counts=MAQC_salmon_gene_rev$MAQC_A_salon_count_rep1, effLen=MAQC_salmon_gene_rev$EFF_LENGTH )
    MAQC_salmon_gene_rev$MAQC_A_salmon_tpm_rep2 = countToTpm( counts=MAQC_salmon_gene_rev$MAQC_A_salon_count_rep2, effLen=MAQC_salmon_gene_rev$EFF_LENGTH )
    MAQC_salmon_gene_rev$MAQC_B_salmon_tpm_rep1 = countToTpm( counts=MAQC_salmon_gene_rev$MAQC_B_salon_count_rep1, effLen=MAQC_salmon_gene_rev$EFF_LENGTH )
    MAQC_salmon_gene_rev$MAQC_B_salmon_tpm_rep2 = countToTpm( counts=MAQC_salmon_gene_rev$MAQC_B_salon_count_rep2, effLen=MAQC_salmon_gene_rev$EFF_LENGTH )

# complete tpm with gene space : kallisto
    MAQC_kallisto_gene_rev = MAQC_kallisto_gene[match(GENE_SPACE, MAQC_kallisto_gene$ens_geneid), ] 
    MAQC_kallisto_gene_rev$MAQC_A_kallisto_tpm_rep1 = countToTpm( counts=MAQC_kallisto_gene_rev$MAQC_A_kallisto_count_rep1, effLen=MAQC_salmon_gene_rev$EFF_LENGTH )
    MAQC_kallisto_gene_rev$MAQC_A_kallisto_tpm_rep2 = countToTpm( counts=MAQC_kallisto_gene_rev$MAQC_A_kallisto_count_rep2, effLen=MAQC_salmon_gene_rev$EFF_LENGTH )
    MAQC_kallisto_gene_rev$MAQC_B_kallisto_tpm_rep1 = countToTpm( counts=MAQC_kallisto_gene_rev$MAQC_B_kallisto_count_rep1, effLen=MAQC_salmon_gene_rev$EFF_LENGTH )
    MAQC_kallisto_gene_rev$MAQC_B_kallisto_tpm_rep2 = countToTpm( counts=MAQC_kallisto_gene_rev$MAQC_B_kallisto_count_rep2, effLen=MAQC_salmon_gene_rev$EFF_LENGTH )

    
# cufflinks fpkm to tpm
    MAQC_cufflinks_rev <- MAQC_cufflinks[match(GENE_SPACE, MAQC_cufflinks$gene_id), ]
    MAQC_cufflinks$MAQC_A_cufflinks_tpm_rep1 <- MAQC_cufflinks$MAQC_A_cufflinks_fpkm_rep1/sum(MAQC_cufflinks$MAQC_A_cufflinks_fpkm_rep1) * 10^6
    MAQC_cufflinks$MAQC_A_cufflinks_tpm_rep2 <- MAQC_cufflinks$MAQC_A_cufflinks_fpkm_rep2/sum(MAQC_cufflinks$MAQC_A_cufflinks_fpkm_rep2) * 10^6
    MAQC_cufflinks$MAQC_B_cufflinks_tpm_rep1 <- MAQC_cufflinks$MAQC_B_cufflinks_fpkm_rep1/sum(MAQC_cufflinks$MAQC_B_cufflinks_fpkm_rep1) * 10^6
    MAQC_cufflinks$MAQC_B_cufflinks_tpm_rep2 <- MAQC_cufflinks$MAQC_B_cufflinks_fpkm_rep2/sum(MAQC_cufflinks$MAQC_B_cufflinks_fpkm_rep2) * 10^6


# mean tpm
    # htseq
    MAQC_HTseq$MAQC_A = apply(MAQC_HTseq[,c("MAQC_A_htseq_tpm_rep1", "MAQC_A_htseq_tpm_rep2")], 1, mean)
    MAQC_HTseq$MAQC_B = apply(MAQC_HTseq[,c("MAQC_B_htseq_tpm_rep1", "MAQC_B_htseq_tpm_rep2")], 1, mean)
    MAQC_HTseq$MAQC_A = round(log2(MAQC_HTseq$MAQC_A +1),3)
    MAQC_HTseq$MAQC_B = round(log2(MAQC_HTseq$MAQC_B +1),3)

    # salmon
    MAQC_salmon_gene_rev$MAQC_A = apply(MAQC_salmon_gene_rev[,c("MAQC_A_salmon_tpm_rep1", "MAQC_A_salmon_tpm_rep2")], 1, mean)
    MAQC_salmon_gene_rev$MAQC_B = apply(MAQC_salmon_gene_rev[,c("MAQC_B_salmon_tpm_rep1", "MAQC_B_salmon_tpm_rep2")], 1, mean)
    MAQC_salmon_gene_rev$MAQC_A = round(log2(MAQC_salmon_gene_rev$MAQC_A +1),3)
    MAQC_salmon_gene_rev$MAQC_B = round(log2(MAQC_salmon_gene_rev$MAQC_B +1),3)

    # kallisto
    MAQC_kallisto_gene_rev$MAQC_A = apply(MAQC_kallisto_gene_rev[,c("MAQC_A_kallisto_tpm_rep1", "MAQC_A_kallisto_tpm_rep2")], 1, mean)
    MAQC_kallisto_gene_rev$MAQC_B = apply(MAQC_kallisto_gene_rev[,c("MAQC_B_kallisto_tpm_rep1", "MAQC_B_kallisto_tpm_rep2")], 1, mean)
    MAQC_kallisto_gene_rev$MAQC_A = round(log2(MAQC_kallisto_gene_rev$MAQC_A +1),3)
    MAQC_kallisto_gene_rev$MAQC_B = round(log2(MAQC_kallisto_gene_rev$MAQC_B +1),3)
  
    # cufflinks
    MAQC_cufflinks$MAQC_A = apply(MAQC_cufflinks[,c("MAQC_A_cufflinks_tpm_rep1", "MAQC_A_cufflinks_tpm_rep2")], 1, mean)
    MAQC_cufflinks$MAQC_B = apply(MAQC_cufflinks[,c("MAQC_B_cufflinks_tpm_rep1", "MAQC_B_cufflinks_tpm_rep2")], 1, mean)
    MAQC_cufflinks$MAQC_A = round(log2(MAQC_cufflinks$MAQC_A +1),3)
    MAQC_cufflinks$MAQC_B = round(log2(MAQC_cufflinks$MAQC_B +1),3)

    MAQC_RNASEQ <- data.frame( 
        gene = GENE_SPACE, 
        htseq = MAQC_HTseq[match(GENE_SPACE, MAQC_HTseq$ID), c("MAQC_A","MAQC_B")],
        salmon = MAQC_salmon_gene_rev[match(GENE_SPACE, MAQC_salmon_gene_rev$ens_geneid), c("MAQC_A","MAQC_B")],
        kallisto = MAQC_kallisto_gene_rev[match(GENE_SPACE, MAQC_kallisto_gene_rev$ens_geneid), c("MAQC_A","MAQC_B")],
        cufflinks = MAQC_cufflinks[match(GENE_SPACE, MAQC_cufflinks$gene_id), c("MAQC_A","MAQC_B")]
    )

#---------------------------------------------------------------------------------------------------

    saveRDS(MAQC_RNASEQ, file="MAQC_RNAseq_log2TPM_matrix.rds")
    saveRDS(MAQC_RT_PCR, file="MAQC_RT-PCR_ct_values.rds")
    saveRDS(MAQC_ARRAY,  file="MAQC_Microarray_Affymetrix.rds")

#---------------------------------------------------------------------------------------------------
    source("/storage/home/kangsm/runScripts/WTS_Fun.AnalysisModules.R")

    gcx_maqc_seqid =c("MAQC_A_rep1", "MAQC_A_rep2", "MAQC_B_rep1", "MAQC_B_rep2", "MAQC_II_A_rep1","MAQC_II_A_rep2","MAQC_II_A_rep3","MAQC_II_B_rep1","MAQC_II_B_rep2","MAQC_II_B_rep3")
    
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
    EFF_LEN$MAQC_A = apply(EFF_LEN[,c("effLen.MAQC_A_rep1","effLen.MAQC_A_rep2")], 1, function(y) mean(y))
    EFF_LEN$MAQC_B = apply(EFF_LEN[,c("effLen.MAQC_B_rep1","effLen.MAQC_B_rep2")], 1, function(y) mean(y))
    EFF_LEN$MAQC_II_A = apply(EFF_LEN[,c("effLen.MAQC_II_A_rep1","effLen.MAQC_II_A_rep2","effLen.MAQC_II_A_rep3")], 1, function(y) mean(y))
    EFF_LEN$MAQC_II_B = apply(EFF_LEN[,c("effLen.MAQC_II_B_rep1","effLen.MAQC_II_B_rep2","effLen.MAQC_II_B_rep3")], 1, function(y) mean(y))
    EFF_LEN$gene_id = sapply(EFF_LEN$gene_id, function(y) unlist(strsplit(y, "\\."))[1])
    EFF_LEN = EFF_LEN[match(GENE_SPACE, EFF_LEN$gene_id), ]   

    GCX_MAQC_COUNT <- Create.Expression.Values.Matrix(
        SeqIdList=gcx_maqc_seqid,
        ValueType="count",
        DataSource="rds"
    )
    rownames(GCX_MAQC_COUNT) = sapply( rownames(GCX_MAQC_COUNT), function(y) unlist(strsplit(y, "\\."))[1] )

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

    ggsave(corr_barplot1, file = "maqc_corr_barplot_type1_rev.png", width=12, height=8, unit='in', dpi=150)










