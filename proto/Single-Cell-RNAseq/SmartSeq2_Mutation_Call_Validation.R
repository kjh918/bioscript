

library(plyr)
library(dplyr)
library(data.table)
library("BSgenome.Hsapiens.UCSC.hg38")
GENOME <- BSgenome.Hsapiens.UCSC.hg38
indexing = function(X) { sapply(unique(X), function(y) list(y)) }

sample_list=fread("/data/scRNAseq/02.Public-Datasets/004.GSE151334.PRJNA635600.cell-line-smartseq2/meta/variant.call.SmartSeq2.sample.list")
VCF_DIR="/data/scRNAseq/02.Public-Datasets/004.GSE151334.PRJNA635600.cell-line-smartseq2/analysis/SmartSeq2_variant_call/vcf"

mut_hek293 = fread("/data/scRNAseq/02.Public-Datasets/004.GSE151334.PRJNA635600.cell-line-smartseq2/meta/HEK293_mutations.csv")
mut_mcf7   = fread("/data/scRNAseq/02.Public-Datasets/004.GSE151334.PRJNA635600.cell-line-smartseq2/meta/MCF7_mutations.csv")

sce = readRDS("/data/scRNAseq/02.Public-Datasets/004.GSE151334.PRJNA635600.cell-line-smartseq2/analysis/SmartSeq2/smartseq2_analyzed.SCE.objects.RDS")

pbcounts <- readRDS(sprintf("/data/scRNAseq/02.Public-Datasets/004.GSE151334.PRJNA635600.cell-line-smartseq2/analysis/SmartSeq2/smartseq2_pseudo.bulk.counts.df.RDS"))

maf_list = lapply(indexing(sample_list$IndivID), function(idx) {
    maf = fread(sprintf("%s/%s.mutect2.Tonly.filtered.norm.vep.annotated.maf", VCF_DIR, idx), skip=1) 
    #maf = maf %>% filter( FILTER == "PASS")
    maf$REF <- apply(maf[,c("Chromosome","vcf_pos","Variant_Type","Reference_Allele")], 1, function(v) {
        if( v[3] == "INS" ){
            ref_allele = as.character(GENOME[[v[1]]][as.numeric(v[2]):as.numeric(v[2])])
        }else if( v[3] == "DEL" ){
            ref_allele = paste0(as.character(GENOME[[v[1]]][as.numeric(v[2]):as.numeric(v[2])]), v[4])
        }else{
            ref_allele = v[4]
        }
        return(ref_allele)
    })
    maf$ALT <- apply(maf[,c("Chromosome","vcf_pos","Variant_Type","Tumor_Seq_Allele2")], 1, function(v) {
        if( v[3] == "INS" ){
            alt_allele = paste0(as.character(GENOME[[v[1]]][as.numeric(v[2]):as.numeric(v[2])]), v[4])
        }else if( v[3] == "DEL" ){
            alt_allele = as.character(GENOME[[v[1]]][as.numeric(v[2]):as.numeric(v[2])])
        }else{
            alt_allele = v[4]
        }
        return(alt_allele)
    })
    return(maf)
})


MCF7_VARS <- mut_mcf7 %>% dplyr::select(c(
    "GENE"="Gene","CHROM"="Chromosome","POS"="Position","REF"="Ref Allele","ALT"="Alt Allele",
    "HGVSp"="Protein Change","DP"="DP","VAF"="Allele Fraction"
)) %>% as.data.frame()
MCF7_VARS_VALID = data.frame(MCF7_VARS, expr=pbcounts[match(MCF7_VARS$GENE, pbcounts$gene), -1])
MCF7_VARS_VALID = MCF7_VARS_VALID %>% mutate(index=as.character(1:nrow(.)))

mcf7_res_valid <- data.frame()
for( i in 1:nrow(MCF7_VARS_VALID) ){

    data_dp <- lapply(maf_list, function(y) {
        mut <- y %>% filter(Chromosome==MCF7_VARS_VALID[i, "CHROM"], vcf_pos==MCF7_VARS_VALID[i, "POS"] ) %>% as.data.frame()
        if( nrow(mut) == 0 ) {
            dp = 0
        }else{
            dp = mut$DP
        }
        return(dp)
    }) %>% ldply() %>% dplyr::rename(cell=1, dp=2) %>% mutate(group="dp")

    data_vaf <- lapply(maf_list, function(y) {
        mut <- y %>% filter(Chromosome==MCF7_VARS_VALID[i, "CHROM"], vcf_pos==MCF7_VARS_VALID[i, "POS"] ) %>% as.data.frame()
        if( nrow(mut) == 0 ) {
            vaf = 0
        }else{
            vaf = round(mut$t_alt_count/(mut$t_alt_count+mut$t_ref_count), 3)
        }
        return(vaf)
    }) %>% ldply() %>% dplyr::rename(cell=1, vaf=2) %>% mutate(group="vaf")

    data_dna <- lapply(maf_list, function(y) {
        mut <- y %>% filter(Chromosome==MCF7_VARS_VALID[i, "CHROM"], vcf_pos==MCF7_VARS_VALID[i, "POS"] ) %>% as.data.frame()
        if( nrow(mut) == 0 ) {
            dna = ""
        }else{
            dna = paste0( mut$REF ,">", mut$ALT )
        }
        return(dna)
    }) %>% ldply() %>% dplyr::rename(cell=1, change=2) %>% mutate(group="change")

    res <- cbind(
        data.frame(index=MCF7_VARS_VALID[i, "index"]),
        reshape(data_dp, idvar='group', timevar='cell', direction="wide")[,-1],
        reshape(data_vaf, idvar='group', timevar='cell', direction="wide")[,-1],
        reshape(data_dna, idvar='group', timevar='cell', direction="wide")[,-1]
    ) 
    mcf7_res_valid = rbind(mcf7_res_valid, res)
}
MCF7_VARS_VALID = MCF7_VARS_VALID %>% left_join(mcf7_res_valid, by="index")
MCF7_VARS_VALID = MCF7_VARS_VALID[, 
c("GENE","CHROM","POS","REF","ALT","HGVSp","DP","VAF", 
"expr.MCF7_SE","dp.MCF7_SE","change.MCF7_SE","vaf.MCF7_SE",
"expr.HEK293_PE","dp.HEK293_PE","change.HEK293_PE","vaf.HEK293_PE",
"expr.HEK293_SE","dp.HEK293_SE","change.HEK293_SE","vaf.HEK293_SE",
"expr.HDFn_SE","dp.HDFn_SE","change.HDFn_SE","vaf.HDFn_SE"            
)
]

MCF7_VARS_VALID %>% filter( dp.MCF7_SE > 0 )


HEK293_VARS <- mut_hek293 %>% dplyr::select(c(
    "GENE"="Gene","CHROM"="Chromosome","POS"="Position","REF"="Ref Allele","ALT"="Alt Allele",
    "HGVSp"="Protein Change","DP"="DP","VAF"="Allele Fraction"
)) %>% as.data.frame()

HEK293_VARS_VALID = data.frame(HEK293_VARS, expr=pbcounts[match(HEK293_VARS$GENE, pbcounts$gene), -1])
HEK293_VARS_VALID = HEK293_VARS_VALID %>% mutate(index=as.character(1:nrow(.)))
hek293_res_valid <- data.frame()
for( i in 1:nrow(HEK293_VARS_VALID) ){

    data_dp <- lapply(maf_list, function(y) {
        mut <- y %>% filter(Chromosome==HEK293_VARS_VALID[i, "CHROM"], vcf_pos==HEK293_VARS_VALID[i, "POS"] ) %>% as.data.frame()
        if( nrow(mut) == 0 ) {
            dp = 0
        }else{
            dp = mut$DP
        }
        return(dp)
    }) %>% ldply() %>% dplyr::rename(cell=1, dp=2) %>% mutate(group="dp")

    data_vaf <- lapply(maf_list, function(y) {
        mut <- y %>% filter(Chromosome==HEK293_VARS_VALID[i, "CHROM"], vcf_pos==HEK293_VARS_VALID[i, "POS"] ) %>% as.data.frame()
        if( nrow(mut) == 0 ) {
            vaf = 0
        }else{
            vaf = round(mut$t_alt_count/(mut$t_alt_count+mut$t_ref_count), 3)
        }
        return(vaf)
    }) %>% ldply() %>% dplyr::rename(cell=1, vaf=2) %>% mutate(group="vaf")

    data_dna <- lapply(maf_list, function(y) {
        mut <- y %>% filter(Chromosome==HEK293_VARS_VALID[i, "CHROM"], vcf_pos==HEK293_VARS_VALID[i, "POS"] ) %>% as.data.frame()
        if( nrow(mut) == 0 ) {
            dna = ""
        }else{
            dna = paste0( mut$REF ,">", mut$ALT )
        }
        return(dna)
    }) %>% ldply() %>% dplyr::rename(cell=1, change=2) %>% mutate(group="change")

    res <- cbind(
        data.frame(index=HEK293_VARS_VALID[i, "index"]),
        reshape(data_dp, idvar='group', timevar='cell', direction="wide")[,-1],
        reshape(data_vaf, idvar='group', timevar='cell', direction="wide")[,-1],
        reshape(data_dna, idvar='group', timevar='cell', direction="wide")[,-1]
    ) 
    hek293_res_valid = rbind(hek293_res_valid, res)
}
HEK293_VARS_VALID = HEK293_VARS_VALID %>% left_join(hek293_res_valid, by="index")
HEK293_VARS_VALID = HEK293_VARS_VALID[, 
c("GENE","CHROM","POS","REF","ALT","HGVSp","DP","VAF", 
"expr.HEK293_PE","dp.HEK293_PE","change.HEK293_PE","vaf.HEK293_PE",
"expr.HEK293_SE","dp.HEK293_SE","change.HEK293_SE","vaf.HEK293_SE",
"expr.MCF7_SE","dp.MCF7_SE","change.MCF7_SE","vaf.MCF7_SE",
"expr.HDFn_SE","dp.HDFn_SE","change.HDFn_SE","vaf.HDFn_SE"            
)
]
HEK293_VARS_VALID %>% filter( dp.HEK293_PE > 0 | dp.HEK293_SE > 0 )

openxlsx::write.xlsx(MCF7_VARS_VALID, "SmartSeq2_variant.call.mutect2.validation.mcf7.xlsx")
openxlsx::write.xlsx(HEK293_VARS_VALID, "SmartSeq2_variant.call.mutect2.validation.hek293.xlsx")



valid_data = openxlsx::read.xlsx("SmartSeq2_variant.call.mutect2.validation.hek293.xlsx")
valid_data[which(is.na(valid_data$expr.MCF7_SE)),"expr.HEK293_PE"] = 0
valid_data[which(is.na(valid_data$dp.HEK293_PE)),"dp.HEK293_PE"] = 0
valid_data[which(is.na(valid_data$expr.HEK293_SE)),"expr.HEK293_SE"] = 0
valid_data[which(is.na(valid_data$dp.HEK293_SE)),"dp.HEK293_SE"] = 0

pdata1 = valid_data[,c('expr.HEK293_PE','dp.HEK293_PE','change.HEK293_PE')] %>% mutate(
    log_expr = log10(expr.HEK293_PE+1),
    log_dp   = log10(dp.HEK293_PE+1),
    colors   = case_when(
        change.HEK293_PE == "" ~ "#adadad",
        TRUE ~ "#C73500"
    )
) 
pdata2 = valid_data[,c('expr.HEK293_SE','dp.HEK293_SE','change.HEK293_SE')] %>% mutate(
    log_expr = log10(expr.HEK293_SE+1),
    log_dp   = log10(dp.HEK293_SE+1),
    colors   = case_when(
        change.HEK293_SE == "" ~ "#adadad",
        TRUE ~ "#C73500"
    )
) 

p_mut1 = ggplot(pdata1) + 
    geom_point(aes(x=log_expr, y=log_dp), col=pdata1$colors, size=3, pch=19) + 
    labs(x="log10 Expr-Counts", y="log10 DP") +    
    theme_light() +
    theme(
        axis.text  = element_text(face="bold", size=10, family="Roboto"),
        axis.title = element_text(face="bold", size=14)
 
    )
p_mut2 = ggplot(pdata2) + 
    geom_point(aes(x=log_expr, y=log_dp), col=pdata2$colors, size=3, pch=19) + 
    labs(x="log10 Expr-Counts", y="log10 DP") +    
    theme_light() +
    theme(
        axis.text  = element_text(face="bold", size=10, family="Roboto"),
        axis.title = element_text(face="bold", size=14)
 
    )
p_mut_list = cowplot::plot_grid(p_mut1, p_mut2, nrow=1)
ggsave(p_mut_list, file="hek293_mut_expr_dp.png", width=18, height=7, unit='in', dpi=150 )



