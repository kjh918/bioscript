

cnvkit_nt = read.delim("/data/wes/WES_24_02/WES_24_02_09/cnv/cnvkit/WES_24_02_09.NT.CNVkit.GeneLevel.CopyNumberAlterations.txt")
cnvkit_tonly = read.delim("/data/wes/WES_24_02/WES_24_02_09/cnv/cnvkit/WES_24_02_09.Tonly.CNVkit.GeneLevel.CopyNumberAlterations.txt")
ascat = read.delim("/data/wes/WES_24_02/WES_24_02_09/cnv/ascat/WES_24_02_09.NT.result/WES_24_02_09.NT.ASCAT.gene.level.CNA.raw.tsv")

# ascat =1 
# cnvkit nt = 2
# cnvkit tonly =3

genespace = unique(c(ascat$GENE, cnvkit_nt$gene, cnvkit_tonly$gene))

res = data.frame(
    gene = genespace,
    cn_1 = ascat[match(genespace, ascat$GENE), "CN"],
    cn_2 = cnvkit_nt[match(genespace, cnvkit_nt$gene), "cn"],
    cn_3 = cnvkit_tonly[match(genespace, cnvkit_tonly$gene), "cn"],
    cna_abs_1 = ascat[match(genespace, ascat$GENE), "CNA_ABS"],
    cna_abs_2 = cnvkit_nt[match(genespace, cnvkit_nt$gene), "absolute_cna"],
    cna_abs_3 = cnvkit_tonly[match(genespace, cnvkit_tonly$gene), "absolute_cna"],
    cna_pld_1 = ascat[match(genespace, ascat$GENE), "CNA_PLOIDY"],
    cna_pld_2 = cnvkit_nt[match(genespace, cnvkit_nt$gene), "ploidy_cna"],
    cna_pld_3 = cnvkit_tonly[match(genespace, cnvkit_tonly$gene), "ploidy_cna"]
)

RES = res[complete.cases(res), ]

ggplot(RES) +
    geom_point( aes( x= cn_1, y=cn_2))


cor()




