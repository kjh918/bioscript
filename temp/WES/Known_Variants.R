
## script for checking known variants site
# 1) pathogenic germline variants
# 2) cancer hot spot

library(RMySQL)
dbCon = dbConnect(dbDriver("MySQL"), host='192.168.0.35', user='gcx', port=3306, password='gencurix!!', db='gcx_common' )
PGV   = dbGetQuery(dbCon, "SELECT * FROM gcx_pathogenic_germline_variants")
CHS   = dbGetQuery(dbCon, "SELECT * FROM cancer_hotspot_variants_hg19")
dbDisconnect(dbCon)

GR_CancerHotspots <- GRanges(
    seqnames = CHS$chr,
    ranges   = IRanges( CHS$position, width=1),
    gene     = CHS$gene,
    variant_class = CHS$variant_class,
    changeDNA = CHS$change_DNA,
    changeProtein = CHS$change_Protein,
    assembly = "hg19"
)
save( GR_CancerHotspots, file="/storage/home/kangsm/myDB/rds/hg19_cancerhotspot_variants.Rdata")

PGV = PGV %>% mutate( change_DNA = paste(ref, ">", alt, spe=" "))

GR_PathogenicGermVars1 = GRanges(
    PGV$chr,
    IRanges(PGV$start, width=1),
    gene = PGV$hugo_symbol,
    variant_class = PGV$variant_classification,
    changeDNA = PGV$change_DNA,
    changeProtein = PGV$HGVSp_short,
    level = PGV$level,
    TCGA_cancer = PGV$cancer,
    tissue = PGV$tissue,
    gene_class = PGV$gene_class
)

save( GR_PathogenicGermVars1, file="/storage/home/kangsm/myDB/rds/hg19_pathogenic_germline_variants_01.Rdata")


GR_PathogenicGermVars2 = GRanges(
    cgv$Chr,
    IRanges( cgv$Position, width=1),
    gene = cgv$Gene,
    variant_class = cgv$type,
    changeDNA = cgv$change_DNA,
    eastasia = cgv$MutationInEastAsia,
    clinvar = cgv$ClinVar_ID,
    disease = cgv$disease

)
save( GR_PathogenicGermVars2, file="/storage/home/kangsm/myDB/rds/hg19_pathogenic_germline_variants_02.Rdata")

pgv1 = as.data.frame(GR_PathogenicGermVars1)
pgv2 = as.data.frame(GR_PathogenicGermVars2)
pgv3 = pgv2[-ovlp$subjectHits, ]

pgv3$disease = gsub("\\\n", ", ", pgv3$disease)


ovlp = as.data.frame(findOverlaps(GR_PathogenicGermVars1, GR_PathogenicGermVars2))

GR_PathogenicGermVars2[ovlp$subjectHits, ]


findOverlaps(CHS.gr,GVAR.gr)

subsetByOverlaps(GVAR.gr,CHS.gr)
subsetByOverlaps(SVAR.gr,CHS.gr)

subsetByOverlaps(SVAR.gr,PGV.gr)
subsetByOverlaps(SVAR.gr,CGV.gr)
subsetByOverlaps(PGV.gr,CGV.gr)

findOverlaps(PGV.gr,CGV.gr)

#------------------------------------------------------------------------------#

cgv <- read.csv("~/analysis/temp/cancer_genes/MaGiC-Data.csv")

cgv2 = cgv %>% filter( disease != "." )
cgv2$Ref = gsub("^'", "", cgv2$Ref)
cgv2$Alt = gsub("^'", "", cgv2$Alt)
cgv2$Molecular_Consequence = gsub("^'", "", cgv2$Molecular_Consequence)
cgv2$ClinVar_ID = gsub("^'", "", cgv2$ClinVar_ID)

cgv3 = mclapply(indexing(1:nrow(cgv2)), function(y) 
{
    gv.y = cgv2[y, ]
    splitchanges = unlist(strsplit(gsub("\\\n", ";", gv.y$Change), ";"))
    ds = gsub("\\\n", ", ", gv.y$disease)
    GV = data.frame( gv.y[,c(2,3,5,6,7,9,10,11,12)], changes=splitchanges, diseases=ds)

    GV$tx = sapply(GV$changes, function(z) unlist(strsplit(z, "\\:"))[2])
    GV$exon = sapply(GV$changes, function(z) unlist(strsplit(z, "\\:"))[3])
    GV$HGVSc = sapply(GV$changes, function(z) unlist(strsplit(z, "\\:"))[4])
    GV$HGVSp = sapply(GV$changes, function(z) unlist(strsplit(z, "\\:"))[5])
    
    return(GV)
}, mc.cores = 25)

refseq = read.delim("/storage/DB/ucsc_public_data/hg19/ncbiRefSeqSelect.txt", header=F)
refseqtb = refseq %>% select(c(2,3,13)) %>% dplyr::rename( tx=1, chr=2, gene=3)
refseqtb$refseq = sapply( refseqtb$tx, function(y)  unlist(strsplit(y, "\\."))[1] )

cgv4 = mclapply( cgv3, function(y)
{
    y$refseq = sapply(y$tx, function(z) unlist(strsplit(z, "\\."))[1] )
    y2 = y[which(y$refseq %in% refseqtb$refseq), ]
    if( nrow(y2) == 0 ){ y2 = y }
    return(y2)
}, mc.cores=25)

cgv5 = mclapply( cgv4, function(y) 
{
    y2 =unique(y[,c(1,2,3,4,5,6,7,8,9,11,14,15)])
    return(y2)
}, mc.cores=25)

CGV = ldply(cgv5)

CGV = CGV %>% mutate( change_DNA = paste(Ref, ">", Alt, sep= " "))

GR_PathogenicGermVars3 = GRanges(
    CGV$Chr,
    IRanges( CGV$Position, width=1),
    gene = CGV$Gene,
    variant_class = CGV$type,
    changeDNA = CGV$change_DNA,
    changeProtein = CGV$HGVSp,
    eastasia = CGV$MutationInEastAsia,
    Consequence = CGV$Molecular_Consequence,
    disease = CGV$diseases
)

ovlp = as.data.frame(findOverlaps(GR_PathogenicGermVars1, GR_PathogenicGermVars3))
pgv1 = as.data.frame(GR_PathogenicGermVars1)
pgv2 = as.data.frame(GR_PathogenicGermVars3)
pgv3 = pgv2[-ovlp$subjectHits, ]

pgv3$disease = gsub("\\\n", ", ", pgv3$disease)

openxlsx::write.xlsx(pgv3, "pathgenic_germline_vars_db2.xlsx", overwrite=T, rowNames=F)

CGV_rev = openxlsx::read.xlsx("pathgenic_germline_vars_db2_rev.xlsx")

selCols = c("seqnames","start","gene","variant_class","changeDNA","changeProtein","disease","Lv","DB")

PGV = rbind(
    pgv1 %>% dplyr::rename( disease = TCGA_cancer, Lv=level ) %>% mutate( DB = "TCGA_PathogenicGermlineVariants") %>% select(any_of(selCols)),
    CGV_rev %>% dplyr::rename( Lv = eastasia ) %>% mutate( DB = "COGVIC") %>% select(any_of(selCols))
) 

GR_PathogenicGermVars <-  GRanges(
    PGV$seqnames,
    IRanges( PGV$start, width=1),
    gene = PGV$gene,
    variant_class = PGV$variant_class,
    changeDNA = PGV$changeDNA,
    changeProtein = PGV$changeProtein,
    disease = PGV$disease,
    level = PGV$Lv,
    db_source = PGV$DB
)

save(GR_PathogenicGermVars, file="/storage/home/kangsm/myDB/rds/hg19_Pathogenic.Germline.Variants.Rdata")














