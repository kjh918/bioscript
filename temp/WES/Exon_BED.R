
#--------------------------------------------------------------------------------------------------#
options(stringsAsFactors=FALSE) 
library(plyr)
library(dplyr)

GFF_DIR <- "/storage/home/kangsm/myDB/gtf/hg19"


#---| UCSC hg19 TxDb |-----------------------------------------------------------------------------#
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb.ucsc.hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene
save(txdb.ucsc.hg19, file=sprintf("%s/hg19_ucsc.TxDb.Rdata", GFF_DIR))
#------------------------------------------------------------------------------#
exonic.only.ucsc.hg19 <- exonicParts(txdb.ucsc.hg19, linked.to.single.gene.only = FALSE)
save(txdb.ensembl.gencode.v19, exonic.only.ucsc.hg19, file=sprintf("%s/hg19_ucsc.TxDb_include.exonic.only.Rdata", GFF_DIR))
#--------------------------------------------------------------------------------------------------#

#---| ENSEMBL (GENCODE V19) TxDb |-----------------------------------------------------------------#
library(GenomicFeatures)
txdb.ensembl.gencode.v19 <- makeTxDbFromGFF(sprintf("%s/gencode.v19.annotation.gff3.gz",GFF_DIR))
save(txdb.ensembl.gencode.v19, file=sprintf("%s/GRCh37_ensembl.gencode.v19.TxDb.Rdata", GFF_DIR))
#------------------------------------------------------------------------------#
exonic.only.ensembl.gencode.v19 <- exonicParts(txdb.ensembl.gencode.v19, linked.to.single.gene.only = FALSE)
save(txdb.ensembl.gencode.v19, exonic.only.ensembl.gencode.v19, file=sprintf("%s/GRCh37_ensembl.gencode.v19.TxDb_include.exonic.only.Rdata", GFF_DIR))
#--------------------------------------------------------------------------------------------------#

#---| TMB HARMONIZATION PROJECT PHASE I : EXONIC REGION REFERENCE |--------------------------------#
# BASH >> curl 'https://ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/15/CCDS.20131129.txt' -o CCDS.20131129.txt
#------------------------------------------------------------------------------#
# BASH >> cp CCDS.20131129.txt TMB_REF_exonic_region.txt
#------------------------------------------------------------------------------#
exonRegion <- read.delim("TMB_REF_exonic_region.txt") %>% filter( ccds_status == "Public" )

INDEXING = function(X) { sapply(unique(X), function(y) list(y)) }
idx = INDEXING(1:nrow(exonRegion))

exon.region.ref <- mclapply( idx, function(N) {
    ER  <- exonRegion[N, ]
    cds <- ER$cds_location
    cds <- gsub("^\\[", "",  cds)
    cds <- gsub("]$",   "",  cds)
    cds <- gsub(", ",   ";", cds)
    CDS <- lapply(unlist(strsplit(cds, ";")), function(gr) {
        GR <- unlist(strsplit(gr, "\\-"))
        res = data.frame( start = GR[1], end = GR[2] )
        return(res)
    })
    res = data.frame(ER[,-10], ldply(CDS))
    return(res)
}, mc.cores = 30)

EXON_REGION_REF = ldply(exon.region.ref)
EXON_REGION_REF$start = as.numeric(EXON_REGION_REF$start)
EXON_REGION_REF$end = as.numeric(EXON_REGION_REF$end)
EXON_REGION_REF$cds_from = as.numeric(EXON_REGION_REF$cds_from)
EXON_REGION_REF$cds_to = as.numeric(EXON_REGION_REF$cds_to)
head(EXON_REGION_REF)

EXON_REGION_REF$source = "CCDS"
EXON_REGION_REF$feqture_type = "CDS"
EXON_REGION_REF$score = "."
EXON_REGION_REF$phase = "."
EXON_REGION_REF$attributes <- apply( EXON_REGION_REF[,c("gene","gene_id")], 1, function(y) paste0(y[1], as.numeric(y[2]), sep=";"))

EXON_REGION_REF %>% filter( chromosome == 1 ) %>% arrange(start)
EXON_REGION_REF %>% group_by(chromosome) %>% summarise( ERL=sum(end-start)/1000000 )

EXON_GFF <- EXON_REGION_REF %>% select(c("chromosome","source","feqture_type","start","end","score","cds_strand","phase","attributes"))

CHR_TAG=c(
"##sequence-region 1 1 249250621",
"##sequence-region 2 1 243199373",
"##sequence-region 3 1 198022430",
"##sequence-region 4 1 191154276",
"##sequence-region 5 1 180915260",
"##sequence-region 6 1 171115067",
"##sequence-region 7 1 159138663",
"##sequence-region 8 1 146364022",
"##sequence-region 9 1 141213431",
"##sequence-region 10 1 135534747",
"##sequence-region 11 1 135006516",
"##sequence-region 12 1 133851895",
"##sequence-region 13 1 115169878",
"##sequence-region 14 1 107349540",
"##sequence-region 15 1 102531392",
"##sequence-region 16 1 90354753",
"##sequence-region 17 1 81195210",
"##sequence-region 18 1 78077248",
"##sequence-region 19 1 59128983",
"##sequence-region 20 1 63025520",
"##sequence-region 21 1 48129895",
"##sequence-region 22 1 51304566",
"##sequence-region X 1 155270560",
"##sequence-region Y 1 59373566")

write("#description: TMB Harmonization Project Phase I : Exonic Region (GRCh37)", file="GRCh37_tmb_exonic_region.gff", append=T)
write("#provider: CCDS (https://ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/15/CCDS.20131129.txt)" , file="GRCh37_tmb_exonic_region.gff", append=T)
write("#contact: knagsm@gencurix.com", , file="GRCh37_tmb_exonic_region.gff", append=T)
write("#format: gff3", file="GRCh37_tmb_exonic_region.gff", append=T)
write("#date: 2024-03-15" , file="GRCh37_tmb_exonic_region.gff", append=T)

CHRS <- c(1:22, "X","Y")

for( N in c(1:24) )
{
    write(CHR_TAG[N], file="GRCh37_tmb_exonic_region.gff", append=T)
    chr_exon = EXON_GFF %>% filter( chromosome == CHRS[N] ) %>% arrange(start)
    chr_exon$attributes = gsub("^", "GENE=", chr_exon$attributes)
    chr_exon$attributes = gsub(";", ";ENTREZ=", chr_exon$attributes)
    chr_exon$start = as.character(chr_exon$start)
    chr_exon$end = as.character(chr_exon$end)
    write.table( chr_exon, file="GRCh37_tmb_exonic_region.gff", quote=F, col.names=F, row.names=F, sep="\t", append=T)
    print(N)
}

#---| CCDS TxDb |----------------------------------------------------------------------------------#
library(GenomicFeatures)
GFF_DIR = "/storage/home/kangsm/analysis/exon_bed"
txdb.tmb.ccds.grch37 <- makeTxDbFromGFF(sprintf("%s/GRCh37_tmb_exonic_region.gff",GFF_DIR))
save(txdb.ensembl.gencode.v19, file=sprintf("%s/GRCh37_ensembl.gencode.v19.TxDb.Rdata", GFF_DIR))
#------------------------------------------------------------------------------#
exonic.only.ensembl.gencode.v19 <- exonicParts(txdb.ensembl.gencode.v19, linked.to.single.gene.only = FALSE)
save(txdb.ensembl.gencode.v19, exonic.only.ensembl.gencode.v19, file=sprintf("%s/GRCh37_ensembl.gencode.v19.TxDb_include.exonic.only.Rdata", GFF_DIR))
#--------------------------------------------------------------------------------------------------#
colnames(EXON_GFF)[3] = "type"
colnames(EXON_GFF)[7] = "strand"
EXON_GFF$type = "exon"
EXON_GFF$phase = 0
EXON_GR <- makeGRangesFromDataFrame(EXON_GFF[,c(1,3,4,5,6,7,9)], keep.extra.columns=T)
txdb.tmb.ccds.grch37 <- makeTxDbFromGRanges(EXON_GR)
exonic.only.tmb.ccds.grch37 <- exonicParts(txdb.tmb.ccds.grch37, linked.to.single.gene.only = FALSE)


for( N in c(1:24) )
{
    print(CHR_TAG[N])
    chr_exon = EXON_GFF %>% filter( chromosome == CHRS[N] ) %>% arrange(start)
    print(chr_exon$start)
    print(chr_exon$end)
    
    print(N)
}


W = EXON_REGION_REF$end - EXON_REGION_REF$start

