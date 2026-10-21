
library(RMySQL)
library(Hmisc)
library(plyr)
library(dplyr)


dbCon = dbConnect(dbDriver("MySQL"), host='192.168.0.35', user='gcx', port=3306, password='gencurix!!', db='gcx_common' )
ginfo = dbGetQuery(dbCon, "SELECT * FROM gene_info_hgnc")
dbDisconnect(dbCon)

GENOME_ASSEMBLY="hg38"

refFlat = read.delim(sprintf("/storage/home/kangsm/myDB/annotation/%s/%s.refFlat.txt",GENOME_ASSEMBLY,GENOME_ASSEMBLY), header=F)

refGene = refFlat %>% rename(gene=1, tx=2, chr=3, strand=4) %>%
    group_by(gene, chr) %>% 
    mutate( refseq = paste(tx, collapse=";")) %>% 
    select(c("gene","chr","strand","refseq")) %>%
    unique()

refGene <- data.frame(refGene, ginfo[match( refGene$gene, ginfo$input_symbol), c(5,6,7,8,9,10,11,4,2)] )

table(refGene$locus_group)
table(refGene$locus_type)
refGene %>% filter( locus_group == "protein-coding gene") %>% .$locus_type %>% table()

CHROMOSOMES = c(1:22, "X","Y")
if( GENOME_ASSEMBLY == "hg38" ){ CHROMOSOMES = paste0( "chr", CHROMOSOMES) }

cnvRefGene <- refGene %>% 
    filter( chr %in% CHROMOSOMES ) %>% 
    filter( locus_group %in% c("protein-coding gene", "other") ) %>% 
    select(c(1,2,3,5,6,9,8,10,11,12,13,4)) %>% 
    unique()


CNV_REF_GENES = list()
CNV_REF_GENES$hg38=cnvRefGene
CNV_REF_GENES$hg19=cnvRefGene

save(CNV_REF_GENES, file= "/storage/home/kangsm/myDB/rds/cnv_refGenes_list_hg19_hg38.Rdata")

saveRDS(CNV_REF_GENES, file= "/storage/home/kangsm/myDB/rds/cnv_refGenes_list_hg19_hg38.RDS")
