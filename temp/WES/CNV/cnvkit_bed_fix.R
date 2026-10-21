
library(RMySQL)
library(dplyr)
dbCon = dbConnect(dbDriver("MySQL"), host='192.168.0.35', user='gcx', port=3306, password='gencurix!!', db='gcx_common' )
ginfo = dbGetQuery(dbCon, "SELECT * FROM gene_info_hgnc")
dbDisconnect(dbCon)

GENOME_ASSEMBLY="hg19"

REF_FLAT=sprintf("/storage/home/kangsm/myDB/annotation/%s/%s.refFlat.txt", GENOME_ASSEMBLY, GENOME_ASSEMBLY)
REF = read.delim(REF_FLAT, header=F)

###  sureselect.v6 |-------------------------------------------------------#
    GENOME_ASSEMBLY="hg19"
    WES_LIB_KIT="sureselect.v6"
    WES_BAIT_BED=sprintf("/storage/references_and_index/%s/bed/%s/%s_%s.probe.bed", GENOME_ASSEMBLY, WES_LIB_KIT, GENOME_ASSEMBLY, WES_LIB_KIT )
    BED = read.delim(WES_BAIT_BED, header=F)

    bed_gene = mclapply(BED[,4], function(rf) {
        tags = unlist(strsplit(rf, ",")) 
        ref_tag = grep("ref", tags)
        REF_TAG = tags[ref_tag]
        if( length(ref_tag) == 0 ) { GENE = "" }else{
            ref_gene = REF_TAG[!grepl("NM_|NR_|XM_|XR_|ClinID|ENST|ENSG|Intron|CCDS", REF_TAG)]
            hgnc_gene = ginfo[which(ginfo$input_symbol %in% unique(gsub("ref\\|", "", ref_gene))), "symbol"]
            REF_GENE = unique(gsub("ref\\|", "", ref_gene))
            HGNC_GENE = unique(hgnc_gene)
            GENE = paste(unique(REF_GENE,HGNC_GENE), collapse=",")
        }
        return(GENE)
    }, mc.cores = 30 )

    BED_gene = BED
    BED_gene[,4] = unlist(bed_gene)

    write.table(BED_gene, 
        "/storage/references_and_index/hg19/bed/sureselect.v6/hg19_sureselect.v6.probe.cnvkit.gene.bed",
        quote=F, col.names=F, row.names=F, sep="\t"
        )

#------------------------------------------------------------------------------#

###  twist.exome.2.0 |-----------------------------------------------------#
    GENOME_ASSEMBLY="hg19"
    WES_LIB_KIT="twist.exome.2.0"
    WES_BAIT_BED=sprintf("/storage/references_and_index/%s/bed/%s/%s_%s.probe.bed", GENOME_ASSEMBLY, WES_LIB_KIT, GENOME_ASSEMBLY, WES_LIB_KIT )
    BED = read.delim(WES_BAIT_BED, header=F)

    bed_gene = mclapply( BED[,4], function(bd) 
    {
        tags = unlist(strsplit(bd, ";"))
        ensg = tags[grepl("^ENSG", tags)]
        hgnc_ens = ginfo[which(ginfo$ens_geneid %in% unique(ensg)), "symbol"]
        remainTag  = tags[!grepl("^NM_|^NR_|^XM_|^XR_|^ClinID|^ENST|^ENSG|^Intron|^CCDS", tags)]
        remainHgnc = ginfo[which(ginfo$input_symbol %in% remainTag), "symbol"]
        tagGenes = paste( unique(c(remainTag,remainHgnc,hgnc_ens)) , collapse=",")
        #tagGenes = paste( unique(hgnc_ens) , collapse=",")
        return(tagGenes)
    }, mc.cores = 30 )

    BED_gene = BED
    BED_gene[,4] = unlist(bed_gene)

    write.table(BED_gene, 
        "/storage/references_and_index/hg19/bed/twist.exome.2.0/hg19_twist.exome.2.0.probe.cnvkit.gene.bed",
        quote=F, col.names=F, row.names=F, sep="\t"
        )
#------------------------------------------------------------------------------#

