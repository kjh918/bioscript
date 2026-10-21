
library(Cairo)
grDevices::X11.options(type='cairo')
options(device='x11')
options(stringsAsFactors=FALSE) 

library(ggplot2)
library(reshape2)
library(dplyr)
library(RMySQL)

SEQ_FOLDER = "23WES04"


    dbCon = dbConnect(dbDriver("MySQL"), host='192.168.0.35', user='gcx', port=3306, password='gencurix!!', db='gcx_ngs_service' )
    sinfo = dbGetQuery(dbCon, sprintf("SELECT seq_folder,seq_id,sample_name FROM seqid_info WHERE seq_folder='%s'", SEQ_FOLDER))
    fqScreenRes = dbGetQuery(dbCon, "SELECT * FROM gcx_wes.fastqscreen WHERE seq_folder = '%s'", SEQ_FOLDER)
    dbDisconnect(dbCon)
# fqScreenRes

Draw.FastQ.Screen <- function( FqScreenRes, SampleInfo, nCols )
{
    require(ggplot2)
    require(reshape2)
    require(dplyr)
    genomeLevels <- c("human","mouse","rat","ecoli","phix","adapters","vectors")
    COLORS <- c('1-Hit : 1-Genome'    ="#afb781", 'Multi-Hits : 1-Genome'     ="#7b882e", 
                '1-Hit : Multi-Genome'="#edc164", 'Multi-Hits : Multi-Genome' = "#e6a722" )
    fqs <- melt(FqScreenRes) %>% filter(fastq == "R1") %>% 
        filter( seq_id %in% SampleInfo$seq_id) %>%
        mutate( variable = as.character(variable), ID = SampleInfo[match(seq_id, SampleInfo$seq_id), "sample_name"] ) %>% 
        mutate( ID = factor(ID, levels=SampleInfo$sample_name) ) 
    fqs$genome <- sapply(fqs$variable, function(y) unlist(strsplit(y, "_"))[1])
    fqs$type <- sapply(fqs$variable, function(y) unlist(strsplit(y, "_"))[2])
    fqs$type <- gsub("H1","1-Hit : ", fqs$type)
    fqs$type <- gsub("Hm","Multi-Hits : ", fqs$type)
    fqs$type <- gsub("G1","1-Genome", fqs$type)
    fqs$type <- gsub("Gm","Multi-Genome", fqs$type)
    fqs <- fqs %>% filter( type != "unmap" ) %>% mutate( genome = factor(genome, levels=genomeLevels[7:1])) %>%
        mutate( type = factor(type, levels=c("1-Hit : 1-Genome","Multi-Hits : 1-Genome","1-Hit : Multi-Genome","Multi-Hits : Multi-Genome")[4:1]) ) %>% dplyr::rename( AlignType = type )
    
    ggplot(fqs, aes(x=value, y=genome, fill=AlignType)) + geom_bar(position="stack", stat="identity", width=0.85) +
        theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.line.x  = element_line(colour = "#6f6d6d"), 
            axis.line.y  = element_line(colour = "#6f6d6d"),
            legend.position = "bottom",
            legend.key.size = unit(3, 'mm')
        ) +
        scale_fill_manual(values=COLORS) +
        geom_vline(xintercept=0, color="#6f6d6d") +
        facet_wrap(~ID, ncol=nCols)
}

fqscreenplot = Draw.FastQ.Screen( FqScreenRes = fqScreenRes, SampleInfo = sinfo , nCols=4)




