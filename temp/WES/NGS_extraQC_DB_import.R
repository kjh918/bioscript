
library(dplyr)
library(RMySQL)
#------------------------------------------------------------------------------#
SeqFolder <- "23WES04"
SeqID     <- "23WES04_30"
extraSpecies <- "mouse"
extraQC_file <- "/storage/home/kangsm/analysis/organoid_science/23WES04/23WES04_30/qcfiles/23WES04_30.mouse.sorted.bam.AlignSummary.metrics.txt"

#------------------------------------------------------------------------------#
StartRow <- system(sprintf("cat %s | grep -n '^## METRICS'", extraQC_file), intern=T)
EndRow   <- system(sprintf("cat %s | grep -n '^## HISTOGRAM'", extraQC_file), intern=T)
StartRow <- as.numeric(unlist(strsplit(StartRow, "\\:"))[1])
EndRow   <- as.numeric(unlist(strsplit(EndRow, "\\:"))[1]) - 2
dataRow  <-EndRow - StartRow -1 
#------------------------------------------------------------------------------#
qcRes = read.delim(extraQC_file, skip=StartRow)[1:dataRow, ] %>% 
    filter( CATEGORY == "PAIR" ) %>% 
    select(c("TOTAL_READS","PF_READS","PF_READS_ALIGNED","PCT_PF_READS","PCT_PF_READS_ALIGNED")) %>% 
    dplyr::rename(rds_total=TOTAL_READS, rds_pf=PF_READS, rds_pf_align=PF_READS_ALIGNED, pct_rds_pf=PCT_PF_READS,pct_rds_pf_align=PCT_PF_READS_ALIGNED ) %>%
    mutate( pct_rds_pf= round(pct_rds_pf*100,3), pct_rds_pf_align= round(pct_rds_pf_align*100,3) )
qcRes$seq_folder <- SeqFolder
qcRes$seq_id     <- SeqID
qcRes$species <- extraSpecies
#------------------------------------------------------------------------------#
dbCon = dbConnect(dbDriver("MySQL"), host='192.168.0.35', user='gcx', port=3306, password='gencurix!!', db='gcx_wes' )
pre_clear <- dbGetQuery(dbCon, sprintf("DELETE FROM extraqc WHERE seq_folder = '%s' AND seq_id = '%s'", SeqFolder, SeqID))
new_write <- dbWriteTable(dbCon, name="extraqc", value=qcRes, row.names=F, append=T)
dbDisconnect(dbCon)
#------------------------------------------------------------------------------#
