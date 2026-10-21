
library(GEOquery)
library(dplyr)

###---| DATASET 1 : GSE248177 |--------------------------------------------------------------------#
dataGSE = getGEO("GSE248177")

selectCols = c("title","geo_accession","characteristics_ch1","characteristics_ch1.3","characteristics_ch1.6","characteristics_ch1.7","characteristics_ch1.8","characteristics_ch1.9","characteristics_ch1.11","characteristics_ch1.12","relation.1","cell line:ch1","culture type:ch1","passage number:ch1","tissue:ch1")

GSE248177 = pData(dataGSE[[1]])[,selectCols]

GSE248177 = GSE248177 %>% dplyr::rename( culture_type='culture type:ch1', SRX='relation.1', cell_type='cell line:ch1',passage='passage number:ch1', tissue='tissue:ch1')

GSE248177$passage = sapply(GSE248177$passage, function(y) ifelse( nchar(y) == 2, gsub("P", "P0", y), y )   )

DataSet1_info = rbind(
    GSE248177 %>% filter( tissue == "RMS1 tissue" ),
    GSE248177 %>% slice(grep("RMS1", GSE248177$cell_type)) %>% filter(culture_type=="3D") %>% arrange(passage),
    GSE248177 %>% filter( tissue == "RMS2 tissue" ),
    GSE248177 %>% slice(grep("RMS2", GSE248177$cell_type)) %>% filter(culture_type=="3D") %>% arrange(passage)
)
openxlsx::write.xlsx(DataSet1_info, "dataset_01_GSE248177_info.xlsx")

#--------------------------------------------------------------------------------------------------#

###---| DATASET 2 : GSE244082 |--------------------------------------------------------------------#
GSE = "GSE102335"
dataGSE = getGEO(GSE)
selectCols = c("title","geo_accession","characteristics_ch1","characteristics_ch1.3","characteristics_ch1.6","characteristics_ch1.7","characteristics_ch1.8","characteristics_ch1.9","characteristics_ch1.11","characteristics_ch1.12","relation.1","cell line:ch1","culture type:ch1","passage number:ch1","tissue:ch1")
dinfo = pData(dataGSE[[1]]) %>% select(any_of(selectCols))
openxlsx::write.xlsx(dinfo, sprintf("%s_info.xlsx", GSE))
#--------------------------------------------------------------------------------------------------#

GSE103826
GSE102335
GSE119004





### 
library(SRAdb)
sqlfile <- getSRAdbFile(destfile = paste(Sys.Date(), "SRAmetadb.sqlite.gz", sep= '.'))
sra_con<- dbConnect(SQLite(), sqlfile)

  sra_dbname <- 'SRAmetadb.sqlite'	
  sra_con <- dbConnect(dbDriver("SQLite"), sra_dbname)
  

##  Get exp accession for run accession:
dbGetQuery(sra_con, "select run_accession, experiment_accession from sra where run_accession = 'SRR1294518'")
  run_accession experiment_accession
1    SRR1294518            SRX549132

## Get more info:
dbGetQuery(sra_con, "select * from sra where run_accession = 'SRR1294518'")

## Query for more than one run accession id:
dbGetQuery(sra_con, "select * from sra where run_accession in ('SRR1294518', 'SRR1294519')")