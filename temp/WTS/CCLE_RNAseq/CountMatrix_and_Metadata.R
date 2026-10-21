
library(plyr)
library(dplyr)
library(data.table)

base_dir <- "/data/wts/CELL_LINES"

sinfo <- fread(sprintf("%s/meta/rnaseq_sample_list.txt",base_dir))



