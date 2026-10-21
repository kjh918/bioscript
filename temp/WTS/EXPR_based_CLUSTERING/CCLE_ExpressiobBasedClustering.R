

#---| PACKAGES |-----------------------------------------------------------------------------------#
    options(stringsAsFactors=FALSE)   
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("umap"))
    suppressPackageStartupMessages(library("openxlsx"))
    suppressPackageStartupMessages(library("Hmisc"))
    suppressPackageStartupMessages(library("scales"))
    suppressPackageStartupMessages(library("plyr"))
    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("NMF"))
    suppressPackageStartupMessages(library("ConsensusClusterPlus"))
    suppressPackageStartupMessages(library("parallel"))
#--------------------------------------------------------------------------------------------------#

#---| FOLDERS |------------------------------------------------------------------------------------#
    ccle_data_dir = "/storage/home/kangsm/myDB/ccle_rnaseq"
    res_dir       = "/storage/home/kangsm/analysis/gcx.rnaseq.based.clustering"
#--------------------------------------------------------------------------------------------------#



