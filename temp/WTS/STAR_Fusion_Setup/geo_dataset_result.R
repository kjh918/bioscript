    library(plyr)
    library(dplyr)
    library(chimeraviz)
    indexing = function(X) { sapply(unique(X), function(y) list(y)) }
    source("/storage/home/kangsm/myScripts/R_FUNCTIONAL_MODULES/FunMod_StarFusionResults.R")

### | GSE230773 |----------------------------------------------------------------------------------#
    DEMO_DIR = "/storage/home/kangsm/analysis/star-fusion/gse230773_test_run"
    GENOME   = "hg19"

    IDX = indexing(list.files(DEMO_DIR, pattern="rep"))
    IDX = IDX[-4]
    RES_DIR = sprintf("%s/results", DEMO_DIR) 

    sf_res_list <- lapply( IDX, function(id) 
    {
        input_name = sprintf("%s/%s/star-fusion.fusion_predictions.abridged.tsv", DEMO_DIR, id)
        sf = import_starfusion(filename=input_name, genome_version=GENOME)
        print(id)
        return(sf)
    })

    sf_res_summary = lapply( sf_res_list , function(sf) StarFusionList.Summary(sf, COMPACT=T) )

    openxlsx::write.xlsx(sf_res_summary, sprintf("%s/GSE230773.star.fusion.result.summary.xlsx",RES_DIR), rowNames=F, overwrite=T)
#--------------------------------------------------------------------------------------------------#

### | GSE243279 |----------------------------------------------------------------------------------#
    DEMO_DIR = "/storage/home/kangsm/analysis/star-fusion/gse243279_test_run"
    GENOME   = "hg19"
    IDX = indexing(list.files(DEMO_DIR, pattern="rep"))
    RES_DIR = sprintf("%s/results", DEMO_DIR) 

    sf_res_list <- lapply( IDX, function(id) 
    {
        input_name = sprintf("%s/%s/star-fusion.fusion_predictions.abridged.tsv", DEMO_DIR, id)
        sf = import_starfusion(filename=input_name, genome_version=GENOME)
        print(id)
        return(sf)
    })

    sf_res_summary = lapply( sf_res_list , function(sf) StarFusionList.Summary(sf, COMPACT=T) )
    sf_res_summary = lapply( sf_res_summary, function(y) arrange(y, desc(FFPM_mean)))

    openxlsx::write.xlsx(sf_res_summary, sprintf("%s/GSE243279.star.fusion.result.summary.xlsx",RES_DIR), rowNames=F, overwrite=T)
#--------------------------------------------------------------------------------------------------#

