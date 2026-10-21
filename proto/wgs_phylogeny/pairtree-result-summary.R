

###
    options(stringsAsFactors=FALSE) 
    suppressPackageStartupMessages({
        library("Hmisc")
        library("plyr")
        library("dplyr")
        library("setwidth")
        library("optparse")
        library("GenomicRanges")
        library("data.table")
        library("RMySQL")
        library("igraph")
        library("ggraph")
        library("ggrepel")
        library("ggplot2")
        library("Cairo")
        library("BSgenome.Hsapiens.UCSC.hg19")
    })
    grDevices::X11.options(type='cairo')
    options(device='x11')
    options(bitmapType='cairo')
    indexing = function(X) { sapply(unique(X), function(y) list(y)) }
    GENOME=BSgenome.Hsapiens.UCSC.hg19
###
###
    RemoveClasses  <- c("Silent","Intron","5'Flank","3'Flank","Targeted_Region","IGR","3'UTR","5'UTR","RNA")
    CancerType     <- "BRCA"
    HmfToolsDriver <- read.delim("/storage/references_and_index/hg19/hmftools/resources_links/DRIVER-GENE-PANEL.tsv")

    dbCon = dbConnect(dbDriver("MySQL"), host='192.168.0.34', user='gcx', port=3306, password='gencurix!!', db='public_database' )
    IntogenDriver <- dbGetQuery(dbCon, "SELECT * FROM intogen_cancer_driver_catalog")
    dbDisconnect(dbCon)
    IntogenDriver <- IntogenDriver %>% filter( CANCER_TYPE == CancerType ) %>% group_by(SYMBOL) %>% 
        mutate( ACT_MODE=paste0(unique(ROLE), collapse=",") ) %>% dplyr::select(
        c("SYMBOL","CANCER_TYPE","ACT_MODE","IS_DRIVER")
    ) %>% unique() %>% as.data.frame()

    GENOME_FASTA="/storage/references_and_index/hg19/hmftools/Homo_sapiens.GRCh37.GATK.illumina.fasta"
###
###
    Check.Driver.Mutation <- function(GENE, VAR_CLASS, HMF_DRIVER_PANEL) {
        VarClassMatch <- list( 
            'reportMissense' = c("Missense_Mutation"),
            'reportNonsense' = c("Nonsense_Mutation"),
            'reportSplice'   = c("Splice_Region","Splice_Site"),
            'reportDisruption' = c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins")
        )
        matched_class <- which(names(VarClassMatch) %in% names(HMF_DRIVER_PANEL)[HMF_DRIVER_PANEL[which(HMF_DRIVER_PANEL$gene == GENE ), ] == "true"])
        chk_res <- ifelse(
            VAR_CLASS %in% unlist(VarClassMatch[matched_class]),
            TRUE,
            FALSE
        )
        return(chk_res) 
    }
###
###
    BaseDir      <- "/data/wgs/RUO_Service"
    BatchID      <- "WGS-26-02"
    AnalysisList <- sprintf("%s/%s/meta/clonality.evolution.analysis.list", BaseDir, BatchID )
    
    VarDataDir <- sprintf("%s/%s/analysis/variant-sites-allele-count", BaseDir, BatchID)
    WorkDir    <- sprintf("%s/%s/analysis/cancer-evolution-pairtree2", BaseDir, BatchID)
    system(sprintf("mkdir -p %s", WorkDir))

    sinfo <- read.table(AnalysisList, header=F, sep=" ") %>% dplyr::rename(SeqID=1, Group=2, Order=3)
    sinfo$SampleID <- make.names(sinfo$SeqID)
    tumor_ids      <- sinfo %>% filter( Group != "normal" ) %>% .$SeqID
    tumor_samples  <- sinfo %>% filter( Group != "normal" ) %>% .$SampleID
###
###
    # read tree result
    res_tree <- read.delim(sprintf("%s/%s.best.tree.node.variants.tsv", WorkDir, BatchID))
    res_vars <- read.delim(sprintf("%s/%s.pairtree.ssm.txt", WorkDir, BatchID))
    var_data_list <- readRDS(sprintf("%s/%s.wgs.variants.read.count.list.RDS", VarDataDir, BatchID)) 

    VarClassMap <- var_data_list[[1]][,c("vkey","variant_class")] 
    for( i in 2:length(var_data_list) ){
        VarClassMap <- VarClassMap %>% left_join(.,
            var_data_list[[i]][,c("vkey","variant_class")],
            by="vkey"
        )
    } 
    colnames(VarClassMap)[2:ncol(VarClassMap)] <- paste0("v", 1:(ncol(VarClassMap)-1))
    VarClassMap$varclass <- apply(VarClassMap[,-1], 1, function(v) unique(v[!is.na(v)]) )

    var_data <- var_data_list %>% ldply() %>% mutate( 
        name=paste(CHROM,POS,REF,ALT, sep="_"),
        vkey=paste(CHROM,POS, sep="_")
    ) %>% dplyr::select(c("CHROM","POS","REF","ALT","gene","variant_class","hgvsp","name","vkey")) 
    var_data$variant_class = VarClassMap[match(var_data$name, VarClassMap$vkey), "varclass"]
    var_data <- unique(var_data)

    node_to_var_list <- lapply( indexing(as.character(res_tree$node)), function(nid) {
        nid_var_ids <- res_tree %>% filter( node == nid ) %>% .$variant_ids
        nid_var_ids <- unlist(strsplit(nid_var_ids, ","))

        node_vars <- res_vars %>% filter( id %in% nid_var_ids ) %>% 
            left_join(., var_data, by="name") %>%
            #left_join(., var_data, by="vkey") %>% 
            #filter( Variant_Classification %nin% RemoveClasses ) %>% 
            mutate( node_id = nid ) 
        node_vars$HMFTOOLS_DriverGene <- apply(node_vars[,c("gene","variant_class")], 1, function(v) 
            Check.Driver.Mutation(GENE=v[1], VAR_CLASS=v[2], HMF_DRIVER_PANEL=HmfToolsDriver) 
        )
        # node_vars <- node_vars %>% mutate(
        #     INTOGEN_DriverGene = case_when(
        #         gene %in% IntogenDriver$SYMBOL ~ TRUE,
        #         TRUE ~ FALSE
        #     )
        # )   
        return(node_vars)
    })

    all_vars_of_node    <- node_to_var_list %>% ldply() %>% dplyr::select(!c(".id")) 
    coding_vars_of_node <- all_vars_of_node %>% filter( variant_class %nin% RemoveClasses)
    driver_vars_of_node <- all_vars_of_node %>% filter( HMFTOOLS_DriverGene ) %>% 
        mutate( driver_label = paste( gene, hgvsp, sep=" " )) 
    driver_vars_of_node$driver_type <- HmfToolsDriver[match(driver_vars_of_node$gene, HmfToolsDriver$gene), "likelihoodType"]

    write.table( all_vars_of_node, sprintf("%s/%s.pairtree.result.best.tree.node.varriant.integrated.tsv", WorkDir, BatchID), quote=F, col.names=T, row.names=F, sep="\t")
    write.table( coding_vars_of_node, sprintf("%s/%s.pairtree.result.best.tree.node.varriant.integrated.coding.vars.onlytsv", WorkDir, BatchID), quote=F, col.names=T, row.names=F, sep="\t")
    write.table( driver_vars_of_node, sprintf("%s/%s.pairtree.result.best.tree.node.varriant.integrated.drivers.only.tsv", WorkDir, BatchID), quote=F, col.names=T, row.names=F, sep="\t")

    saveRDS(
        list(   
            best_tree = res_tree,
            all_nodes_vars = all_vars_of_node,
            node_coding_vars = coding_vars_of_node,
            node_driver_vars = driver_vars_of_node
        ),
        file=sprintf("%s/%s.pairtree.best.tree.node.varaints.data.RDS", WorkDir, BatchID)
    )


###

