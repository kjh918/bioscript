# AML_clonal_prev.csv
# timepoint clone_id clonal_prev

# AML_mutations.csv
# Tier chrom coord clone_id timepoint VAF

# AML_tree_edges.csv
# source target

###

    pairtree_dir <- "/data/wgs/WGS-26-01/analysis/pairtree"
    res_pairtree        <- read.delim(sprintf("%s/%s.pairtree.result.best.tree.node.variants.tsv", pairtree_dir, BatchID))
    all_vars_of_node    <- read.delim(sprintf("%s/%s.pairtree.result.best.tree.node.varriant.integrated.tsv", pairtree_dir, BatchID))
    driver_vars_of_node <- read.delim(sprintf("%s/%s.pairtree.result.best.tree.node.varriant.integrated.drivers.only.tsv", pairtree_dir, BatchID))
    pyclone_res         <- read.delim(sprintf("%s/%s.pyclone-vi.result.tsv", pairtree_dir, BatchID))

    pyclone_res[,c("sample_id","cluster_id","cellular_prevalence")] %>% unique() %>% arrange(sample_id, cluster_id)


    pyres1 = pyclone_res %>% filter( sample_id == "WGS-26-01_01" ) 
    pyres2 = pyclone_res %>% filter( sample_id == "WGS-26-01_02" ) 

    all_vars_of_node$PT_clusterID = pyres1[match(all_vars_of_node$var_index,pyres1$mutation_id ), "cluster_id"]
    all_vars_of_node$BM_clusterID = pyres2[match(all_vars_of_node$var_index,pyres2$mutation_id ), "cluster_id"]

    all_vars_of_node[,c("node_id","PT_clusterID","BM_clusterID")] %>% unique() %>% arrange(PT_clusterID)


clone_id    sample  prevalence
C5  PT    0.3705
C5  BM    0.6691
C3  PT    0.5056
C3  BM    0.0000
C2  PT    0.4943
C2  BM    0.0000
C1  PT    0.3326
C1  BM    0.2986
C4  PT    0.0000
C4  BM    0.3309
C0  PT    0.0000
C0  BM    0.6664

parent  child
C5  C3
C3  C2
C5  C1
C5  C4
C4  C0


    timescape(clonal_prev = clon_prev, tree_edges = clone_tree)

    clonal_res <- read.delim("/data/wgs/WGS-26-01/analysis/palimpsest/WGS-26-01.palimpsest.clonality.analysis.result.tsv")

    pdata_pt <- all_vars_of_node %>% left_join(.,
        clonal_res %>% mutate( var_index=paste(CHROM,POS,REF,ALT, sep="_")) %>% filter(Sample == "WGS-26-01_01" ) %>%  dplyr::select(c("var_index", "CCF.adj")),
        by="var_index"   
    ) %>% filter( !is.na(CCF.adj))
    pdata_bm <- all_vars_of_node %>% left_join(.,
        clonal_res %>% mutate( var_index=paste(CHROM,POS,REF,ALT, sep="_")) %>% filter(Sample == "WGS-26-01_02" ) %>%  dplyr::select(c("var_index", "CCF.adj")),
        by="var_index"   
    )  %>% filter( !is.na(CCF.adj))  

    
    pdata_pt <- pdata_pt %>% group_by(node_id) %>% mutate(clone_prev=median(CCF.adj))
    pdata_bm <- pdata_bm %>% group_by(node_id) %>% mutate(clone_prev=median(CCF.adj))

    pdata_pt[,c("node_id","clone_prev")] %>% unique()
    pdata_bm[,c("node_id","clone_prev")] %>% unique()

WGS-26-01.pyclone-vi.result