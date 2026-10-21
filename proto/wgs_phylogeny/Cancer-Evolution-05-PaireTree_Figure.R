
###
    options(stringsAsFactors=FALSE) 
    suppressPackageStartupMessages({
        library("Hmisc")
        library("plyr")
        library("dplyr")
        library("setwidth")
        library("RMySQL")
        library("igraph")
        library("ggraph")
        library("ggrepel")
        library("ggplot2")
        library("optparse")
        library("ggplot2")
        library("Cairo")
    })
    grDevices::X11.options(type='cairo')
    options(device='x11')
    options(bitmapType='cairo')
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
    # paritree res
    pairtree_res_list = readRDS(file=sprintf("%s/%s.pairtree.best.tree.node.varaints.data.RDS", WorkDir, BatchID))
    
    best_tree_data      <- pairtree_res_list[["best_tree"]]
    all_vars_of_node    <- pairtree_res_list[["all_nodes_vars"]]
    driver_vars_of_node <- pairtree_res_list[["node_driver_vars"]]
    driver_vars_of_node <- driver_vars_of_node[c(1,3,13),]
    # mutaiontimeR res
    mtr_dir <- "/data/wgs/WGS-26-01/analysis/mutationTimeR"
    res_cnvtime1 <- read.delim(sprintf("%s/%s.mutationTimeR.segment.result.tsv", mtr_dir, tumor_ids[1]))
    res_cnvtime2 <- read.delim(sprintf("%s/%s.mutationTimeR.segment.result.tsv", mtr_dir, tumor_ids[2]))

    # medicc2 res
    medicc2 <- read.delim("/data/wgs/WGS-26-01/analysis/medicc2/WGS-26-01.medicc2.input_branch_lengths.tsv", header=F)

    # clonality 
    clonal_res <- read.delim("/data/wgs/WGS-26-01/analysis/palimpsest/WGS-26-01.palimpsest.clonality.analysis.result.tsv")

    # hmftools driver mutations
    drivermut = rbind(
        read.delim("/data/wgs/WGS-26-01/WGS-26-01_01/hmftools/purple/WGS-26-01_01.purple.driver.catalog.somatic.tsv") %>% filter( driver == "MUTATION" ),
        read.delim("/data/wgs/WGS-26-01/WGS-26-01_02/hmftools/purple/WGS-26-01_02.purple.driver.catalog.somatic.tsv") %>% filter( driver == "MUTATION" )
    )
    drivermut <- drivermut %>% filter( driverLikelihood > 0.5)
    driver_vars_of_node <- driver_vars_of_node %>% filter( Hugo_Symbol %in% drivermut$gene )

    pinfo <- best_tree_data[,c(4,6,9:ncol(best_tree_data))] %>% mutate( cluster_name = paste0("C", node)) %>%
        mutate( 
            phi_WGS.26.01_01 = round(phi_WGS.26.01_01, 3),
            phi_WGS.26.01_02 = round(phi_WGS.26.01_02, 3)

        )
    pinfo %>% filter( phi_WGS.26.01_02 > 0.05, phi_WGS.26.01_01 < 0.05 ) 

###
###
    # evolution tree
    p_df <- best_tree_data %>% mutate( node=as.character(node), parent=as.character(parent))
    # tree edges
    edge_data <- p_df %>% filter(!is.na(parent), parent != "NA") %>% transmute( from = parent, to = node )
    # end-point nodes
    end_point_data <- p_df %>% filter( node %nin% parent ) %>%dplyr::select( !c("variant_ids") )
    colnames(end_point_data) <- gsub("^phi_", "", colnames(end_point_data))
    for( j in 1:nrow(end_point_data) ){
        end_point_sample_id <- names(end_point_data[, 8:ncol(end_point_data)])[ end_point_data[j, 8:ncol(end_point_data)] != 0 ]
        end_sample <- sinfo %>% filter( SampleID %in% end_point_sample_id ) %>% arrange( Order ) %>% .$Group %>% paste0(collapse="\n")
        end_point_data[j, "end_point"] <- end_sample
    }
    # network objects    
    NET <- graph_from_data_frame( edge_data, directed = TRUE )
    # node properties
    net_nodes <- as_data_frame(NET, 'vertices') %>% mutate( node_label = paste0("C", name) )
    net_nodes$driver_mut <- sapply(net_nodes$name, function(v) 
        driver_vars_of_node %>% filter( node_id == as.character(v) ) %>% .$driver_type %>% unique() %>% paste0(collapse=",") 
    )
    net_nodes$driver_mut_lbl <- apply(net_nodes[,c("driver_mut","name")], 1, function(v) {
        if( v[1] == "" ){ 
            lbl = ""
        }else{
            lbl = driver_vars_of_node %>% filter( node_id == as.character(v[2]) ) %>% .$driver_label %>% paste0(collapse="\n")
        }
        lbl = gsub("p.", "", lbl)
        return(lbl)
    })
    meta_node <- c("C8","C10")
    pm_node   <- c("C6","C9")
    net_nodes <- net_nodes %>% mutate(
        node_color = case_when(
            node_label %in% meta_node & driver_mut == "" ~ "#f0dcdc",
            node_label %in% meta_node & driver_mut != "" ~ "#B85450",
            node_label %in% pm_node & driver_mut == "" ~ "#e1e8f2",
            node_label %in% pm_node & driver_mut != "" ~ "#6C8EBF",
            node_label %nin% c(meta_node,pm_node) & driver_mut == "" ~ "#e6efe0",
            node_label %nin% c(meta_node,pm_node) & driver_mut != ""~ "#82b366"
        ),
        node_lbl_color = case_when(
           driver_mut == "" ~ "#2d2d2d",
           driver_mut != "" ~ "#FFFFFF"
        ),
        node_shape  = 21,
        node_border = "#FFFFFF",
        node_stroke = 0.25
    )
    for ( v in colnames(net_nodes)[-1] ){ vertex_attr(NET, v) <- net_nodes[[v]] }
   


    #  network layout
    net_layout <- create_layout(NET, layout = "tree" )

    ggraph(net_layout) + geom_edge_fan() + geom_node_point() + geom_node_text(aes(label=node_label), size=10) 


    new_axis <- net_layout %>% arrange(y) %>% dplyr::select(c("x","y","node_label")) %>% as.data.frame()
    new_axis[c(18:15,13:12,10:8), "x_new"] = seq(0,6,  length.out=9)
    new_axis[c(14),    "x_new"] = 3.2
    new_axis[c(11),    "x_new"] = 4.7
    new_axis[c(6,5,2), "x_new"] = seq(6.5,9,  length.out=4)[-1]
    new_axis[c(7),     "x_new"] = seq(6.5,9,  length.out=4)[2]
    new_axis[c(4),     "x_new"] = seq(6.5,10, length.out=3)[2]
    new_axis[c(1,3),   "x_new"] = 10
    new_axis[18:8,     "y_new"] = 0
    new_axis[c(11,14), "y_new"] = 0.15
    new_axis[c(6,5,2), "y_new"] = 0.3
    new_axis[c(7,4),   "y_new"] = -0.3
    new_axis[c(3),     "y_new"] = -0.15
    new_axis[c(1),     "y_new"] = -0.3
    net_layout$x = new_axis[match(net_layout$node_label, new_axis$node_label), "x_new"]
    net_layout$y = new_axis[match(net_layout$node_label, new_axis$node_label), "y_new"]
 
    pdata_dist = net_layout %>% filter( node_label %in% c("C0","C7","C14","C16"))
    CNV_ENVENT1=c("WGD")
    CNV_ENVENT2=c("11q13-17q12 Translocation\nERBB2 Amplification\nCCND1 Amplification")
    CNV_ENVENT3=c("PTEN Loss")
###
###
    # network plot
    p_evol <- ggraph(net_layout) +
        geom_edge_fan( color = "#4c4c4c", width=0.4, 
            arrow = arrow(type="closed", length = unit(3, "mm")), 
            end_cap = circle(11,"mm"), start_cap = circle(11, 'mm')
        ) +
        coord_cartesian(xlim=c(-0.2, 10.2), ylim=c(-0.55, 0.45) ) +
        geom_node_point( 
            size = 21, pch=net_layout$node_shape, fill=net_layout$node_color, color=net_layout$node_border, stroke=net_layout$node_stroke
        ) +
        geom_node_text(
            aes(label = node_label), fontface="bold", size=8, col=net_layout$node_lbl_color, family="Roboto Condensed"
        ) +
        annotate("text", 
            x = 1.5, y = 0.1, label = "MEN1 D177Mfs*13", color = "#0f0f0f", size=8, vjust=1, hjust=0.85,
            family="Roboto Condensed"
        ) +
        annotate("text", 
            x = 2.25, y = 0.1, label = "KMT2C L3261*", color = "#0f0f0f", size=8, vjust=1, hjust=0.4,
            family="Roboto Condensed"
        ) +
        annotate("text", 
            x = 10, y = -0.375, label = "SUZ12 E584*", color = "#0f0f0f", size=8, vjust=1, hjust=0.5,
            family="Roboto Condensed"
        ) +
        # geom_text_repel(data=net_layout, aes(x=x, y=y, label=driver_mut_lbl), size=6,
        #     box.padding=1, point.padding=1, min.segment.length=0, segment.size=0.2, color="#0f0f0f", force=3
        # ) +
        annotate("text", 
            x = 9, y = 0.3, label = "Primary", color = "#2D7600", size=9, vjust=0.5, hjust=-0.7, fontface="bold",
            family="Roboto Condensed"
        ) +
        annotate("text", 
            x = 10, y = -0.225, label = "Metastasis", color = "#C73500", size=9, vjust=0.5, hjust=-0.2, fontface="bold",
            family="Roboto Condensed"
        ) +
        theme(
            panel.background = element_rect(fill="#FFFFFF"), 
            panel.grid.major = element_blank(),
            plot.margin      = unit(c(1, 1, 1, 1), "cm")
        ) + 
        geom_segment( x=2.5, xend=2.5, y=-0.02, yend=-0.10, color="#006EAF", linewidth=0.5,
            arrow = arrow(length = unit(3, "mm"), type = "closed")
        ) +
        geom_point( x=2.5, y=-0.02, size=2.5, pch=19, color="#006EAF" ) +
        # medicc2 distance
        # geom_segment(x=0, xend=4,  y=0.5,   yend=0.5,  color="#82b366", linewidth=0.3) +
        # geom_segment(x=4, xend=4,  y=0.5,   yend=0.55, color="#82b366", linewidth=0.3) +
        # geom_segment(x=4, xend=6,  y=0.55,  yend=0.55, color="#82b366", linewidth=0.3) +
        # geom_segment(x=4, xend=4,  y=0.45,  yend=0.5,  color="#B85450", linewidth=0.3) +
        # geom_segment(x = 4,    xend = 7.6,  y = 0.45, yend = 0.45, color="#B85450", linewidth=0.3) +
        # geom_segment(x = 7.65, xend = 10,   y = 0.45, yend = 0.45, color="#B85450", linewidth=0.3) +
        # geom_segment(x = 7.55, xend = 7.65, y = 0.43, yend = 0.47, color="#B85450", linewidth=0.4) +
        # geom_segment(x = 7.60, xend = 7.70, y = 0.43, yend = 0.47, color="#B85450", linewidth=0.4) +
        # geom_point( x=0,  y=0.5,  size=6, col="#82b366") + 
        # geom_point( x=4,  y=0.5,  size=6, col="#82b366") + 
        #     geom_point( x=6,  y=0.55, size=6, col="#82b366") + 
        # geom_point( x=10, y=0.45, size=6, col="#B85450") +
        # annotate("text", x = 2, y = 0.5,  label = "101", color = "#2D7600", size=8, vjust=-0.5, hjust=0.5, fontface="bold") +
        # annotate("text", x = 5, y = 0.55, label = "70",  color = "#2D7600", size=8, vjust=-0.5, hjust=0.5, fontface="bold") +
        # annotate("text", x = 8.825, y = 0.45, label = "712", color = "#C73500", size=8, vjust=-0.5, hjust=0.5, fontface="bold") +
        # cnv event
        # geom_segment(x=0, xend=10,  y=-0.5,   yend=-0.5,  color="#a3a3a3", linewidth=0.5,
        #     arrow = arrow(length = unit(5, "mm"), type = "closed")
        # ) +
        # geom_segment(x=0, xend=0, y=-0.05, yend=-0.5, linewidth=0.4, color="#d6d6d6", alpha=0.8) +
        # geom_segment(x=2.75, xend=0.01, y=-0.25, yend=-0.5, linewidth=0.4, color="#d6d6d6", alpha=0.8) +
        # geom_segment(x=2.75, xend=2.75, y=-0.025, yend=-0.25, linewidth=0.4, color="#d6d6d6", alpha=0.8) +
        # geom_segment(x=6, xend=4, y=-0.25, yend=-0.5, linewidth=0.4, color="#d6d6d6", alpha=0.8) +
        # geom_segment(x=6, xend=6, y=-0.05, yend=-0.25, linewidth=0.4, color="#d6d6d6", alpha=0.8) +
        # geom_point(x=0.005, y=-0.5, size=3, pch=19, color="#a3a3a3" ) +
        # geom_point(x=4, y=-0.5, size=3, pch=19, color="#a3a3a3" ) +
        annotate("text", 
            x = 1.5 , y = -0.2, label = CNV_ENVENT1, color = "#006EAF", size=10, vjust=0.5, hjust=0.1, fontface="bold",
            family="Roboto Condensed"
        ) +
        annotate("text", 
            x = 2.5 , y = -0.2, label = CNV_ENVENT2, color = "#006EAF", size=8.5, vjust=0.5, hjust=0.5, fontface="bold",
            family="Roboto Condensed"
        ) +
        annotate("text", 
            x = 3.8 , y = -0.2, label = CNV_ENVENT3, color = "#006EAF", size=9, vjust=0.5, hjust=0.9, fontface="bold",
            family="Roboto Condensed"
        ) +
        geom_segment(x=1.85, xend=1.9, y=-0.2, yend=-0.2, linewidth=0.5, color="#006EAF", 
            arrow = arrow(length = unit(5, "mm"), type = "closed")
        ) +
        geom_segment(x=3.2, xend=3.25, y=-0.2, yend=-0.2, linewidth=0.5, color="#006EAF", 
            arrow = arrow(length = unit(5, "mm"), type = "closed")
        ) #+
        # annotate("text", 
        #     x = -0.5 , y = 0.55, label = "A. CNA-based Phylogenetic Reconstruction (MEDICC2)", 
        #         color = "#0f0f0f", size=7, vjust=0.1, hjust=0, fontface="bold", family="Roboto Condensed"
        # ) +
        # annotate("text", 
        #     x = -0.5 , y = 0.275, label = "B. Inferred Clonal Evolutionary Trajectory (PaireTree)", 
        #         color = "#0f0f0f", size=7, vjust=0.5, hjust=0, fontface="bold", family="Roboto Condensed"
        # ) +
        # annotate("text", 
        #     x = -0.5 , y = -0.35, label = "C. Temporal Ordering of Copy Number Gains \nand Structural Variants (mutationTimeR)", 
        #         color = "#0f0f0f", size=7, vjust=0.5, hjust=0, fontface="bold", family="Roboto Condensed"
        # ) 
###
###
    p_evol

    ggsave(p_evol, file="WGS-26-01.evolution.analysis.summary.raw.png", width=32, height=12, unit='in', dpi=150 )


    net_layout$x[8] = 1.25
    net_layout$y[8] = 2
    net_layout$x[11] = 1.25
    net_layout$y[11] = 1
    net_layout$x[9] = 0.25
    net_layout$y[9] = 2
    net_layout$x[12] = 0.25
    net_layout$y[12] = 1


    net_layout[8, c("x","y")] = c(1.25, 2)
    net_layout[11, c("x","y")] = c(1.25, 1)
    net_layout[9, c("x","y")] = c(0.25, 2)
    net_layout[8, c("x","y")] = c(0.25, 1)

    p_evol <- ggraph(net_layout) +
        geom_edge_fan( color = "#4c4c4c", width=0.4, 
            arrow = arrow(type="closed", length = unit(3, "mm")), 
            end_cap = circle(11,"mm"), start_cap = circle(11, 'mm')
        ) +
        geom_node_point( 
            size = 21, pch=net_layout$node_shape, fill=net_layout$node_color, color=net_layout$node_border, stroke=net_layout$node_stroke
        ) +
        geom_node_text(
            aes(label = node_label), fontface="bold", size=8, col=net_layout$node_lbl_color, family="Roboto Condensed"
        ) +
        geom_text_repel(data=net_layout, aes(x=x, y=y, label=driver_mut_lbl), size=6,
            box.padding=1.1, point.padding=3, min.segment.length=0, segment.size=0, color="#0f0f0f", force=3
        ) +
        coord_flip() + scale_y_reverse() + theme(
            panel.background = element_rect(fill="#FFFFFF"), 
            panel.grid.major = element_blank(),
            plot.margin      = unit(c(1, 1, 1, 1), "cm")
        )
  
    
    ggsave(p_evol, file="WGS-26-01.wgs.tso.evolution.tree.png", width=32, height=12, unit='in', dpi=150 )



