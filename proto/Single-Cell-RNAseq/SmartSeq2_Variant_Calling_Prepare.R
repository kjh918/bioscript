####
    library(data.table)
    library(plyr)
    library(dplyr)
    library(optparse)
    indexing = function(X) { sapply(unique(X), function(y) list(y)) }
####

####
    option_list = list( 
        make_option(c("--SMARTSEQ2_METADATA"), action="store", default=NULL, type="character", help="[REQUIRED] meta data file (process config)"),
        make_option(c("--MOUSE_ID_COLUMN"), action="store", default=NULL, type="character", help="[REQUIRED] mouse id column name"),   
        make_option(c("--GCX_ID_COLUMN"), action="store", default=NULL, type="character", help="[REQUIRED] gcx-id column name")        
    )
    #--------------------------------------------------------------------------#  
    ARGS = parse_args(OptionParser(option_list=option_list))
    #--------------------------------------------------------------------------#
    SMARTSEQ2_METADATA <- ARGS$SMARTSEQ2_METADATA
    MOUSE_ID_COLUMN    <- ARGS$MOUSE_ID_COLUMN
    GCX_ID_COLUMN      <- ARGS$GCX_ID_COLUMN
####

#### MANUAL PARAMS
    # SMARTSEQ2_METADATA <- "/data/scRNAseq/02.Public-Datasets/002.Tabula-Muris-Senis/meta/metadata_SmartSeq2.tsv"
    # MOUSE_ID_COLUMN    <- "mouse.id"
    # GCX_ID_COLUMN      <- "GCX_ID"
####

### Note
    # SMARTSEQ2_METADATA 
    # column1 : 개체 ID
    # column2 : Seq ID
    # column3 : tissue (or cell type) : annotation info <-- CellMarker annotation 기준이 맞을것으로 보임
###

####
    sinfo       <- fread(SMARTSEQ2_METADATA)
    indiv_sinfo <- data.frame(
        MOUSE_ID = sinfo[[MOUSE_ID_COLUMN]],
        GCX_ID   = sinfo[[GCX_ID_COLUMN]]
    )
    indiv_sinfo_list <- lapply(indexing(indiv_sinfo$MOUSE_ID), function(idx) indiv_sinfo %>% filter( MOUSE_ID == idx ) )
    indiv_list       <- data.frame( indiv_id = unique(indiv_sinfo$MOUSE_ID) )
    output_preform   <- gsub("metadata_SmartSeq2.tsv", "smartseq2_individual_sample_list.INDIV_ID.txt", SMARTSEQ2_METADATA)

    for( j in 1:length(indiv_sinfo_list) )
    {
        output_name <- gsub("INDIV_ID", indiv_list[j, "indiv_id"], output_preform)
        output_list <- indiv_sinfo_list[[ indiv_list[j, "indiv_id"] ]]
        write.table(output_list, output_name, quote=F, col.names=T, row.names=F, sep=" ")
        message(indiv_list[j, "indiv_id"])
    }

    indiv_list_output_name <- gsub("metadata_SmartSeq2.tsv", "smartseq2_individual_samples.txt", SMARTSEQ2_METADATA)
    write.table(indiv_list, indiv_list_output_name, quote=F, col.names=F, row.names=F, sep="\t")
####







