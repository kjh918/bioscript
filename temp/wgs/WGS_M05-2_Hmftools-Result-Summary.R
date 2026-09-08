###
    options(stringsAsFactors=FALSE) 
    suppressPackageStartupMessages({
        library("Hmisc")
        library("plyr")
        library("dplyr")
        library("setwidth")
        library("RMySQL")
        library("optparse")
        library("ggplot2")
        library("Cairo")
    })
    grDevices::X11.options(type='cairo')
    options(device='x11')
    options(bitmapType='cairo')
###
### Params
    BaseDir        = "/data/wgs/RUO_Service"
    BatchID        = "WGS-26-02"
    GenomeAssembly = "hg19"
    TumorSeqID     = "WGS-26-02_02"
    NormalSeqID    = "WGS-26-02_01"
### 
### Functions and Resources
    Add.Clonality.To.MAF <- function( MAF, MUTATIONTIMER_DATA, PALIMPSEST_DATA ) {
        clonal_convert = c("clonal [late]"="clonal", "subclonal"="subclonal", "clonal [early]"="clonal", "clonal [NA]"="clonal")
        timing_convert = c("clonal [late]"="late", "subclonal"="", "clonal [early]"="early", "clonal [NA]"="")
        RemoveClass <- c("Silent","Intron","5'Flank","3'Flank","Targeted_Region","IGR","3'UTR","5'UTR","RNA")
        MAF <- MAF %>% filter(FILTER %in% c("PASS","common_variant"),Variant_Classification %nin% RemoveClass) %>% 
            filter( t_alt_count >=2 , t_alt_count+t_ref_count >= 10 ) %>% 
            filter( TIER != "LOW_CONFIDENCE" )
        MAF$REF <- apply(MAF[,c("Chromosome","vcf_pos","Variant_Type","Reference_Allele")], 1, function(v) {
            if( v[3] == "INS" ){
                ref_allele = as.character(GENOME[[paste0("chr",v[1])]][as.numeric(v[2]):as.numeric(v[2])])
            }else if( v[3] == "DEL" ){
                ref_allele = paste0(as.character(GENOME[[paste0("chr",v[1])]][as.numeric(v[2]):as.numeric(v[2])]), v[4])
            }else{
                ref_allele = v[4]
            }
            return(ref_allele)
        } )
        MAF$ALT <- apply(MAF[,c("Chromosome","vcf_pos","Variant_Type","Tumor_Seq_Allele2")], 1, function(v) {
            if( v[3] == "INS" ){
                alt_allele = paste0(as.character(GENOME[[paste0("chr",v[1])]][as.numeric(v[2]):as.numeric(v[2])]), v[4])
            }else if( v[3] == "DEL" ){
                alt_allele = as.character(GENOME[[paste0("chr",v[1])]][as.numeric(v[2]):as.numeric(v[2])])
            }else{
                alt_allele = v[4]
            }
            return(alt_allele)
        } )
        MAF_REV <- MAF %>% mutate(index=paste0(Chromosome,":",vcf_pos,":",REF,"_",ALT)) %>% 
            mutate( VAF = Rounding(t_alt_count/(t_alt_count+t_ref_count), digits=3) ) %>%
            dplyr::select(c("index","Hugo_Symbol","Chromosome","vcf_pos","REF","ALT","VAF","HGVSp_Short","Variant_Classification","REPORTED") ) %>% 
            left_join( .,
                MUTATIONTIMER_DATA %>% mutate(
                    index = paste0(CHROM,":",CHROM,":",REF,"_",ALT),
                    clonality_1 = clonal_convert[CLS],
                    time_1 = timing_convert[CLS]
                #) %>% dplyr::select(c( "index","PURPLE_AF","REPORTED","SUBCL","TIER","clonality_1","time_1")),
                ) %>% dplyr::select(c( "index","vaf","clonal_frequency","pSub","TIER","clonality_1","time_1")),
                by="index"
            ) %>% left_join( .,
                PALIMPSEST_DATA %>% mutate( CHROM = gsub("^chr", "", CHROM) ) %>% 
                    mutate( index = paste0(CHROM,":",POS,":",REF,"_",ALT) ) %>% 
                    dplyr::select(c("index","Clonality","Timing","CCF.adj")),
                by="index"
            ) %>% mutate(
                CLONALITY = case_when(
                    clonality_1 == "clonal" & Clonality == "clonal" ~ "clonal",
                    TRUE ~ "subclonal"
                )
            ) %>% dplyr::select(
                gene=Hugo_Symbol, chrom=Chromosome, pos=vcf_pos, hgvsp=HGVSp_Short, variantClass=Variant_Classification,
                purpleReport=REPORTED, purpleTier=TIER, purpleSubclonal=pSub, ccfAdjusted=CCF.adj, gainTiming=time_1, mutationTiming=Timing, 
                ref=REF, alt=ALT, vaf=VAF, purple_af=clonal_frequency, clonality=CLONALITY
            )
        MAF_REV[which(MAF_REV$gainTiming == ""), "gainTiming"] <- NA
        MAF_REV[which(MAF_REV$mutationTiming == ""), "mutationTiming"] <- NA
        return(MAF_REV)
    }
    Rounding <- function(x, digits = 0) {   
        posneg <- sign(x)
        z <- abs(x) * 10^digits
        z <- floor(z + 0.5)
        v <- z / 10^digits * posneg
        return(v)
    }
    RemoveClass <- c("Silent","Intron","5'Flank","3'Flank","Targeted_Region","IGR","3'UTR","5'UTR","RNA")
    HMF_DRIVER_PANALES <- data.table::fread("/storage/references_and_index/hg19/hmftools/resources_links/DRIVER-GENE-PANEL.tsv")
    HMF_DRV_COLUMNS <- c("gene","chromosomeBand","transcript","isCanonical","driver","category","driverLikelihood","biallelic")
    VarClass <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Splice_Region","Splice_Site",
        "Intron","5'UTR","3'UTR","Silent","IGR","5'Flank","3'Flank","RNA","Targeted_Region")
    maf_required_fields <- c(
        'Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', "Strand","Reference_Allele", 
        "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2",
        'Variant_Classification', 'Variant_Type', 'Tumor_Sample_Barcode','Matched_Norm_Sample_Barcode'
    )
    rename_fields <- c(
        'NCBI_Build'='Reference_Genome',
        'FILTER'='FILTER',
        'Existing_variation'='dbSNP',
        't_depth'='Tumor_DP',
        't_ref_count'='Tumor_Ref_Count',
        't_alt_count'='Tumor_Alt_Count',
        'n_depth'='Normal_DP',
        'n_ref_count'='Normal_Ref_Count',
        'n_alt_count'='Normal_Alt_Count',
        'VAF'='VAF',
        'Consequence'='Consequence',
        'SYMBOL'='Gene_Symbol',
        'Transcript_ID'='Ensembl_Transcript_ID',
        'Gene'='Ensembl_Gene_ID',
        'Refseq_Feature'='RefSeq_Transcript_ID',
        'Refseq_Gene'='Entrez_ID',
        'Refseq_HGVSc'='RefSeq_HGVSc',
        'Refseq_HGVSp'='RefSeq_HGVSp',
        'HGVSc'='HGVSc',
        'HGVSp'='HGVSp',
        'HGVSp_Short'='HGVSp_Short',
        'HGVSg'='HGVSg',
        'HGVS_OFFSET'='HGVS_Offset',
        'EXON'='Exon',
        'INTRON'='Intron',
        'cDNA_position'='cDNA_Position',
        'CDS_position'='CDS_Position',
        'Amino_acids'='AminoAcid_Change',
        'Codons'='Codon_Change',
        'Protein_position'='Protein_Position',
        'CANONICAL'='Canonical',
        'AF'='gnomAD_AF',
        'EAS_AF'='gnomAD_EAS_AF',
        'gnomADe_AF'='gnomADe_AF',
        'gnomADe_EAS_AF'='gnomADe_EAS_AF',
        'MAX_AF'='Max_AF',
        'MAX_AF_POPS'='Max_AF_Pops',
        'SIFT'='Sift',
        'PolyPhen'='PolyPhen',
        'GENE_PHENO'='Gene_Pheno',
        'COSMIC'='COSMIC_ID',
        'COSMIC_HGVSC'='COSMIC_HGVSc',
        'COSMIC_LEGACY_ID'='COSMIC_LegacyID',
        'COSMIC_TIER'='COSMIC_TIER',
        'ClinVar_CLNSIG'='Clinvar_Significance',
        'ClinVar_CLNDN'='Clinvar_Disease',
        'ClinVar_ORIGIN'='Clinvar_Origin',
        'am_class'='AlphaMissense_Class',
        'am_pathogenicity'='AlphaMissense_Pathogenicity',
        'SYMBOL_SOURCE'='Symbol_Source',
        'HGNC_ID'='HGNC_ID',
        'BIOTYPE'='BioType',
        'CCDS'='CCDS',
        'SWISSPROT'='SWISSPROT',
        'UNIPARC'='UNIPARC',
        'DISTANCE'='Distance',
        'all_effects'='All_Effects',
        'DOMAINS'='Domains'  , 
        'vcf_pos'='VCF_Position',
        'TIER'='Purple_Tier',
        'HOTSPOT'='Purple_HotSpot',
        'BLACKLIST_VCF'='Purple_Blacklist_VCF',
        'BLACKLIST_BED'='Purple_Blacklist_BED'
    )
    maf_output_fields <- c(maf_required_fields,rename_fields)
    names(maf_output_fields) <- NULL
    if( GenomeAssembly == "hg19" ){
        suppressPackageStartupMessages({library("BSgenome.Hsapiens.UCSC.hg19")})
        GENOME=BSgenome.Hsapiens.UCSC.hg19
    }else{
        suppressPackageStartupMessages({library("BSgenome.Hsapiens.UCSC.hg38")})
        GENOME=BSgenome.Hsapiens.UCSC.hg38
    }
### 
### Folders
    OUTPUT_DIR        <- sprintf("%s/%s/analysis/summary_and_report", BaseDir, BatchID)
    HMF_RES_BASE_DIR  <- sprintf("%s/%s/%s/hmftools", BaseDir, BatchID, TumorSeqID)
    ANALYSIS_BASE_DIR <- sprintf("%s/%s/analysis", BaseDir, BatchID)
###
### Sample Info
    dbCon <- dbConnect(dbDriver("MySQL"), host='192.168.0.34', user='gcx', port=3306, password='gencurix!!', db='gcx_ngs_service' )
    sinfo <- dbGetQuery(dbCon, sprintf("SELECT seq_folder,seq_id,sample_name FROM seqid_info WHERE seq_folder = '%s'", BatchID)) 
    dbDisconnect(dbCon)
###
### READ DATA
    annot_maf   <- data.table::fread(sprintf("%s/purple/%s.purple.somatic.vep.annotation.maf", HMF_RES_BASE_DIR, TumorSeqID), skip=1)
    hmf_drv_res <- data.table::fread(sprintf("%s/purple/%s.purple.driver.catalog.somatic.tsv", HMF_RES_BASE_DIR, TumorSeqID))  
    hmf_cnv_res <- suppressWarnings(data.table::fread(sprintf("%s/purple/%s.purple.cnv.gene.tsv", HMF_RES_BASE_DIR, TumorSeqID)))
    hmf_qc      <- data.table::fread(sprintf("%s/purple/%s.purple.qc", HMF_RES_BASE_DIR, TumorSeqID), header=F) %>% dplyr::rename( category=1, value=2 )
    hmf_purity  <- data.table::fread(sprintf("%s/purple/%s.purple.purity.tsv", HMF_RES_BASE_DIR, TumorSeqID))
    hmf_hrd     <- data.table::fread(sprintf("%s/chord/%s.chord.prediction.tsv", HMF_RES_BASE_DIR, TumorSeqID))
    mutationtime_res <- read.delim(sprintf("%s/clonality-mutationTimeR/%s.mutationTimeR.result.tsv", ANALYSIS_BASE_DIR, TumorSeqID))
    palimpsest_res   <- read.delim(sprintf("%s/clonality-palimpsest/%s.palimpsest.clonality.analysis.result.tsv", ANALYSIS_BASE_DIR, BatchID)) %>% 
        filter( Sample == TumorSeqID )
    annot_maf <- annot_maf %>% mutate( REPORTED = case_when(Hugo_Symbol %in% hmf_drv_res$gene ~ TRUE, TRUE ~ FALSE))
    maf_data  <- Add.Clonality.To.MAF( MAF=annot_maf, MUTATIONTIMER_DATA=mutationtime_res, PALIMPSEST_DATA=palimpsest_res )
###
### SV
    PURPLE_SV_VCF   <- sprintf("%s/purple/%s.purple.sv.vcf.gz", HMF_RES_BASE_DIR, TumorSeqID)
    FILTERED_SV_OUT <- sprintf("%s/purple/%s.purple.filtered.sv.tsv", HMF_RES_BASE_DIR, TumorSeqID)
    OUTPUT_FORAMT   <- "%CHROM\t%POS\t%REF\t%ALT[\t%AD{0}\t%AD{1}\t%AF\t%DP]\t%INFO/SVTYPE\t%INFO/SVID"
    FILTER_CMD      <- sprintf("bcftools view -f PASS -s %s %s | bcftools query -f '%s' > %s", TumorSeqID, PURPLE_SV_VCF, OUTPUT_FORAMT, FILTERED_SV_OUT)
    system(FILTER_CMD)

    sv_res <- data.table::fread(FILTERED_SV_OUT, header=F)

    if( nrow(sv_res) > 0 ){
        SV_RES <- sv_res %>% dplyr::rename( 
            chrom=1, pos=2, ref=3, alt=4, ref_count=5, alt_count=6, af=7, dp=8, sv_type=9, sv_id=10
        ) %>% mutate( 
            BatchID=BatchID, SeqID=TumorSeqID 
        ) %>% mutate( 
            breakpos = paste0(chrom,":",pos)
        ) %>% group_by(sv_id) %>% mutate( 
            breakend = paste(breakpos, collapse="-")
        ) %>% dplyr::select(!c("breakpos"))

        SV_STAT = SV_RES %>% group_by(sv_type) %>% reframe( found=length(unique(sv_id))) %>% mutate( sv_found = sum(found) )
    }else{
        SV_STAT = data.frame( sv_found=0, sv_type="", found=0 )
    }
    
###
### SUMMARY
    # driver mutations
    DRIVERS_MUTATIONS <- rbind(
        maf_data %>% filter( 
            variantClass == "Missense_Mutation",
            gene %in% (HMF_DRIVER_PANALES %>% filter(reportMissense) %>% .$gene)
        ) %>% left_join(., hmf_drv_res %>% dplyr::select(all_of(HMF_DRV_COLUMNS)), by="gene") %>% mutate( driverType = "mutation" ),
        
        maf_data %>% filter( 
            variantClass == "Nonsense_Mutation",
            gene %in% (HMF_DRIVER_PANALES %>% filter(reportNonsense) %>% .$gene)
        ) %>% left_join(., hmf_drv_res %>% dplyr::select(all_of(HMF_DRV_COLUMNS)), by="gene") %>% mutate( driverType = "mutation" ),
    
        maf_data %>% filter( 
            variantClass %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del","In_Frame_Ins"),
            gene %in% (HMF_DRIVER_PANALES %>% filter(reportDisruption) %>% .$gene)
        ) %>% left_join(., hmf_drv_res %>% dplyr::select(all_of(HMF_DRV_COLUMNS)), by="gene") %>% mutate( driverType = "mutation" )
    ) %>% mutate( 
        driverLikelihood = case_when( is.na(driverLikelihood) ~ 0, TRUE ~ driverLikelihood )
    ) %>% mutate(
        HighConf = case_when( driverLikelihood > 0.5 ~ TRUE, TRUE ~ FALSE )
    ) %>% mutate( BatchID=BatchID, SeqID=TumorSeqID ) 

    # driver cnv
    DRIVERS_CNV <- rbind(
        hmf_cnv_res %>% filter( 
            gene %in% ( HMF_DRIVER_PANALES %>% filter(reportAmplification) %>% .$gene),
            maxCopyNumber >= hmf_purity$ploidy * 3
        ) %>% left_join(., HMF_DRIVER_PANALES %>% dplyr::select(c("gene","likelihoodType")), by="gene") %>% mutate(
            driverType="cnv", driver="AMP" 
        ),
        hmf_cnv_res %>% filter( 
            gene %in% ( HMF_DRIVER_PANALES %>% filter(reportDeletion) %>% .$gene),
            minCopyNumber <  0.5,
            maxCopyNumber <= 1
        ) %>% left_join(., HMF_DRIVER_PANALES %>% dplyr::select(c("gene","likelihoodType")), by="gene") %>% mutate(
            driverType="cnv", driver="DEL" 
        )
    )  %>% mutate( BatchID=BatchID, SeqID=TumorSeqID ) %>% dplyr::select(!c("V20","V21"))
    
    # analysis summary
    HMF_SUMMARY <- hmf_purity %>% dplyr::select(
        c("purity","ploidy","gender","wholeGenomeDuplication","msIndelsPerMb","msStatus","tml","tmlStatus","tmbPerMb","tmbStatus","svTumorMutationalBurden")
    ) %>% mutate( 
        Contamination = hmf_qc %>% filter(category=="Contamination") %>% .$value,
        LOH_Rate      = hmf_qc %>% filter(category=="LohPercent") %>% .$value
    ) %>% cbind(., hmf_hrd[,-1]) %>% mutate( BatchID=BatchID, SeqID=TumorSeqID )


    somatic_vars <- annot_maf %>% filter( FILTER == "PASS" ) %>% 
        filter( Variant_Classification %nin% RemoveClass ) %>% 
        filter( t_alt_count >= 2 ) %>% 
        filter( t_alt_count + t_ref_count >= 10 ) %>% 
        filter( TIER != "LOW_CONFIDENCE")
    cnv_amp_genes <- hmf_cnv_res %>% filter( maxCopyNumber >= hmf_purity$ploidy * 3 ) 
    cnv_del_genes <- hmf_cnv_res %>% filter( minCopyNumber < 0.5, maxCopyNumber <= 1 ) 

    VARAINT_STATS <- data.frame(
        raw_variants = nrow(annot_maf),
        filtered_somatic_variants = nrow(somatic_vars),
        driver_mutations = nrow(DRIVERS_MUTATIONS),
        drivder_mutations_highconf = DRIVERS_MUTATIONS %>% filter( HighConf ) %>% nrow(.),
        cnv_amp = cnv_amp_genes  %>% .$gene %>% unique() %>% length(),
        cnv_del = cnv_del_genes  %>% .$gene %>% unique() %>% length(),
        driver_cnv = DRIVERS_CNV %>% .$gene %>% unique() %>% length(),
        sv_found = SV_STAT %>% .$sv_found %>% unique() %>% length()
    ) %>% mutate( BatchID=BatchID, SeqID=TumorSeqID )

    # openxlsx::write.xlsx(
    #     cnv_amp_genes[,c(1,2,3,10,4,5,6,17)] %>% dplyr::rename(chrom=1, start=2, end=3, cytoband=4, gene=5, min_CN=6, max_CN=7, min_minor_allele_CN=8),
    #     sprintf("%s/%s.cnv.amp.genes.result.xlsx", OUTPUT_DIR, TumorSeqID),
    #     rowNames=F, overwrite=T
    # )
    # openxlsx::write.xlsx(
    #     cnv_del_genes[,c(1,2,3,10,4,5,6,17)] %>% dplyr::rename(chrom=1, start=2, end=3, cytoband=4, gene=5, min_CN=6, max_CN=7, min_minor_allele_CN=8),
    #     sprintf("%s/%s.cnv.del.genes.result.xlsx", OUTPUT_DIR, TumorSeqID),
    #     rowNames=F, overwrite=T
    # )
    openxlsx::write.xlsx(
        list(
            AMPLIFICATION=cnv_amp_genes[,c(1,2,3,10,4,5,6,17)] %>% dplyr::rename(chrom=1, start=2, end=3, cytoband=4, gene=5, min_CN=6, max_CN=7, min_minor_allele_CN=8),
            DELETION=cnv_del_genes[,c(1,2,3,10,4,5,6,17)] %>% dplyr::rename(chrom=1, start=2, end=3, cytoband=4, gene=5, min_CN=6, max_CN=7, min_minor_allele_CN=8)
        ),
        sprintf("%s/%s.cnv.genes.result.xlsx", OUTPUT_DIR, (sinfo %>% filter( seq_id == TumorSeqID  ) %>% .$sample_name)),
        rowNames=F, overwrite=T
    )

    somatic_vars_output = somatic_vars %>% mutate( VAF= Rounding(t_alt_count/(t_alt_count+t_ref_count), digits=3) ) %>% 
        dplyr::rename(!!!setNames(names(rename_fields), rename_fields)) %>% 
        dplyr::select(all_of(maf_output_fields)) 
    somatic_vars_output$Tumor_Sample_Barcode        <- sinfo %>% filter( seq_id == TumorSeqID  ) %>% .$sample_name
    somatic_vars_output$Matched_Norm_Sample_Barcode <- sinfo %>% filter( seq_id == NormalSeqID ) %>% .$sample_name

    write.table( somatic_vars_output, 
        sprintf("%s/%s.somatic.variants.filtered.maf", OUTPUT_DIR, (sinfo %>% filter( seq_id == TumorSeqID  ) %>% .$sample_name)),
        quote=F, col.names=T, row.names=F, sep="\t"
    )

    VAR_CLASS_STAT = data.frame( variant_class=VarClass ) %>% 
        left_join(., as.data.frame(table(annot_maf$Variant_Classification)) %>% dplyr::rename(variant_class=1, unfiltered=2), by="variant_class") %>% 
        left_join(., as.data.frame(table(somatic_vars$Variant_Classification)) %>% dplyr::rename(variant_class=1, somatic_filtered=2), by="variant_class") %>% 
        mutate( BatchID=BatchID, SeqID=TumorSeqID ) %>% 
        mutate( variant_class=factor(variant_class, levels=VarClass)) %>% 
        arrange(variant_class)
    VAR_CLASS_STAT[is.na(VAR_CLASS_STAT)] = 0 
### 
### IMPORT INTO DB
    dbCon = dbConnect(dbDriver("MySQL"), host='192.168.0.34', user='gcx', port=3306, password='gencurix!!', db='gcx_wgs' )
    
    prev_clear_stat     <- dbGetQuery(dbCon, sprintf("DELETE FROM hmf_varaint_stats WHERE BatchID = '%s' AND SeqID = '%s'", BatchID, TumorSeqID))
    prev_clear_varclass <- dbGetQuery(dbCon, sprintf("DELETE FROM hmf_var_class_stat WHERE BatchID = '%s' AND SeqID = '%s'", BatchID, TumorSeqID))
    prev_clear_summary  <- dbGetQuery(dbCon, sprintf("DELETE FROM hmf_summary WHERE BatchID = '%s' AND SeqID = '%s'", BatchID, TumorSeqID))
    prev_clear_drvmut   <- dbGetQuery(dbCon, sprintf("DELETE FROM hmf_drivers_mutation WHERE BatchID = '%s' AND SeqID = '%s'", BatchID, TumorSeqID))
    prev_clear_drvcnv   <- dbGetQuery(dbCon, sprintf("DELETE FROM hmf_drivers_cnv WHERE BatchID = '%s' AND SeqID = '%s'", BatchID, TumorSeqID))
    
    update_stat     <- dbWriteTable(dbCon, name="hmf_varaint_stats", value=VARAINT_STATS, row.names=F, append=T)
    update_varclass <- dbWriteTable(dbCon, name="hmf_var_class_stat", value=VAR_CLASS_STAT, row.names=F, append=T)
    update_summary  <- dbWriteTable(dbCon, name="hmf_summary", value=HMF_SUMMARY, row.names=F, append=T)
    update_drvmut   <- dbWriteTable(dbCon, name="hmf_drivers_mutation", value=DRIVERS_MUTATIONS, row.names=F, append=T)
    update_drvcnv   <- dbWriteTable(dbCon, name="hmf_drivers_cnv", value=DRIVERS_CNV, row.names=F, append=T)

    if( nrow(SV_RES) > 0 ){
        prev_clear_svres <- dbGetQuery(dbCon, sprintf("DELETE FROM hmf_sv WHERE BatchID = '%s' AND SeqID = '%s'", BatchID, TumorSeqID))
        update_svres     <- dbWriteTable(dbCon, name="hmf_sv", value=SV_RES, row.names=F, append=T)
    }

    dbDisconnect(dbCon)
###




