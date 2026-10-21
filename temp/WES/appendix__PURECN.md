
## PureCN INTERVAL FILE GENERATION

#### 1. Prepare Mappability BigWig Files

[ <b>hg38</b> ]  
- wget https://s3.amazonaws.com/purecn/GCA_000001405.15_GRCh38_no_alt_analysis_set_76.bw  
- wget https://s3.amazonaws.com/purecn/GCA_000001405.15_GRCh38_no_alt_analysis_set_100.bw

[ <b>hg19</b> ]
- wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign75mer.bigWig
- wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig


#### 2. Mappability BigWig File Modification
- USCS hg19 mappability bigwig 파일의 chromosome contig가 'chr'을 포함하고 있기 때문에 이부분을 제거하는 과정을 거쳐야함
- BigWig &rarr; Convert into Wig &rarr; Contig Fix &rarr; Convert into BigWig

<pre style="font-family:consolas">
    # BigWig to Wig
    bigWigToWig [<i>input.BigWig</i>] [<i>output.Wig</i>]

    # contig fix
    sed -i 's/^chr//' [<i>Wig</i>]

    # Wig to BigWig
    wigToBigWig.pl [<i>input.wig</i>] [<i>chrom.size</i>] [<i>output.bw</i>]
</pre>

| File | Location |
|:---:|:---|
| chrom.size | /storage/home/kangsm/myDB/wig/hg19/hg19.chrom.sizes |
| hg19 Wig | /storage/home/kangsm/myDB/cnv_purecn/mappability_files/GRCh37_wgEncodeCrgMapabilityAlign100mer.wig |
| hg19 BigWig | /storage/home/kangsm/myDB/cnv_purecn/mappability_files/GRCh37_wgEncodeCrgMapabilityAlign100mer.bw |


#### 2. Prepare WES Capture Library <b>BAIT</b> BED Files

- Target 이 아닌 Bait 영역의 BED 파일을 사용해야 함.

| Genome | Library | Location |
|:---:|:---:|:---|
| hg19 | Twist.Exome.2.0 | /storage/references_and_index/hg19/bed/twist.exome.2.0/hg19_twist.exome.2.0.probe.bed |
| hg19 | SureSelect.V6 | /storage/references_and_index/hg19/bed/sureselect.v6/hg19_sureselect.v6.probe.bed |

#### 3. Run PureCN Interval File Generation Script

<pre style="font-family:consolas">
    PURECN=/storage/apps/R/4.3.1/lib/R/library/PureCN/extdata
    BAIT_BED=[ <i>WES.Library.Bait.bed</i> ]
    GENOME_FASTA=[ <i>Reference.Genome.FASTA</i> ]

    Rscript ${RUN_PURECN}/IntervalFile.R \
        --in-file ${BAIT_BED} \
        --fasta ${GENOME_FASTA} \
        --out-file [ <i>output.intervals.txt</i> ] \
        --off-target \
        --genome [ <i>Genome.Version</i> ] \
        --export [ <i>output.optimized.bed</i> ] \
        --mappability [ <i> Mappability.bw</i> ]
</pre>

<hr>

<pre style="font-family:consolas">
#---| hg19 bigWig Download |---------------------------------------------------#
    wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig
    UCSC_BIGWIG="wgEncodeCrgMapabilityAlign100mer.bigWig"

#---| Chromosome Contig Modification |-----------------------------------------#
    HG19_WIG="GRCh37_wgEncodeCrgMapabilityAlign100mer.wig"
    HG19_BW="GRCh37_wgEncodeCrgMapabilityAlign100mer.bw"
    CHROM_SIZE=/storage/home/kangsm/myDB/wig/hg19/hg19.chrom.sizes

    bigWigToWig $UCSC_BIGWIG $HG19_WIG

    sed -i 's/^chrM/MT/' $HG19_WIG
    sed -i 's/^chr//' $HG19_WIG
    grep -E '1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|22|22|X|Y|MT' $HG19_WIG > $HG19_WIG

    wigToBigWig.pl $HG19_WIG $CHROM_SIZE $HG19_BW
#------------------------------------------------------------------------------#
</pre>

## PureCN Panel-Of-Normals DB ( PON-RDS ) CREATION

| Step | Process |
|:---:|:---|
| 1 | Generate Normal Sample's Coverage |
| 2 | Create Normal-Coverage List |
| 3 | Create PureCN-PON-VCF file |
| 4 | CREATE NORMALDB RDS |

#### 1. Normal Sample's Coverages

<pre style="font-family:consolas">
#---| NORMAL SAMPLE COVERAGES |------------------------------------------------#
    PURECN=/storage/apps/R/4.3.1/lib/R/library/PureCN/extdata

    Rscript $PURECN/Coverage.R \
        --out-dir [ <i>Output.DIR</i> ] \
        --bam [ <i>Input.BAM</i> ] \
        --intervals [ <i>PureCN.Intervals</i> ]
#------------------------------------------------------------------------------#
</pre>

#### 2. Normal Coverage List 
- Normal Sample Coverage File의 위치를 절대경로로 입력 (한줄에 한개씩)
  
#### 3. PureCN-PON-VCF

- A. Create Individual Sample's PON-VCF

<pre style="font-family:consolas">
#
#---| PACKAGES |---------------------------------------------------------------#
    suppressPackageStartupMessages(library("vcfR"))
    suppressPackageStartupMessages(library("dplyr"))
#---| FILES |------------------------------------------------------------------#
    VCF_FILE     <- [ <i>INPUT.VCF</i> ]
    PON_VCF_FILE <- [ <i>OUTPUT.PureCN.PON.VCF</i> ]

#---| READ VCF AND REFORMAT |---------------------------------------------------#
    read_vcf <- read.vcfR(VCF_FILE)
    vcf      <- cbind(as.data.frame(getFIX(read_vcf)), extract.gt( read_vcf, element="AD"))
    vcf_pf   <- vcf[which(vcf$FILTER=="PASS"),]
    #--------------------------------------------------------------------------#
    vcf_pf$INFO   <- "."
    vcf_pf$ID     <- "."
    vcf_pf$QUAL   <- "."
    vcf_pf$FORMAT <- "AD"
    #--------------------------------------------------------------------------#
    vcf_out <- vcf_pf %>% select(c(1,2,3,4,5,6,7,9,10,8)) %>% dplyr::rename("#CHROM"=CHROM) 
#------------------------------------------------------------------------------#

#---| WRITE NEW VCF |----------------------------------------------------------#
    system(sprintf("grep ^##fileformat %s > %s",    VCF_FILE, PON_VCF_FILE))
    system(sprintf("grep ^##FORMAT %s >> %s",       VCF_FILE, PON_VCF_FILE))
    system(sprintf("grep ^##INFO %s >>%s",          VCF_FILE, PON_VCF_FILE))
    system(sprintf("grep ^##contig %s >> %s",       VCF_FILE, PON_VCF_FILE))
    system(sprintf("grep ^##tumor_sample %s >> %s", VCF_FILE, PON_VCF_FILE))
    #--------------------------------------------------------------------------#
    suppressWarnings(write.table(vcf_out, sprintf("%s", PON_VCF_FILE), quote=F, col.names=T, row.names=F, sep="\t", append=T))
#------------------------------------------------------------------------------#
</pre>

- B. Create Individual Sample's PON-VCF

<pre style="font-family:consolas">
#---| Copy & Index VCF |-------------------------------------------------------#    
    cp [ <i>PureCN.PON.vcf</i> ] [ <i>PureCN.PON.VCF.DIR</i> ]
    cd [ <i>PureCN.PON.VCF.DIR</i> ]
    bgzip -c [ <i>PureCN.PON.vcf</i> ] > [ <i>PureCN.PON.vcf.gz</i> ]
    tabix -p vcf [ <i>PureCN.PON.vcf.gz</i> ]
#------------------------------------------------------------------------------#
</pre>

<pre style="font-family:consolas">
#---| VCF MERGE & INDEX |------------------------------------------------------#
    bcftools merge [ <i>PureCN.PON.VCF.DIR</i> ]/*.vcf.gz -o [ <i>PureCN.PON.VCF.DIR</i> ]/[ <i>PureCN.PON.Merged.vcf</i> ]
    bgzip -c [ <i>PureCN.PON.VCF.DIR</i> ]/[ <i>PureCN.PON.Merged.vcf</i> ] > [ <i>PureCN.PON.VCF.DIR</i> ]/[ <i>PureCN.PON.Merged.vcf.gz</i> ]
    tabix -p vcf [ <i>PureCN.PON.VCF.DIR</i> ]/[ <i>PureCN.PON.Merged.vcf.gz</i> ]
#------------------------------------------------------------------------------#
</pre>

#### 4. Create PureCN NormalDB RDS

<pre style="font-family:consolas">
#---| CREATE NORMAL DB RDS |---------------------------------------------------#
    PURECN=/storage/apps/R/4.3.1/lib/R/library/PureCN/extdata

    Rscript $PURECN/NormalDB.R \
        --out-dir [ <i>NORMALDB.DIR</i> ] \
        --coverage-files [ <i>Coverage.List.File</i> ] \
        --normal-panel [ <i>PureCN.PON.Merged.vcf.gz</i> ] \
        --genome [ <i>Reference.Genome.Version</i> ] \
        --assay [ <i>RDS.Name</i> ] \
        --force
#------------------------------------------------------------------------------#
</pre>        


## References Sites/Docs
- https://bioconductor.org/packages/devel/bioc/vignettes/PureCN/inst/doc/Quick.html#41_Coverage

- https://bioconductor.statistik.tu-dortmund.de/packages/3.6/bioc/vignettes/PureCN/inst/doc/Quick.pdf

- https://rstudio-pubs-static.s3.amazonaws.com/444222_0505f4eecb74455ca9fda115989dad5f.html



