

# 1. Convert Purple-Somatic-VCF to Table

#### Params
    BaseDir="/data/wgs/RUO_Service"
    BatchID="WGS-26-02"
    AnalysisList="${BaseDir}/${BatchID}/meta/clonality.evolution.analysis.list"
    OutputDir="${BaseDir}/${BatchID}/analysis/variant-sites-allele-count"
####

#### WGS DATA
    # work dir
    mkdir -p "${OutputDir}"

    # genome fasta
    GenomeFasta="/storage/references_and_index/hg19/hmftools/Homo_sapiens.GRCh37.GATK.illumina.fasta"

    # tumor-samples list    
    TumorList=($(cat ${AnalysisList} | grep -v "normal" | awk '{print $1}'))

    # filter and normalization of vcfs : tumors
    for TID in ${TumorList[@]}
    do
        hmftools_purple_vcf=${BaseDir}/${BatchID}/${TID}/hmftools/purple/${TID}.purple.somatic.vcf.gz

        bcftools norm -m -any -f "${GenomeFasta}" "${hmftools_purple_vcf}" \
        | bcftools view -f PASS -s ${TID} \
            -Oz -o ${OutputDir}/${TID}.wgs.purple.somatic.filtered.vcf.gz
        tabix -p vcf ${OutputDir}/${TID}.wgs.purple.somatic.filtered.vcf.gz
    done

    # merge vcfs and extract AD only
    VCFS=()
    for TID in ${TumorList[@]} 
    do
        VCFS+=("${OutputDir}/${TID}.wgs.purple.somatic.filtered.vcf.gz")
    done
    
    # vcf merge and extract variants sites only as table
    bcftools merge -m none -Oz -o "${OutputDir}/${BatchID}.wgs.tumors.variants.vcf.gz" "${VCFS[@]}"
    {
        printf "CHROM\tPOS\tREF\tALT\tFILTER"
        for ID in ${TumorList[@]}; do
            printf "\t%s" "$ID"
        done
        printf "\n"

        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%AD]\n' \
            "${OutputDir}/${BatchID}.wgs.tumors.variants.vcf.gz"
    } > "${OutputDir}/${BatchID}.wgs.tumors.variants.vaf.tsv" 

    # create variants site bed 
    awk -F '\t' 'BEGIN{OFS="\t"}
    NR>1 {
        chrom=$1
        pos=$2
        ref=$4
        alt=$5

        start = pos - 1
        end   = pos + 1 
        
        print chrom, start, end
    }' ${OutputDir}/${BatchID}.wgs.tumors.variants.vaf.tsv \
    > ${OutputDir}/${BatchID}.wgs.variants.site.bed

    # count alleles of all variant sites
    RUN_BAMREADCOUNT_COMMADNS=()
    for TID in ${TumorList[@]}
    do
        CMD="singularity exec -B /storage,/data/ /storage/images/bam-readcount.sif bam-readcount -f ${GenomeFasta} -l ${OutputDir}/${BatchID}.wgs.variants.site.bed -w 1 -q 20 -b 20 ${BaseDir}/${BatchID}/${TID}/hmftools/redux/${TID}.redux.bam > ${OutputDir}/${TID}.wgs.variant.sites.bam.counts.txt"
        RUN_BAMREADCOUNT_COMMADNS+=("$CMD")
    done
    printf '%s\n' "${RUN_BAMREADCOUNT_COMMADNS[@]}" | parallel -j "${#RUN_BAMREADCOUNT_COMMADNS[@]}" -k --line-buffer --halt soon,fail=1

    # summarize bam-readcount results
    for TID in ${TumorList[@]}
    do
        awk 'BEGIN{
            OFS="\t";
            print "chrom","pos","ref","dp","A","T","G","C","ref_count","snv_alt_counts","indel_counts","all_alt_counts","indel_detail"
        }
        {
            chrom=$1
            pos=$2
            ref=$3
            dp=$4

            A=0; T=0; G=0; C=0
            indel_counts=0
            indel_detail=""

            for (i=5; i<=NF; i++) {
                split($i, x, ":")
                base=x[1]
                count=x[2]

                if (base=="A") A=count
                else if (base=="T") T=count
                else if (base=="G") G=count
                else if (base=="C") C=count
                else if (base ~ /^[+-]/) {
                    indel_counts += count
                    if (indel_detail == "") {
                        indel_detail = base ":" count
                    } else {
                        indel_detail = indel_detail "," base ":" count
                    }
                }
            }

            ref_count=0
            if (ref=="A") ref_count=A
            else if (ref=="T") ref_count=T
            else if (ref=="G") ref_count=G
            else if (ref=="C") ref_count=C

            snv_alt_counts = A + T + G + C - ref_count
            all_alt_counts = snv_alt_counts + indel_counts

            print chrom,pos,ref,dp,A,T,G,C,ref_count,snv_alt_counts,indel_counts,all_alt_counts,indel_detail
        }' ${OutputDir}/${TID}.wgs.variant.sites.bam.counts.txt > ${OutputDir}/${TID}.wgs.variant.sites.bam.counts.summary.tsv 
    done
####

#### TSO500 DATA
    # work dir
    mkdir -p "${OutputDir}"

    # genome fasta
    GenomeFasta="/storage/references_and_index/hg19/hmftools/Homo_sapiens.GRCh37.GATK.illumina.fasta"

    # tumor-samples list    
    TumorList=($(cat ${AnalysisList} | grep -v "normal" | awk '{print $1}'))

    # filter and normalization of vcfs : tumors
    for TID in ${TumorList[@]}
    do
        hmftools_purple_vcf=${BaseDir}/${BatchID}/${TID}/tso500/tso_hmftools/08.purple/${TID}.purple.somatic.vcf.gz

        bcftools norm -m -any -f "${GenomeFasta}" "${hmftools_purple_vcf}" \
        | bcftools view -f PASS -s ${TID} \
            -Oz -o ${OutputDir}/${TID}.tso.purple.somatic.filtered.vcf.gz
        tabix -p vcf ${OutputDir}/${TID}.tso.purple.somatic.filtered.vcf.gz
    done

    # merge vcfs and extract AD only
    VCFS=()
    for TID in ${TumorList[@]} 
    do
        VCFS+=("${OutputDir}/${TID}.tso.purple.somatic.filtered.vcf.gz")
    done
    
    # vcf merge and extract variants sites only as table
    bcftools merge -m none -Oz -o "${OutputDir}/${BatchID}.tso.tumors.variants.vcf.gz" "${VCFS[@]}"
    {
        printf "CHROM\tPOS\tREF\tALT\tFILTER"
        for ID in ${TumorList[@]}; do
            printf "\t%s" "$ID"
        done
        printf "\n"

        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%AD]\n' \
            "${OutputDir}/${BatchID}.tso.tumors.variants.vcf.gz"
    } > "${OutputDir}/${BatchID}.tso.tumors.variants.vaf.tsv" 

    # create variants site bed 
    awk -F '\t' 'BEGIN{OFS="\t"}
    NR>1 {
        chrom=$1
        pos=$2
        ref=$4
        alt=$5

        start = pos - 1
        end   = pos + 1 

        print chrom, start, end
    }' ${OutputDir}/${BatchID}.tso.tumors.variants.vaf.tsv \
    > ${OutputDir}/${BatchID}.tso.variants.site.bed

    # count alleles of all variant sites
    RUN_TSO_BAMREADCOUNT_COMMADNS=()
    for TID in ${TumorList[@]}
    do
        CMD="singularity exec -B /storage,/data/ /storage/images/bam-readcount.sif bam-readcount -f ${GenomeFasta} -l ${OutputDir}/${BatchID}.tso.variants.site.bed -w 1 -q 20 -b 20 ${BaseDir}/${BatchID}/${TID}/tso500/tso_hmftools/02.redux/${TID}.tso.redux.bam > ${OutputDir}/${TID}.tso.variant.sites.bam.counts.txt"
        RUN_TSO_BAMREADCOUNT_COMMADNS+=("$CMD")
    done
    printf '%s\n' "${RUN_TSO_BAMREADCOUNT_COMMADNS[@]}" | parallel -j "${#RUN_TSO_BAMREADCOUNT_COMMADNS[@]}" -k --line-buffer --halt soon,fail=1

    # summarize bam-readcount results
    for TID in ${TumorList[@]}
    do
        awk 'BEGIN{
            OFS="\t";
            print "chrom","pos","ref","dp","A","T","G","C","ref_count","snv_alt_counts","indel_counts","all_alt_counts","indel_detail"
        }
        {
            chrom=$1
            pos=$2
            ref=$3
            dp=$4

            A=0; T=0; G=0; C=0
            indel_counts=0
            indel_detail=""

            for (i=5; i<=NF; i++) {
                split($i, x, ":")
                base=x[1]
                count=x[2]

                if (base=="A") A=count
                else if (base=="T") T=count
                else if (base=="G") G=count
                else if (base=="C") C=count
                else if (base ~ /^[+-]/) {
                    indel_counts += count
                    if (indel_detail == "") {
                        indel_detail = base ":" count
                    } else {
                        indel_detail = indel_detail "," base ":" count
                    }
                }
            }

            ref_count=0
            if (ref=="A") ref_count=A
            else if (ref=="T") ref_count=T
            else if (ref=="G") ref_count=G
            else if (ref=="C") ref_count=C

            snv_alt_counts = A + T + G + C - ref_count
            all_alt_counts = snv_alt_counts + indel_counts

            print chrom,pos,ref,dp,A,T,G,C,ref_count,snv_alt_counts,indel_counts,all_alt_counts,indel_detail
        }' ${OutputDir}/${TID}.tso.variant.sites.bam.counts.txt \
        > ${OutputDir}/${TID}.tso.variant.sites.bam.counts.summary.tsv
    done
```

#### WGS allele count of TSO variant sites

    # tumor-samples list    
    TumorList=($(cat ${AnalysisList} | grep -v "normal" | awk '{print $1}'))

    # count alleles of all variant sites
    RUN_TSOWGS_BAMREADCOUNT_COMMADNS=()
    for TID in ${TumorList[@]}
    do
        CMD="singularity exec -B /storage,/data/ /storage/images/bam-readcount.sif bam-readcount -f ${GenomeFasta} -l ${OutputDir}/${BatchID}.tso.variants.site.bed -w 1 -q 20 -b 20 ${BaseDir}/${BatchID}/${TID}/hmftools/redux/${TID}.redux.bam > ${OutputDir}/${TID}.wgs.variant.tso.sites.bam.counts.txt"
        RUN_TSOWGS_BAMREADCOUNT_COMMADNS+=("$CMD")
    done
    printf '%s\n' "${RUN_TSOWGS_BAMREADCOUNT_COMMADNS[@]}" | parallel -j "${#RUN_TSOWGS_BAMREADCOUNT_COMMADNS[@]}" -k --line-buffer --halt soon,fail=1

    # summarize bam-readcount results
    for TID in ${TumorList[@]}
    do
        awk 'BEGIN{
            OFS="\t";
            print "chrom","pos","ref","dp","A","T","G","C","ref_count","snv_alt_counts","indel_counts","all_alt_counts","indel_detail"
        }
        {
            chrom=$1
            pos=$2
            ref=$3
            dp=$4

            A=0; T=0; G=0; C=0
            indel_counts=0
            indel_detail=""

            for (i=5; i<=NF; i++) {
                split($i, x, ":")
                base=x[1]
                count=x[2]

                if (base=="A") A=count
                else if (base=="T") T=count
                else if (base=="G") G=count
                else if (base=="C") C=count
                else if (base ~ /^[+-]/) {
                    indel_counts += count
                    if (indel_detail == "") {
                        indel_detail = base ":" count
                    } else {
                        indel_detail = indel_detail "," base ":" count
                    }
                }
            }

            ref_count=0
            if (ref=="A") ref_count=A
            else if (ref=="T") ref_count=T
            else if (ref=="G") ref_count=G
            else if (ref=="C") ref_count=C

            snv_alt_counts = A + T + G + C - ref_count
            all_alt_counts = snv_alt_counts + indel_counts

            print chrom,pos,ref,dp,A,T,G,C,ref_count,snv_alt_counts,indel_counts,all_alt_counts,indel_detail
        }' ${OutputDir}/${TID}.wgs.variant.tso.sites.bam.counts.txt > ${OutputDir}/${TID}.wgs.variant.tso.sites.bam.counts.summary.tsv 
    done
####


