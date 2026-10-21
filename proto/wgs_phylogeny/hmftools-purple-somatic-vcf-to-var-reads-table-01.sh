

#! /bin/bash

### manual params
    BaseDir="/data/wgs/RUO_Service"
    BatchID="WGS-26-02"
    AnalysisList="${BaseDir}/${BatchID}/meta/clonality.evolution.analysis.list"
    OutputDir="${BaseDir}/${BatchID}/analysis/variant-sites-allele-count"
###

### Run : WGS
    # work dir
    mkdir -p "${OutputDir}"

    # genome fasta (hmftools fasta)
    GenomeFasta="/storage/references_and_index/hg19/hmftools/Homo_sapiens.GRCh37.GATK.illumina.fasta"

    # tumor-samples list    
    TumorList=($(cat ${AnalysisList} | grep -v "normal" | awk '{print $1}'))

    # filter and normalization of vcfs : tumors
    for TID in ${TumorList[@]}
    do
        hmftools_purple_vcf=${BaseDir}/${BatchID}/${TID}/hmftools/purple/${TID}.purple.somatic.vcf.gz
        bcftools norm -m -any -f "${GenomeFasta}" "${hmftools_purple_vcf}" \
        | bcftools view -f PASS -s ${TID} \
            -Oz -o ${OutputDir}/${TID}.purple.somatic.filtered.vcf.gz
        tabix -p vcf ${OutputDir}/${TID}.purple.somatic.filtered.vcf.gz
    done

    # merge vcfs and extract AD only
    VCFS=()
    for TID in ${TumorList[@]} 
    do
        VCFS+=("${OutputDir}/${TID}.purple.somatic.filtered.vcf.gz")
    done
    
    # vcf merge and extract variants sites only as table
    bcftools merge -m none -Oz -o "${OutputDir}/${BatchID}.tumors.variants.vcf.gz" "${VCFS[@]}"
    {
        printf "CHROM\tPOS\tREF\tALT\tFILTER"
        for ID in ${TumorList[@]}; do
            printf "\t%s" "$ID"
        done
        printf "\n"

        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%AD]\n' \
            "${OutputDir}/${BatchID}.tumors.variants.vcf.gz"
    } > "${OutputDir}/${BatchID}.tumors.variants.vaf.tsv" 

    # create variants site bed 
    cat ${OutputDir}/${BatchID}.tumors.variants.vaf.tsv | awk -F '\t' '{print $1"\t"$2"\t"$2}' | sed '1d' \
    >  ${OutputDir}/${BatchID}.variants.site.bed

    # count alleles of all variant sites
    RUN_BAMREADCOUNT_COMMADNS=()
    for TID in ${TumorList[@]}
    do
        CMD="singularity exec -B /storage,/data/ /storage/images/bam-readcount.sif bam-readcount -f ${GenomeFasta} -l ${OutputDir}/${BatchID}.variants.site.bed -w 1 -q 20 -b 20 ${BaseDir}/${BatchID}/${TID}/hmftools/redux/${TID}.redux.bam > ${OutputDir}/${TID}.variant.sites.bam.counts.txt"
        RUN_BAMREADCOUNT_COMMADNS+=("$CMD")
    done
    printf '%s\n' "${RUN_BAMREADCOUNT_COMMADNS[@]}" | parallel -j "${#RUN_BAMREADCOUNT_COMMADNS[@]}" -k --line-buffer --halt soon,fail=1

    # summarize bam-readcount results
    for TID in ${TumorList[@]}
    do
        awk 'BEGIN{
        OFS="\t";
        print "chrom","pos","ref","dp","A","T","G","C","ref_count","all_alt_counts"
        }
        {
        chrom=$1
        pos=$2
        ref=$3
        dp=$4

        A=0; T=0; G=0; C=0

        for (i=5; i<=NF; i++) {
            split($i, x, ":")
            base=x[1]
            count=x[2]

            if (base=="A") A=count
            else if (base=="T") T=count
            else if (base=="G") G=count
            else if (base=="C") C=count
        }

        ref_count=0
        if (ref=="A") ref_count=A
        else if (ref=="T") ref_count=T
        else if (ref=="G") ref_count=G
        else if (ref=="C") ref_count=C

        all_alt_counts = A + T + G + C - ref_count

        print chrom,pos,ref,dp,A,T,G,C,ref_count,all_alt_counts
        }' ${OutputDir}/${TID}.variant.sites.bam.counts.txt > ${OutputDir}/${TID}.variant.sites.bam.counts.summary.tsv 
    done
###

### Run : TSO
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
    bcftools merge -m none -Oz -o "${OutputDir}/${BatchID}.tumors.tso.variants.vcf.gz" "${VCFS[@]}"
    {
        printf "CHROM\tPOS\tREF\tALT\tFILTER"
        for ID in ${TumorList[@]}; do
            printf "\t%s" "$ID"
        done
        printf "\n"

        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%AD]\n' \
            "${OutputDir}/${BatchID}.tumors.tso.variants.vcf.gz"
    } > "${OutputDir}/${BatchID}.tumors.tso.variants.vaf.tsv" 

    # create variants site bed 
    cat ${OutputDir}/${BatchID}.tumors.tso.variants.vaf.tsv | awk -F '\t' '{print $1"\t"$2"\t"$2}' | sed '1d' \
    >  ${OutputDir}/${BatchID}.tso.variants.site.bed

    # count alleles of all variant sites
    for TID in ${TumorList[@]}
    do
        bam_file=${BaseDir}/${BatchID}/${TID}/tso500/tso_hmftools/02.redux/${TID}.tso.redux.bam
        singularity exec -B /storage,/data/ /storage/images/bam-readcount.sif bam-readcount \
            -f ${GenomeFasta} \
            -l ${OutputDir}/${BatchID}.tso.variants.site.bed \
            -w 1 -q 20 -b 20 \
            ${bam_file} > ${OutputDir}/${TID}.tso.variant.sites.bam.counts.txt
    done

    # summarize bam-readcount results
    for TID in ${TumorList[@]}
    do
        awk 'BEGIN{
        OFS="\t";
        print "chrom","pos","ref","dp","A","T","G","C","ref_count","all_alt_counts"
        }
        {
        chrom=$1
        pos=$2
        ref=$3
        dp=$4

        A=0; T=0; G=0; C=0

        for (i=5; i<=NF; i++) {
            split($i, x, ":")
            base=x[1]
            count=x[2]

            if (base=="A") A=count
            else if (base=="T") T=count
            else if (base=="G") G=count
            else if (base=="C") C=count
        }

        ref_count=0
        if (ref=="A") ref_count=A
        else if (ref=="T") ref_count=T
        else if (ref=="G") ref_count=G
        else if (ref=="C") ref_count=C

        all_alt_counts = A + T + G + C - ref_count

        print chrom,pos,ref,dp,A,T,G,C,ref_count,all_alt_counts
        }' ${OutputDir}/${TID}.tso.variant.sites.bam.counts.txt > ${OutputDir}/${TID}.tso.variant.sites.bam.counts.summary.tsv 
    done
###
