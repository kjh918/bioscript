

#### PARAMS
    BaseDir="/data/scRNAseq/02.Public-Datasets"
    BatchID="004.GSE151334.PRJNA635600.cell-line-smartseq2"
    GenomeAssembly="hg38_ercc.plus"
####
    
####
    DataDir="${BaseDir}/${BatchID}/SmartSeq2"
    WorkDir="${BaseDir}/${BatchID}/analysis/SmartSeq2_variant_call"
    mkdir -p ${WorkDir}

    SampleListFile="${BaseDir}/${BatchID}/meta/variant.call.SmartSeq2.sample.list"
    SampleInfoFile="${BaseDir}/${BatchID}/meta/variant.call.SmartSeq2.sample.info"
    SampleList=$(tail -n +2 ${SampleListFile} | awk '{print $1}')
####

####
    GATK_VERSION="4.6.2.0"
    RunSingularity="singularity exec -B /storage,/data /storage/images"
    RunGatk4="${RunSingularity}/gatk-${GATK_VERSION}.sif gatk"
    RunScAllele="${RunSingularity}/scAllele.sif scAllele"
    RunVep="singularity exec -B /storage,/data /storage/images/vep-110.1.sif /opt/vep/src/ensembl-vep/vep"
    RunVcf2Maf="perl /storage/apps/vcf2maf-1.6.22/vcf2maf.pl"

####

#### RESOURCES
    GenomeFasta="/storage/references_and_index/RNAseq/${GenomeAssembly}/fasta/genome.fa"
    GeneExonBed="/storage/references_and_index/RNAseq/${GenomeAssembly}/genes/gene.exon.bed"
####

#### 전처리 부분

    #### READ-GROUP ADD TO BAM  
        
        tail -n +2 "${SampleInfoFile}" | parallel -j 8 --line-buffer --colsep ' ' \
        "
            ${RunGatk4} AddOrReplaceReadGroups \
                -I ${BaseDir}/${BatchID}/SmartSeq2/{2}/bam/{2}.Aligned.sortedByCoord.out.bam \
                -O ${BaseDir}/${BatchID}/SmartSeq2/{2}/bam/{2}.readgroup.added.bam \
                --RGID {2} --RGLB {2} --RGSM {1} --RGPL {4} --RGPU {5}  --RGCN {6} \
                --CREATE_INDEX true
        "

    ####

    #### MERGE BAMS INTO  ONE-BAM-FILE
        MergedBamDir="${WorkDir}/bam"
        mkdir -p ${MergedBamDir}/tmp

        for SampleID in ${SampleList}
        do
            IndivSampleListFile="${BaseDir}/${BatchID}/meta/variant.call.SmartSeq2.sample.${SampleID}.list"
            
            IndivSeqIDList=$(tail -n +2 ${IndivSampleListFile} | awk '{print $2}')
            
            GATK_BAM_MERGE_ARGS=()
            for SeqID in ${IndivSeqIDList}
            do
                GATK_BAM_MERGE_ARGS+=("-I ${BaseDir}/${BatchID}/SmartSeq2/${SeqID}/bam/${SeqID}.readgroup.added.bam")
            done
            GATK_BAM_MERGE_ARGS+=("-O ${MergedBamDir}/${SampleID}.smartseq2.merged.bam")
            GATK_BAM_MERGE_ARGS+=("--CREATE_INDEX true")
            GATK_BAM_MERGE_ARGS+=("--REFERENCE_SEQUENCE ${GenomeFasta}")
            GATK_BAM_MERGE_ARGS+=("--SORT_ORDER coordinate")
            GATK_BAM_MERGE_ARGS+=("--USE_THREADING true")
            GATK_BAM_MERGE_ARGS+=("--VALIDATION_STRINGENCY LENIENT")
            GATK_BAM_MERGE_ARGS+=("--TMP_DIR ${MergedBamDir}/tmp")

            ${RunGatk4} --java-options "-Xmx24g" MergeSamFiles ${GATK_BAM_MERGE_ARGS[@]}
        done
    ####

    #### MARK DUPLICATES

        # for SampleID in ${SampleList}
        # do
        #     ${RunGatk4} --java-options "-Xmx32g" MarkDuplicates \
        #         --INPUT "${MergedBamDir}/${SampleID}.smartseq2.merged.bam" \
        #         --METRICS_FILE "${MergedBamDir}/${SampleID}.mark.duplicates.metrics.txt" \
        #         --OUTPUT "${MergedBamDir}/${SampleID}.smartseq2.merged.dedup.bam" \
        #         --CREATE_INDEX true \
        #         --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
        #         --REMOVE_DUPLICATES false \
        #         --TAGGING_POLICY All \
        #         --MAX_RECORDS_IN_RAM 2000000 
        # done

        parallel -j 4 -k "
            ${RunGatk4} --java-options \"-Xmx32g\" MarkDuplicates \
                --INPUT "${MergedBamDir}/{}.smartseq2.merged.bam" \
                --METRICS_FILE "${MergedBamDir}/{}.mark.duplicates.metrics.txt" \
                --OUTPUT "${MergedBamDir}/{}.smartseq2.merged.dedup.bam" \
                --CREATE_INDEX true \
                --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
                --REMOVE_DUPLICATES false \
                --TAGGING_POLICY All \
                --MAX_RECORDS_IN_RAM 2000000 
        " ::: ${SampleList}


    ####

    #### SplitNCigarReads

        # for SampleID in ${SampleList}
        # do
        #     ${RunGatk4} --java-options "-Xmx16g" SplitNCigarReads \
        #         --reference ${GenomeFasta} \
        #         --input "${MergedBamDir}/${SampleID}.smartseq2.merged.dedup.bam" \
        #         --output "${MergedBamDir}/${SampleID}.smartseq2.merged.dedup.split.bam" \
        #         --create-output-bam-index true 
        # done

        parallel -j 4 -k "
            ${RunGatk4} --java-options "-Xmx32g" SplitNCigarReads \
                --reference ${GenomeFasta} \
                --input "${MergedBamDir}/{}.smartseq2.merged.dedup.bam" \
                --output "${MergedBamDir}/{}.smartseq2.merged.dedup.split.bam" \
                --create-output-bam-index true 
        " ::: ${SampleList}

    ####
####

#### Run HaplotypeCaller (GERMLINE VARIANTS)

    VcfDir="${WorkDir}/vcf"
    mkdir -p ${VcfDir}/tmp

    # run gatk haplotypecaller by individual sample level
    HAPLOTYPECALLER_COMMANDS=()
    for SID in ${SampleList}
    do
        CMD="${RunGatk4} --java-options \"-Xmx16g\" HaplotypeCaller --reference ${GenomeFasta} --input ${MergedBamDir}/${SID}.smartseq2.merged.dedup.split.bam --output ${VcfDir}/${SID}.smartseq2.haplotypecaller.gvcf.gz --intervals ${GeneExonBed} -ERC GVCF --dont-use-soft-clipped-bases true"
        HAPLOTYPECALLER_COMMANDS+=("$CMD")
    done
    printf '%s\n' "${HAPLOTYPECALLER_COMMANDS[@]}" | parallel -j "${#HAPLOTYPECALLER_COMMANDS[@]}" -k --line-buffer --halt soon,fail=1
 
    # sample map file 
    # >> sample <tab> sample's gvcf path
    for gvcf in ${VcfDir}/*.smartseq2.haplotypecaller.gvcf.gz
    do
        sample=$(bcftools query -l "${gvcf}")
        printf "%s\t%s\n" "${sample}" "$(readlink -f "${gvcf}")"
    done > ${VcfDir}/sample_map_gvcfs.txt
    
    # chromosome list (remove ERCC)
    ChromosomeList=$(grep "^ERCC-" -v ${GenomeFasta}.fai | cut -f1 )
    
    # genomicsDB 생성 (chrom 별- 대용량 데이터 대비)
    mkdir -p ${VcfDir}/genomicsDB
    rm -r ${VcfDir}/genomicsDB/*

    GENOMICDBIMPORT_COMMANDS=()
    for CHROM in ${ChromosomeList}
    do
        CMD="${RunGatk4} --java-options \"-Xmx32g\" GenomicsDBImport --reference ${GenomeFasta} --sample-name-map ${VcfDir}/sample_map_gvcfs.txt --genomicsdb-workspace-path ${VcfDir}/genomicsDB/${CHROM} --intervals ${CHROM} --reader-threads 4 --tmp-dir ${VcfDir}/tmp"
        GENOMICDBIMPORT_COMMANDS+=("$CMD")
    done
    printf '%s\n' "${GENOMICDBIMPORT_COMMANDS[@]}" | parallel -j 10 -k --line-buffer --halt soon,fail=1

    # joint genotype call (chrom 별- 대용량 데이터 대비)
    mkdir -p ${VcfDir}/jointGenotype

    GENOTYPEGVCF_COMMMANDS=()
    for CHROM in ${ChromosomeList}
    do
        CMD="${RunGatk4} --java-options \"-Xmx24g\" GenotypeGVCFs -R ${GenomeFasta} -V gendb://${VcfDir}/genomicsDB/${CHROM} -L ${CHROM} -O ${VcfDir}/jointGenotype/smartseq2.haplotypecaller.genotype.${CHROM}.raw.vcf.gz"
        GENOTYPEGVCF_COMMMANDS+=("$CMD")
    done
    printf '%s\n' "${GENOTYPEGVCF_COMMMANDS[@]}" | parallel -j "${#GENOTYPEGVCF_COMMMANDS[@]}" -k --line-buffer --halt soon,fail=1

    # chrom별 genotype vcf 를 하나로 합침 (cohort level로 변환)
    GATHERVCFS_ARGS=()
    for CHROM in ${ChromosomeList}
    do
        ${RunGatk4} ValidateVariants \
            -R ${GenomeFasta} \
            -V ${VcfDir}/jointGenotype/smartseq2.haplotypecaller.genotype.${CHROM}.raw.vcf.gz \
            --validation-type-to-exclude ALLELES

        GATHERVCFS_ARGS+=("-I ${VcfDir}/jointGenotype/smartseq2.haplotypecaller.genotype.${CHROM}.raw.vcf.gz")
    done

    mkdir -p ${VcfDir}/cohortVCF
    ${RunGatk4} --java-options "-Xmx16g" GatherVcfs ${GATHERVCFS_ARGS[@]} \
        -O ${VcfDir}/cohortVCF/smartseq2.haplotypecaller.genotype.cohort.raw.vcf.gz
    tabix -p vcf ${VcfDir}/cohortVCF/smartseq2.haplotypecaller.genotype.cohort.raw.vcf.gz

    # cohort-level-genotype-vcf 정규화
    bcftools norm -f ${GenomeFasta} -m -any -Oz \
        -o ${VcfDir}/cohortVCF/smartseq2.haplotypecaller.genotype.cohort.raw.norm.vcf.gz \
        ${VcfDir}/cohortVCF/smartseq2.haplotypecaller.genotype.cohort.raw.vcf.gz
    tabix -p vcf ${VcfDir}/cohortVCF/smartseq2.haplotypecaller.genotype.cohort.raw.norm.vcf.gz

    # cohort-level-genotype-vcf SNP만 추출
    bcftools view -v snps -m2 -M2 -Oz \
        -o ${VcfDir}/cohortVCF/smartseq2.haplotypecaller.genotype.cohort.snp.raw.vcf.gz \
        ${VcfDir}/cohortVCF/smartseq2.haplotypecaller.genotype.cohort.raw.norm.vcf.gz
    tabix -p vcf ${VcfDir}/cohortVCF/smartseq2.haplotypecaller.genotype.cohort.snp.raw.vcf.gz

    # cohort-level-genotype-vcf 필터 조건을 잡기위한 stats 추출
    ${RunGatk4} --java-options "-Xmx16g" VariantsToTable \
        -V ${VcfDir}/cohortVCF/smartseq2.haplotypecaller.genotype.cohort.snp.raw.vcf.gz \
        -O ${VcfDir}/cohortVCF/smartseq2.haplotypecaller.genotype.cohort.snp.annotations.tsv \
        -F CHROM -F POS -F REF -F ALT -F QUAL -F FILTER \
        -F QD -F FS -F SOR -F MQ -F MQRankSum -F ReadPosRankSum

    # QD : quality by depth : QUAL / DP
    # FS : fisher strand : balanced alt-read strand. 0 = best
    # SOR : strand odds ratio : FS 보완 
    # MQ : mapping quality 
    # MQRankSum : ref-MQ alt-MQ 차이 wilcox test
    # ReadPosRankSum : alt read position

    VCF_STATS_RSCRIPT="/storage/home/kangsm/runScripts/00_Pipeline_Development/Single-Cell-RNAseq/SmartSeq2_HaplotypeCaller_GVCF_Stats.R"
    INPUT_STATS="${VcfDir}/cohortVCF/smartseq2.haplotypecaller.genotype.cohort.snp.annotations.tsv" 
    OUTPUT_STATS="${VcfDir}/cohortVCF/smartseq2.haplotypecaller.genotype.cohort.snp.vcf.stats.summary.tsv" 
    Rscript ${VCF_STATS_RSCRIPT} ${INPUT_STATS} ${OUTPUT_STATS}
    cat ${OUTPUT_STATS} | column -t -s $'\t'

    # stat 결과보고 아래 컷오프 설정
    ###
    Cutoff_QD="2.0"                # default = 2
    Cutoff_FS="60.0"               # default = 60
    Cutoff_SOR="5.0"               # default = 3
    Cutoff_MQ="40.0"               # default = 40
    Cutoff_MQRankSum="-12.5"       # default = -12.5  
    Cutoff_ReadPosRankSum="-8.0"   # default = -8

    # cohort-level-genotype-vcf 'FILTER' 태그 작성
    ${RunGatk4} VariantFiltration \
        --reference ${GenomeFasta} \
        --variant ${VcfDir}/cohortVCF/smartseq2.haplotypecaller.genotype.cohort.snp.raw.vcf.gz \
        --output ${VcfDir}/cohortVCF/smartseq2.haplotypecaller.genotype.cohort.snp.filtered.vcf.gz \
        --filter-name "FFILTER_QD" --filter-expression "QD < ${Cutoff_QD}" \
        --filter-name "FILTER_FS" --filter-expression "FS > ${Cutoff_FS}" \
        --filter-name "FILTER_SOR" --filter-expression "SOR > ${Cutoff_SOR}" \
        --filter-name "FILTER_MQ" --filter-expression "MQ < ${Cutoff_MQ}" \
        --filter-name "FILTER_MQRankSum" --filter-expression "MQRankSum < ${Cutoff_MQRankSum}" \
        --filter-name "FILTER_ReadPosRankSum" --filter-expression "ReadPosRankSum < ${Cutoff_ReadPosRankSum}"

    # cohort-level-genotype-vcf genotype masking
    bcftools +setGT ${VcfDir}/cohortVCF/smartseq2.haplotypecaller.genotype.cohort.snp.filtered.vcf.gz \
        -Oz -o ${VcfDir}/cohortVCF/smartseq2.haplotypecaller.genotype.cohort.snp.GT.masked.vcf.gz \
        -- \
        -t q -n . -i 'FMT/DP<5 || FMT/GQ<20'
    tabix -p vcf ${VcfDir}/cohortVCF/smartseq2.haplotypecaller.genotype.cohort.snp.GT.masked.vcf.gz

    # cohort-level-genotype-vcf low-reliable call masking
    bcftools +setGT ${VcfDir}/cohortVCF/smartseq2.haplotypecaller.genotype.cohort.snp.GT.masked.vcf.gz \
        -Oz -o ${VcfDir}/cohortVCF/smartseq2.haplotypecaller.genotype.cohort.snp.GT.masked2.vcf.gz \
        -- \
        -t q -n . -i 'GT="alt" && FMT/AD[:1]<3'
    tabix -p vcf ${VcfDir}/cohortVCF/smartseq2.haplotypecaller.genotype.cohort.snp.GT.masked2.vcf.gz

    bcftools view -f PASS -Oz -o ${VcfDir}/cohortVCF/smartseq2.haplotypecaller.genotype.cohort.snp.masked.filtered.vcf.gz \
        ${VcfDir}/cohortVCF/smartseq2.haplotypecaller.genotype.cohort.snp.GT.masked2.vcf.gz
    tabix -p vcf ${VcfDir}/cohortVCF/smartseq2.haplotypecaller.genotype.cohort.snp.masked.filtered.vcf.gz

    # 최종 cohort-level-genotype-vcf 복사 (메인 vcf 폴더로)
    COHORT_VCF_FINAL=${VcfDir}/cohortVCF/smartseq2.haplotypecaller.genotype.cohort.snp.masked.filtered.vcf.gz
    COHORT_VCF_COPY=${VcfDir}/Cohort.SmartSeq2.HaplotypeCaller.Cohort.Joint.Genotype.SNP.vcf.gz

    cp ${COHORT_VCF_FINAL} ${COHORT_VCF_COPY}
    tabix -p vcf ${COHORT_VCF_COPY}

    # sample 별 genotype-VCF 생성 
    while read -r sample
    do
        bcftools view -f PASS -s ${sample} -Ou ${COHORT_VCF_COPY} |
        bcftools view -i 'GT="alt"' -Oz -o ${VcfDir}/${sample}.SmartSeq2.HaplotypeCaller.Cohort.Joint.Genotype.SNP.vcf.gz
        tabix -p vcf ${VcfDir}/${sample}.SmartSeq2.HaplotypeCaller.Cohort.Joint.Genotype.SNP.vcf.gz
    done < <(bcftools query -l ${COHORT_VCF_COPY})

    # 공통 SNP 리스트
    PROP=1.0
    CALLABLE=0.7

    # 일정 비율 이상 공통 SNP 
    bcftools view -f PASS -i "F_PASS(GT=\"alt\")>=${PROP}" \
        -Oz -o ${VcfDir}/Cohort.SmartSeq2.HaplotypeCaller.Cohort.Joint.Genotype.Prop.${PROP}.Common.SNP.vcf.gz \
        ${COHORT_VCF_COPY}
    tabix -p vcf ${VcfDir}/Cohort.SmartSeq2.HaplotypeCaller.Cohort.Joint.Genotype.Prop.${PROP}.Common.SNP.vcf.gz
    
    # call rate (=callable) 일정비율 이상이면서 callable 개체의 일정비율 (=prop) 이상인 공통 SNP
    bcftools view \
        -f PASS -i "F_PASS(GT!=\"mis\")>=${CALLABLE} && N_PASS(GT=\"alt\")/N_PASS(GT!=\"mis\")>=${PROP}" \
        -Oz -o ${VcfDir}/Cohort.SmartSeq2.HaplotypeCaller.Cohort.Joint.Genotype.Callable.${CALLABLE}.Prop.${PROP}.Common.SNP.vcf.gz \
        ${COHORT_VCF_COPY}
    tabix -p vcf ${VcfDir}/Cohort.SmartSeq2.HaplotypeCaller.Cohort.Joint.Genotype.Callable.${CALLABLE}.Prop.${PROP}.Common.SNP.vcf.gz

    # MASTER SNP SUMMARY TABLE
    COHORT_SNP_VCF=${VcfDir}/Cohort.SmartSeq2.HaplotypeCaller.Cohort.Joint.Genotype.SNP.vcf.gz
    SUMMARY_OUTPUT=${VcfDir}/Cohort.SmartSeq2.HaplotypeCaller.Cohort.Joint.Genotype.SNP.Summary.tsv
    N_SAMPLE=$(bcftools query -l "${COHORT_SNP_VCF}" | wc -l)
    {
        printf '%s\n' \
            $'CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tQD\tFS\tSOR\tMQ\tReadPosRankSum\tN_SAMPLE\tN_CALLED\tN_MISSING\tN_HOM_REF\tN_HET\tN_HOM_ALT\tN_ALT_CARRIER\tREF_ALLELE_COUNT\tALT_ALLELE_COUNT\tCALL_RATE\tALT_CARRIER_RATE\tAF_CALLED'

        bcftools query \
            -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%QD\t%FS\t%SOR\t%MQ\t%ReadPosRankSum[\t%GT]\n' \
            "${COHORT_SNP_VCF}" |
        awk -v OFS='\t' -v n_sample="${N_SAMPLE}" '
        {
            n_called      = 0
            n_missing     = 0
            n_hom_ref     = 0
            n_het         = 0
            n_hom_alt     = 0
            n_alt_carrier = 0
            ref_ac        = 0
            alt_ac        = 0

            for (i = 13; i <= NF; i++) {
                gt = $i

                # Normalize phased genotypes.
                gsub(/\|/, "/", gt)

                if (gt == "." || gt == "./.") {
                    n_missing++
                }
                else if (gt == "0/0") {
                    n_called++
                    n_hom_ref++
                    ref_ac += 2
                }
                else if (gt == "0/1" || gt == "1/0") {
                    n_called++
                    n_het++
                    n_alt_carrier++
                    ref_ac++
                    alt_ac++
                }
                else if (gt == "1/1") {
                    n_called++
                    n_hom_alt++
                    n_alt_carrier++
                    alt_ac += 2
                }
                else {
                    # Unexpected genotype: treat conservatively as missing.
                    n_missing++
                }
            }
            # Remove sites where no reliable ALT genotype remains
            # after genotype masking.
            if (n_alt_carrier == 0)
                next

            call_rate = (n_sample > 0) \
                ? n_called / n_sample \
                : 0

            carrier_rate = (n_called > 0) \
                ? n_alt_carrier / n_called \
                : 0

            total_ac = ref_ac + alt_ac
            af_called = (total_ac > 0) \
                ? alt_ac / total_ac \
                : 0

            printf "%s", $1
            for (i = 2; i <= 12; i++) {
                printf OFS "%s", $i
            }

            printf OFS "%d", n_sample
            printf OFS "%d", n_called
            printf OFS "%d", n_missing
            printf OFS "%d", n_hom_ref
            printf OFS "%d", n_het
            printf OFS "%d", n_hom_alt
            printf OFS "%d", n_alt_carrier
            printf OFS "%d", ref_ac
            printf OFS "%d", alt_ac
            printf OFS "%.6f", call_rate
            printf OFS "%.6f", carrier_rate
            printf OFS "%.6f\n", af_called
        }'
    } > "${SUMMARY_OUTPUT}"

    # 단일 개체 genoype-vcf
    INDIV_GENOTYPEVCF_COMMANDS=()
    for SID in ${SampleList}
    do
        CMD="${RunGatk4} --java-options \"-Xmx16g\" GenotypeGVCFs -R ${GenomeFasta} -V ${VcfDir}/${SID}.smartseq2.haplotypecaller.gvcf.gz -O ${VcfDir}/${SID}.smartseq2.haplotypecaller.individual.genotype.raw.vcf.gz"
        INDIV_GENOTYPEVCF_COMMANDS+=("$CMD")
    done
    printf '%s\n' "${INDIV_GENOTYPEVCF_COMMANDS[@]}" | parallel -j "${#INDIV_GENOTYPEVCF_COMMANDS[@]}" -k --line-buffer --halt soon,fail=1

    for SID in ${SampleList}
    do
        # cohort-level-genotype-vcf 정규화
        bcftools norm -f ${GenomeFasta} -m -any -Oz \
            -o ${VcfDir}/${SID}.smartseq2.haplotypecaller.individual.genotype.raw.norm.vcf.gz \
            ${VcfDir}/${SID}.smartseq2.haplotypecaller.individual.genotype.raw.vcf.gz
        tabix -p vcf ${VcfDir}/${SID}.smartseq2.haplotypecaller.individual.genotype.raw.norm.vcf.gz

        # cohort-level-genotype-vcf SNP만 추출
        bcftools view -v snps -m2 -M2 -Oz \
            -o ${VcfDir}/${SID}.smartseq2.haplotypecaller.individual.genotype.raw.norm.snp.vcf.gz \
            ${VcfDir}/${SID}.smartseq2.haplotypecaller.individual.genotype.raw.norm.vcf.gz
        tabix -p vcf ${VcfDir}/${SID}.smartseq2.haplotypecaller.individual.genotype.raw.norm.snp.vcf.gz

        ${RunGatk4} VariantFiltration \
            --reference ${GenomeFasta} \
            --variant ${VcfDir}/${SID}.smartseq2.haplotypecaller.individual.genotype.raw.norm.snp.vcf.gz \
            --output ${VcfDir}/${SID}.smartseq2.haplotypecaller.individual.genotype.norm.snp.filtered.vcf.gz \
            --filter-name "FFILTER_QD" --filter-expression "QD < ${Cutoff_QD}" \
            --filter-name "FILTER_FS" --filter-expression "FS > ${Cutoff_FS}" \
            --filter-name "FILTER_SOR" --filter-expression "SOR > ${Cutoff_SOR}" \
            --filter-name "FILTER_MQ" --filter-expression "MQ < ${Cutoff_MQ}" \
            --filter-name "FILTER_MQRankSum" --filter-expression "MQRankSum < ${Cutoff_MQRankSum}" \
            --filter-name "FILTER_ReadPosRankSum" --filter-expression "ReadPosRankSum < ${Cutoff_ReadPosRankSum}"

        # cohort-level-genotype-vcf genotype masking
        bcftools +setGT ${VcfDir}/${SID}.smartseq2.haplotypecaller.individual.genotype.norm.snp.filtered.vcf.gz \
            -Oz -o ${VcfDir}/${SID}.smartseq2.haplotypecaller.individual.genotype.norm.snp.filtered.GT.masked.vcf.gz \
            -- \
            -t q -n . -i 'FMT/DP<5 || FMT/GQ<20'
        tabix -p vcf ${VcfDir}/${SID}.smartseq2.haplotypecaller.individual.genotype.norm.snp.filtered.GT.masked.vcf.gz

        # cohort-level-genotype-vcf low-reliable call masking
        bcftools +setGT ${VcfDir}/${SID}.smartseq2.haplotypecaller.individual.genotype.norm.snp.filtered.GT.masked.vcf.gz \
            -Oz -o ${VcfDir}/${SID}.smartseq2.haplotypecaller.individual.genotype.norm.snp.filtered.GT.masked2.vcf.gz \
            -- \
            -t q -n . -i 'GT="alt" && FMT/AD[:1]<3'
        tabix -p vcf ${VcfDir}/${SID}.smartseq2.haplotypecaller.individual.genotype.norm.snp.filtered.GT.masked2.vcf.gz

        bcftools view -f PASS -i 'FMT/DP>=10' -Oz -o ${VcfDir}/${SID}.SmartSeq2.HaplotypeCaller.Individual.Genotype.SNP.vcf.gz \
            ${VcfDir}/${SID}.smartseq2.haplotypecaller.individual.genotype.norm.snp.filtered.GT.masked2.vcf.gz
        tabix -p vcf ${VcfDir}/${SID}.SmartSeq2.HaplotypeCaller.Individual.Genotype.SNP.vcf.gz

        rm ${VcfDir}/${SID}.smartseq2.haplotypecaller.individual.genotype.raw.*
        rm ${VcfDir}/${SID}.smartseq2.haplotypecaller.individual.genotype.norm.snp.filtered*
    done

####

#### SOMATIC VARIANTS VALIDATION 

    VcfDir="${WorkDir}/vcf"
    mkdir -p ${VcfDir}/tmp

    # run mutect2 tumor-only
    MUTECT2_RUN_COMMANDS=()
    for TASK in $(seq 1 $(tail -n +2 ${SampleListFile} | wc -l))
    do
        RunInfo=$(tail -n +2 ${SampleListFile} | sed -n -e "${TASK} p" - )
        SampleID=$(echo ${RunInfo} | awk '{print $1}')
        GenomeAssembly=$(echo ${RunInfo} | awk '{print $2}')
        Paired=$(echo ${RunInfo} | awk '{print $3}')

        if [[ ${GenomeAssembly} = "hg38_ercc.plus" ]]; then
            MutectGenomeRef="hg38"
        elif [[ ${GenomeAssembly} = "hg19_ercc.plus" ]]; then
            MutectGenomeRef="hg19"
        elif [[ ${GenomeAssembly} = "hg38" ]]; then
            MutectGenomeRef="hg38"
        elif [[ ${GenomeAssembly} = "hg19" ]]; then
            MutectGenomeRef="hg19"
        fi

        VCF_gnomad=/storage/references_and_index/${MutectGenomeRef}/vcf/${MutectGenomeRef}.gnomad.af.only.vcf.gz
        VCF_gnomad_mutect_annot=/storage/references_and_index/${MutectGenomeRef}/vcf/${MutectGenomeRef}.gnomad.split.gnomad-af.vcf.gz
        VCF_pon=/storage/references_and_index/${MutectGenomeRef}/vcf/${MutectGenomeRef}.panel.of.normals.vcf.gz
        VCF_common=/storage/references_and_index/${MutectGenomeRef}/vcf/${MutectGenomeRef}.common.biallelic.snp.vcf.gz 
        VCF_dbsnp=/storage/references_and_index/${MutectGenomeRef}/vcf/${MutectGenomeRef}.dbsnp138.vcf.gz 
        TARGET_INTERVALS=${BaseDir}/${BatchID}/meta/mutation_gene_regions.interval_list

        MUTECT2_ARGS=()
        MUTECT2_ARGS+=("-tumor ${SampleID}")
        MUTECT2_ARGS+=("--input ${MergedBamDir}/${SampleID}.smartseq2.merged.dedup.split.bam")
        MUTECT2_ARGS+=("--output ${VcfDir}/${SampleID}.mutect2.Tonly.unfiltered.vcf")
        if [[ "${Paired}" = "True" ]]; then
            MUTECT2_ARGS+=("--f1r2-tar-gz ${MergedBamDir}/${SampleID}.f1r2.tar.gz")
            MUTECT2_ARGS+=("--f1r2-max-depth 200")
        fi

        CMD="${RunGatk4} --java-options \"-Xmx32g\" Mutect2 ${MUTECT2_ARGS[@]} --max-reads-per-alignment-start 0 -pairHMM AVX_LOGLESS_CACHING_OMP --native-pair-hmm-threads 8 --reference ${GenomeFasta} --intervals ${TARGET_INTERVALS} --germline-resource ${VCF_gnomad} --af-of-alleles-not-in-resource 2.5e-06 --panel-of-normals ${VCF_pon} --max-mnp-distance 0"

        MUTECT2_RUN_COMMANDS+=("$CMD")
   
    done

    printf '%s\n' "${MUTECT2_RUN_COMMANDS[@]}" | parallel -j "${#MUTECT2_RUN_COMMANDS[@]}" -k --line-buffer --halt soon,fail=1

    # run strand-bias modeling (paired-sample only)
    PairedSampleList=$(awk 'NR>1 && $3=="True" {print $1}' ${SampleListFile})
    STRAND_BIAS_MODELING_COMMANDS=()
    for pSeqID in ${PairedSampleList}
    do
        CMD="${RunGatk4} --java-options \"-Xmx8g\" LearnReadOrientationModel -I ${MergedBamDir}/${pSeqID}.f1r2.tar.gz -O ${MergedBamDir}/${pSeqID}.read-orientation-model.tar.gz"
        STRAND_BIAS_MODELING_COMMANDS+=("$CMD")
    done

    printf '%s\n' "${STRAND_BIAS_MODELING_COMMANDS[@]}" | parallel -j "${#STRAND_BIAS_MODELING_COMMANDS[@]}" -k --line-buffer --halt soon,fail=1

    # pileup summary
    PILEUP_SUMMARY_COMMANDS=()
    for TASK in $(seq 1 $(tail -n +2 ${SampleListFile} | wc -l))
    do
        RunInfo=$(tail -n +2 ${SampleListFile} | sed -n -e "${TASK} p" - )
        SampleID=$(echo ${RunInfo} | awk '{print $1}')
        GenomeAssembly=$(echo ${RunInfo} | awk '{print $2}')
        Paired=$(echo ${RunInfo} | awk '{print $3}')

        if [[ ${GenomeAssembly} = "hg38_ercc.plus" ]]; then
            MutectGenomeRef="hg38"
        elif [[ ${GenomeAssembly} = "hg19_ercc.plus" ]]; then
            MutectGenomeRef="hg19"
        elif [[ ${GenomeAssembly} = "hg38" ]]; then
            MutectGenomeRef="hg38"
        elif [[ ${GenomeAssembly} = "hg19" ]]; then
            MutectGenomeRef="hg19"
        fi
        VCF_common=/storage/references_and_index/${MutectGenomeRef}/vcf/${MutectGenomeRef}.common.biallelic.snp.vcf.gz

        CMD="${RunGatk4} --java-options \"-Xmx4g\" GetPileupSummaries -R ${GenomeFasta} -I ${MergedBamDir}/${SampleID}.smartseq2.merged.dedup.split.bam -V ${VCF_common} -L ${TARGET_INTERVALS} -O ${MergedBamDir}/${SampleID}.pileups.table"

        PILEUP_SUMMARY_COMMANDS+=("$CMD")
    done

    printf '%s\n' "${PILEUP_SUMMARY_COMMANDS[@]}" | parallel -j "${#PILEUP_SUMMARY_COMMANDS[@]}" -k --line-buffer --halt soon,fail=1

    # calculate contamination\
    CAL_COMTAM_COMMANDS=()
    for SampleID in ${SampleList}
    do
        CMD="${RunGatk4} CalculateContamination -I ${MergedBamDir}/${SampleID}.pileups.table -O ${MergedBamDir}/${SampleID}.contamination.table --tumor-segmentation ${MergedBamDir}/${SampleID}.segments.table"
        CAL_COMTAM_COMMANDS+=("$CMD")
    done

    printf '%s\n' "${CAL_COMTAM_COMMANDS[@]}" | parallel -j "${#CAL_COMTAM_COMMANDS[@]}" -k --line-buffer --halt soon,fail=1

    # filter mutect2 calls
    FILER_MUTECT2_COMMANDS=()
    for TASK in $(seq 1 $(tail -n +2 ${SampleListFile} | wc -l))
    do
        RunInfo=$(tail -n +2 ${SampleListFile} | sed -n -e "${TASK} p" - )
        SampleID=$(echo ${RunInfo} | awk '{print $1}')
        GenomeAssembly=$(echo ${RunInfo} | awk '{print $2}')
        Paired=$(echo ${RunInfo} | awk '{print $3}')

        if [[ "${Paired}"="True" ]]; then 
            CMD="${RunGatk4} FilterMutectCalls -R ${GenomeFasta} -V ${VcfDir}/${SampleID}.mutect2.Tonly.unfiltered.vcf --contamination-table ${MergedBamDir}/${SampleID}.contamination.table --tumor-segmentation ${MergedBamDir}/${SampleID}.segments.table --ob-priors ${MergedBamDir}/${pSeqID}.read-orientation-model.tar.gz -O ${VcfDir}/${SampleID}.mutect2.Tonly.filtered.vcf"
        else
            CMD="${RunGatk4} FilterMutectCalls -R ${GenomeFasta} -V ${VcfDir}/${SampleID}.mutect2.Tonly.unfiltered.vcf --contamination-table ${MergedBamDir}/${SampleID}.contamination.table --tumor-segmentation ${MergedBamDir}/${SampleID}.segments.table -O ${VcfDir}/${SampleID}.mutect2.Tonly.filtered.vcf"
        fi

        FILER_MUTECT2_COMMANDS+=("$CMD")
    done

    printf '%s\n' "${FILER_MUTECT2_COMMANDS[@]}" | parallel -j "${#FILER_MUTECT2_COMMANDS[@]}" -k --line-buffer --halt soon,fail=1

    # vep annotation

    for SampleID in ${SampleList}
    do
        bcftools norm -m -any -f ${GenomeFasta} -Oz -o ${VcfDir}/${SampleID}.mutect2.Tonly.filtered.norm.vcf.gz ${VcfDir}/${SampleID}.mutect2.Tonly.filtered.vcf
        tabix -p vcf ${VcfDir}/${SampleID}.mutect2.Tonly.filtered.norm.vcf.gz
    done

    VEP_REFSEQ_COMMANDS=()
    for TASK in $(seq 1 $(tail -n +2 ${SampleListFile} | wc -l))
    do
        RunInfo=$(tail -n +2 ${SampleListFile} | sed -n -e "${TASK} p" - )
        SampleID=$(echo ${RunInfo} | awk '{print $1}')
        GenomeAssembly=$(echo ${RunInfo} | awk '{print $2}')
        Paired=$(echo ${RunInfo} | awk '{print $3}')

        if [[ ${GenomeAssembly} = "hg38_ercc.plus" ]]; then
            MutectGenomeRef="hg38"
        elif [[ ${GenomeAssembly} = "hg19_ercc.plus" ]]; then
            MutectGenomeRef="hg19"
        elif [[ ${GenomeAssembly} = "hg38" ]]; then
            MutectGenomeRef="hg38"
        elif [[ ${GenomeAssembly} = "hg19" ]]; then
            MutectGenomeRef="hg19"
        fi

        VepCacheDir=/storage/references_and_index/${MutectGenomeRef}/vep
        VepPluginDir=${VepCacheDir}/Plugins/110
        VepBam=${VepCacheDir}/bam/${MutectGenomeRef}.vep.bam
        VepTopLevelFasta=${VepCacheDir}/toplevel_fasta/${MutectGenomeRef}.toplevel.fa.gz
        if [ "${MutectGenomeRef}" = "hg38" ]; then
            VepGenome="GRCh38"
        else
            VepGenome="GRCh37"
        fi

        CMD="${RunVep} --force_overwrite --offline --input_file ${VcfDir}/${SampleID}.mutect2.Tonly.filtered.norm.vcf.gz --vcf --output_file ${VcfDir}/${SampleID}.mutect2.Tonly.filtered.norm.refseq.vcf --stats_text --show_ref_allele --species homo_sapiens --assembly ${VepGenome} --cache --dir_cache ${VepCacheDir} --dir_plugins ${VepPluginDir} --refseq --use_transcript_ref --hgvs --canonical --pick --mane_select --bam ${VepBam} --fasta ${GenomeFasta} --fork 2 --buffer_size 50000"

        VEP_REFSEQ_COMMANDS+=("$CMD")
    done
    printf '%s\n' "${VEP_REFSEQ_COMMANDS[@]}" | parallel -j "${#VEP_REFSEQ_COMMANDS[@]}" -k --line-buffer --halt soon,fail=1

    VEP_ENSEMBL_COMMANDS=()
    for TASK in $(seq 1 $(tail -n +2 ${SampleListFile} | wc -l))
    do
        RunInfo=$(tail -n +2 ${SampleListFile} | sed -n -e "${TASK} p" - )
        SampleID=$(echo ${RunInfo} | awk '{print $1}')
        GenomeAssembly=$(echo ${RunInfo} | awk '{print $2}')
        Paired=$(echo ${RunInfo} | awk '{print $3}')

        if [[ ${GenomeAssembly} = "hg38_ercc.plus" ]]; then
            MutectGenomeRef="hg38"
        elif [[ ${GenomeAssembly} = "hg19_ercc.plus" ]]; then
            MutectGenomeRef="hg19"
        elif [[ ${GenomeAssembly} = "hg38" ]]; then
            MutectGenomeRef="hg38"
        elif [[ ${GenomeAssembly} = "hg19" ]]; then
            MutectGenomeRef="hg19"
        fi

        VepCacheDir=/storage/references_and_index/${MutectGenomeRef}/vep
        VepPluginDir=${VepCacheDir}/Plugins/110
        VepBam=${VepCacheDir}/bam/${MutectGenomeRef}.vep.bam
        VepTopLevelFasta=${VepCacheDir}/toplevel_fasta/${MutectGenomeRef}.toplevel.fa.gz
        if [ "${MutectGenomeRef}" = "hg38" ]; then
            VepGenome="GRCh38"
        else
            VepGenome="GRCh37"
        fi

        CMD="${RunVep} --force_overwrite --offline --input_file ${VcfDir}/${SampleID}.mutect2.Tonly.filtered.norm.vcf.gz --vcf --output_file ${VcfDir}/${SampleID}.mutect2.Tonly.filtered.norm.ensembl.vcf --stats_text --species homo_sapiens --assembly ${VepGenome} --cache --dir_cache ${VepCacheDir} --dir_plugins ${VepPluginDir} --variant_class --canonical --hgvs --hgvsg --max_af --af_1kg --af_gnomade --symbol --numbers --exclude_predicted --pick --mane_select --bam ${VepBam} --fasta ${GenomeFasta} --fork 2 --buffer_size 50000"

        VEP_ENSEMBL_COMMANDS+=("$CMD")
    done
    printf '%s\n' "${VEP_ENSEMBL_COMMANDS[@]}" | parallel -j "${#VEP_ENSEMBL_COMMANDS[@]}" -k --line-buffer --halt soon,fail=1

    for TASK in $(seq 1 $(tail -n +2 ${SampleListFile} | wc -l))
    do
        RunInfo=$(tail -n +2 ${SampleListFile} | sed -n -e "${TASK} p" - )
        SampleID=$(echo ${RunInfo} | awk '{print $1}')
        GenomeAssembly=$(echo ${RunInfo} | awk '{print $2}')

        if [[ ${GenomeAssembly} = "hg38_ercc.plus" ]]; then
            MutectGenomeRef="hg38"
        elif [[ ${GenomeAssembly} = "hg19_ercc.plus" ]]; then
            MutectGenomeRef="hg19"
        elif [[ ${GenomeAssembly} = "hg38" ]]; then
            MutectGenomeRef="hg38"
        elif [[ ${GenomeAssembly} = "hg19" ]]; then
            MutectGenomeRef="hg19"
        fi
        
        if [ "${MutectGenomeRef}" = "hg38" ]; then
            VepGenome="GRCh38"
        else
            VepGenome="GRCh37"
        fi

        bcftools +split-vep -c Feature,Gene,EXON,INTRON,HGVSc,HGVSp -p Refseq_ -Ou ${VcfDir}/${SampleID}.mutect2.Tonly.filtered.norm.refseq.vcf | bcftools annotate -x INFO/CSQ -Oz -o ${VcfDir}/${SampleID}.refseq.add.vcf.gz
        bcftools index -t ${VcfDir}/${SampleID}.refseq.add.vcf.gz
        bgzip -c ${VcfDir}/${SampleID}.mutect2.Tonly.filtered.norm.ensembl.vcf > ${VcfDir}/${SampleID}.mutect2.Tonly.filtered.norm.ensembl.vcf.gz
        bcftools index -t ${VcfDir}/${SampleID}.mutect2.Tonly.filtered.norm.ensembl.vcf.gz
        bcftools annotate \
            -a ${VcfDir}/${SampleID}.refseq.add.vcf.gz \
            -c INFO/Refseq_Feature,INFO/Refseq_Gene,INFO/Refseq_EXON,INFO/Refseq_INTRON,INFO/Refseq_HGVSc,INFO/Refseq_HGVSp \
            -Ov -o ${VcfDir}/${SampleID}.mutect2.Tonly.filtered.norm.vep.vcf \
            ${VcfDir}/${SampleID}.mutect2.Tonly.filtered.norm.ensembl.vcf.gz
        bgzip -c ${VcfDir}/${SampleID}.mutect2.Tonly.filtered.norm.vep.vcf > ${VcfDir}/${SampleID}.mutect2.Tonly.filtered.norm.vep.vcf.gz
        tabix -p vcf ${VcfDir}/${SampleID}.mutect2.Tonly.filtered.norm.vep.vcf.gz
        bcftools norm -f ${GenomeFasta} -m -any ${VcfDir}/${SampleID}.mutect2.Tonly.filtered.norm.vep.vcf.gz -Ov -o ${VcfDir}/${SampleID}.mutect2.Tonly.filtered.norm.vep.annotated.vcf

        ${RunVcf2Maf} --inhibit-vep \
            --input-vcf ${VcfDir}/${SampleID}.mutect2.Tonly.filtered.norm.vep.annotated.vcf \
            --output-maf ${VcfDir}/${SampleID}.mutect2.Tonly.filtered.norm.vep.annotated.maf \
            --tumor-id ${SampleID} \
            --ref-fasta ${GenomeFasta} \
            --retain-info "Refseq_Feature,Refseq_Gene,Refseq_EXON,Refseq_INTRON,Refseq_HGVSc,Refseq_HGVSp,DP,FILTER" \
            --ncbi-build ${VepGenome} \
            --retain-ann "MAX_AF,MAX_AF_POPS"
        
        rm ${VcfDir}/${SampleID}.refseq.add.vcf.gz ${VcfDir}/${SampleID}.refseq.add.vcf.gz.tbi ${VcfDir}/${SampleID}.mutect2.Tonly.filtered.norm.ensembl.vcf.gz ${VcfDir}/${SampleID}.mutect2.Tonly.filtered.norm.ensembl.vcf.gz.tbi
    done
####



