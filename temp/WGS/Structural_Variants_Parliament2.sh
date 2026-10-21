
## Test Sample : cbNIPT_PubDB_04_07 (coriell GM21075) - chr9 inverison RM

    SeqID=cbNIPT_PubDB_04_07

    TrimFastqDir=/data/cbNIPT/fastq_trimmed
    BamDir=/data/cbNIPT/processed_dataset/${SeqID}/bam
    mkdir -p BamDir=/data/cbNIPT/processed_dataset/${SeqID}/bam
    BwaIndex=/storage/references_and_index/hg38/fasta/cbNIPT/hg38.fa
    ReferenceFasta=/storage/references_and_index/hg38/fasta/cbNIPT/hg38.fa
    KnownSnp=/storage/references_and_index/hg38/vcf/Homo_sapiens_assembly38.dbsnp138.vcf.gz
    KnownIndel1=/storage/references_and_index/hg38/vcf/Homo_sapiens_assembly38.known_indels.vcf.gz
    KnownIndel2=/storage/references_and_index/hg38/vcf/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

    singularity exec -B /storage,/data /storage/images/bwa-0.7.17.sif bwa mem -M -t 32 -Y -L 50,50 \
        -R "@RG\tID:GM21075\tPL:ILLUMINA\tLB:WGS\tSM:${SeqID}\tCN:GCX" ${BwaIndex} \
        ${TrimFastqDir}/${SeqID}.trimmed_R1.fastq.gz ${TrimFastqDir}/${SeqID}.trimmed_R2.fastq.gz | samtools view -bS -o ${BamDir}/${SeqID}.primary.bam -

    samtools sort -@ ${Threads} -o ${BamDir}/${SeqID}.sorted.bam ${BamDir}/${SeqID}.primary.bam
    samtools index -b -@ ${Threads} ${BamDir}/${SeqID}.sorted.bam

    singularity exec -B /storage,/data /storage/images/gatk-4.4.0.0.sif gatk MarkDuplicates --java-options "-XX:ParallelGCThreads=14 -Xmx16384m" \
        --INPUT ${BamDir}/${SeqID}.sorted.bam --OUTPUT ${BamDir}/${SeqID}.sorted.dedup.bam \
        --METRICS_FILE ${BamDir}/${SeqID}.mark.duplicates.metrics.txt --CREATE_INDEX true --REMOVE_SEQUENCING_DUPLICATES true
    
    RecalInputBAM=${BamDir}/${SeqID}.sorted.dedup.bam
    singularity exec -B /storage,/data /storage/images/gatk-4.4.0.0.sif gatk BaseRecalibrator --java-options "-XX:ParallelGCThreads=14 -Xmx16384m" \
        --input ${RecalInputBAM} --reference ${ReferenceFasta} --output ${BamDir}/${SeqID}.recal.table.txt \
        --known-sites ${KnownSnp} --known-sites ${KnownIndel1} --known-sites ${KnownIndel2} 
    # run bqsr
    singularity exec -B /storage,/data /storage/images/gatk-4.4.0.0.sif gatk ApplyBQSR --java-options "-XX:ParallelGCThreads=14 -Xmx16384m" \
        --input ${RecalInputBAM} --bqsr-recal-file ${BamDir}/${SeqID}.recal.table.txt --output ${BamDir}/${SeqID}.recal.bam 

###

### 
    # cbNIPT_PubDB_04_07 cbNIPT_PubDB_04_07_DS

    BaseDir=/storage/home/kangsm/analysis/gcx.cbNIPT/04.ANALYSIS_DEV/Structural_variations

   
    InputDir=/storage/home/kangsm/analysis/gcx.cbNIPT/04.ANALYSIS_DEV/Structural_variations/wgs_bam
    OutputDir=/storage/home/kangsm/analysis/gcx.cbNIPT/04.ANALYSIS_DEV/Structural_variations/results/wgs_bam
    InputPrefix=GM21075.WGS
    OutputPrefix=GM21075.WGS
    
    singularity exec -B ${InputDir}:/home/dnanexus/in,${OutputDir}:/home/dnanexus/out -B /storage,/data /storage/images/parliament2-0.1.11.sif parliament2.py \
        --bam ${InputPrefix}.bam --bai ${InputPrefix}.bam.bai --ref_genome ref.fa --fai ref.fa.fai --prefix ${OutputPrefix} \
        --breakdancer --delly_insertion --delly_inversion 
