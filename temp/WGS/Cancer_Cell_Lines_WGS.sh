
#! /bin/bash

###
    RawFastqDir=/storage/DB/cancer_cell_lines_data/wgs_fastq
    TrimFastqDir=/data/wgs/cancer_cell_lines/trimmed_fastq

    for SeqID in SRR1639624_SK.BR.3_WGS SRR1639625_SK.BR.3_C03.WGS SRR1639626_SK.BR.3_C10.WGS SRR1639627_SK.BR.3_C13.WGS SRR1639628_SK.BR.3_C20.WGS SRR1639629_SK.BR.3_C25.WGS SRR1639630_SK.BR.3_C26.WGS SRR8639147_HCC.827_WGS SRR8639152_DU145_WGS SRR8639173_A549_WGS SRR8639189_Capan.1_WGS SRR8652060_KATO.III_WGS SRR8652105_MCF.7_WGS SRR8652114_HT.29_WGS SRR8652122_MIA.PaCa.2_WGS SRR8652123_MKN.1_WGS SRR8670696_NCI.H596_WGS
    do
        singularity exec -B /storage,/data /storage/images/fastp-0.23.4.sif fastp --thread 16 \
            --in1 ${RawFastqDir}/${SeqID}_R1.fastq.gz --in2 ${RawFastqDir}/${SeqID}_R2.fastq.gz \
            --out1 ${TrimFastqDir}/${SeqID}.trimmed_R1.fastq.gz --out2 ${TrimFastqDir}/${SeqID}.trimmed_R2.fastq.gz \
            --json ${TrimFastqDir}/${SeqID}.fastp.json --html ${TrimFastqDir}/${SeqID}.fastp.html \
            --trim_poly_g --detect_adapter_for_pe \
            --length_required 70 --average_qual 10 --qualified_quality_phred 15 
    done

    SEQ_ID_LIST="SRR1639629_SK.BR.3_C25.WGS SRR1639630_SK.BR.3_C26.WGS"
    
    parallel -j 2 -k "singularity exec -B /storage,/data /storage/images/fastp-0.23.4.sif fastp --thread 8 --in1 ${RawFastqDir}/{}_R1.fastq.gz --in2 ${RawFastqDir}/{}_R2.fastq.gz --out1 ${TrimFastqDir}/{}.trimmed_R1.fastq.gz --out2 ${TrimFastqDir}/{}.trimmed_R2.fastq.gz --json ${TrimFastqDir}/{}.fastp.json --html ${TrimFastqDir}/{}.fastp.html --trim_poly_g --detect_adapter_for_pe --length_required 35 --average_qual 10 --qualified_quality_phred 15" ::: $SEQ_ID_LIST
###

###
    qsub -q all.q@ngsmaster /storage/home/kangsm/runScripts/cbNIPT/GCX.cbNIPT_Run.SGE.Sample.Processing.sh


### Manual WGS Processing

    SeqID=cbNIPT_PubDB_05_10_DS2
    RGID=A549WGS20X
    TrimFastqDir=/data/cbNIPT/fastq_trimmed
    BamDir=/data/cbNIPT/processed_dataset/${SeqID}/bam
    mkdir -p /data/cbNIPT/processed_dataset/${SeqID}/bam


    BwaIndex=/storage/references_and_index/hg38/fasta/cbNIPT/hg38.fa
    ReferenceFasta=/storage/references_and_index/hg38/fasta/cbNIPT/hg38.fa
    KnownSnp=/storage/references_and_index/hg38/vcf/Homo_sapiens_assembly38.dbsnp138.vcf.gz
    KnownIndel1=/storage/references_and_index/hg38/vcf/Homo_sapiens_assembly38.known_indels.vcf.gz
    KnownIndel2=/storage/references_and_index/hg38/vcf/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

    singularity exec -B /storage,/data /storage/images/bwa-0.7.17.sif bwa mem -M -t 32 -Y -L 50,50 \
        -R "@RG\tID:${RGID}\tPL:ILLUMINA\tLB:WGS\tSM:${SeqID}\tCN:GCX" ${BwaIndex} \
        ${TrimFastqDir}/${SeqID}.trimmed_R1.fastq.gz ${TrimFastqDir}/${SeqID}.trimmed_R2.fastq.gz | samtools view -bS -o ${BamDir}/${SeqID}.primary.bam -

    samtools sort -@ ${Threads} -o ${BamDir}/${SeqID}.sorted.bam ${BamDir}/${SeqID}.primary.bam
    samtools index -b -@ ${Threads} ${BamDir}/${SeqID}.sorted.bam

    singularity exec -B /storage,/data /storage/images/gatk-4.4.0.0.sif gatk MarkDuplicates --java-options "-XX:ParallelGCThreads=14 -Xmx16384m" \
        --INPUT ${BamDir}/${SeqID}.sorted.bam --OUTPUT ${BamDir}/${SeqID}.sorted.dedup.bam \
        --METRICS_FILE ${BamDir}/${SeqID}.mark.duplicates.metrics.txt --CREATE_INDEX true --REMOVE_SEQUENCING_DUPLICATES true
    
    RecalInputBAM=${BamDir}/${SeqID}.sorted.dedup.bam
    singularity exec -B /storage,/data /storage/images/gatk-4.4.0.0.sif gatk BaseRecalibrator --java-options "-XX:ParallelGCThreads=14 -Xmx16384m" \
        --input ${BamDir}/${SeqID}.sorted.dedup.bam --reference ${ReferenceFasta} --output ${BamDir}/${SeqID}.recal.table.txt \
        --known-sites ${KnownSnp} --known-sites ${KnownIndel1} --known-sites ${KnownIndel2} 
    # run bqsr
    singularity exec -B /storage,/data /storage/images/gatk-4.4.0.0.sif gatk ApplyBQSR --java-options "-XX:ParallelGCThreads=14 -Xmx16384m" \
        --input ${BamDir}/${SeqID}.sorted.dedup.bam --bqsr-recal-file ${BamDir}/${SeqID}.recal.table.txt --output ${BamDir}/${SeqID}.recal.bam 

########




    /storage/home/kangsm/runScripts/cbNIPT/GCX.cbNIPT_run.Analysis.STR.sh \
        --SeqID cbNIPT_PubDB_05_10_DS \
        --BamFile /data/cbNIPT/processed_dataset/cbNIPT_PubDB_05_10_DS/bam/cbNIPT_PubDB_05_10_DS.analysisReady.bam \
        --OutputDir /data/cbNIPT/processed_dataset/cbNIPT_PubDB_05_10_DS/str \
        --SampleID cbnipt.pubdb.05.10_A549.ds.wgs \
        --SampleName A549

    /storage/home/kangsm/runScripts/cbNIPT/GCX.cbNIPT_run.Analysis.STR.sh \
        --SeqID cbNIPT_PubDB_05_10_DS2 \
        --BamFile /data/cbNIPT/processed_dataset/cbNIPT_PubDB_05_10_DS2/bam/cbNIPT_PubDB_05_10_DS2.analysisReady.bam \
        --OutputDir /data/cbNIPT/processed_dataset/cbNIPT_PubDB_05_10_DS2/str \
        --SampleID cbnipt.pubdb.05.10_A549.20x.ds2.wgs \
        --SampleName A549


#------------------------------------------------------------------------------#

    




