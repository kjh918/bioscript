

#! /bin/bash

    Threads=16
    GenomeIndex=/storage/references_and_index/hg38/fasta/hg38.primary_decoy_HLA_no.ALT.fa

    BaseDir="/data/wgs"
    BatchID="cancer_cell_lines"
    SeqID="SRR8652105_MCF.7"

    FastqR1=${BaseDir}/${BatchID}/${SeqID}/fastq/${SeqID}.trimmed_R1.fastq.gz
    FastqR2=${BaseDir}/${BatchID}/${SeqID}/fastq/${SeqID}.trimmed_R2.fastq.gz

    /storage/apps/circle_finder/circle_finder-pipeline-bwa-mem-samblaster.sh \
        ${Threads} \
        ${GenomeIndex} \
        ${FastqR1} \
        ${FastqR2} \
        10 \
        ${SeqID} \
        "hg38"

