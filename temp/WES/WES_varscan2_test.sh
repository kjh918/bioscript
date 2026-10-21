

REF_DIR=/storage/references_and_index/gcx.hg19
GENOME_FASTA=${REF_DIR}/gcx.hg19.fasta
samtools mpileup --min-MQ 1 --fasta-ref ${GENOME_FASTA}\
    --output /data/wes/23WES02/17-100X/bam/17-100X.merged.dup.marked.recal.bam.mpileup \
    /data/wes/23WES02/17-100X/bam/17-100X.merged.dup.marked.recal.bam

java -jar VarScan.jar somatic normal.pileup tumor.pileup output.basename