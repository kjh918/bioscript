
BamFile1=/data/cbNIPT/processed_dataset/cbNIPT_25_03_01_DS/bam/cbNIPT_25_03_01_DS.analysisReady.bam
BamFile2=/data/cbNIPT/processed_dataset/cbNIPT_25_03_02_DS/bam/cbNIPT_25_03_02_DS.analysisReady.bam

/storage/apps/HipSTR-0.7/HipSTR \
    --bams ${BamFile1},${BamFile2} \
    --fasta /storage/references_and_index/hg38/fasta/cbNIPT/hg38.fa \
    --regions /storage/home/kangsm/myDB/STR_references/HipSTR/hg38.hipstr_reference.bed \
    --str-vcf /storage/home/kangsm/analysis/gcx.cbNIPT/04.ANALYSIS_DEV/STR_analysis/HipSTR/HTR8.HipSTR.calls.vcf.gz \
    --min-reads 2 --def-stutter-model

/storage/apps/gangstr-2.5.0/bin/GangSTR \
    --bam ${BamFile1},${BamFile2} \
    --ref /storage/references_and_index/hg38/fasta/cbNIPT/hg38.fa \
    --regions /storage/home/kangsm/myDB/STR_references/GangSTR/hg38_ver13.bed \
    --out /storage/home/kangsm/analysis/gcx.cbNIPT/04.ANALYSIS_DEV/STR_analysis/HipSTR/HTR8.GangSTR.calls

