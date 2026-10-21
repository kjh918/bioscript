


NEW_DIR=/data/wes/23WES02
OLD_DIR=/data/wes/23WES02_old
touch variantNumberCompare.txt
echo -e "Sample\trelaign-no-vcf\trealign-yes-vcf\trealign-no-filtered-vcf\trealign-yes-filtered-vcf" >> /data/wes/23WES02/meta/VCF-Compare2.txt


for seqID in "3-100X" "3-200X" "3-400X" "17-100X" "17-200X" "17-400X" "tissue-200X"
do
    VCF_DIR=/data/wes/23WES02/${seqID}/vcf
    # old vcf (no-local-realign, no-filtered)
    oldVcfNoFilterSkip=$(cat ${VCF_DIR}/${seqID}.mutect2.vcf | grep ^# | wc -l)
    oldVcfVariants=$(sed 1,${oldVcfNoFilterSkip}d ${VCF_DIR}/${seqID}.mutect2.vcf | wc -l)

    # new vcf (yes-local-realign, no-filtered)
    newVcfNoFilterSkip=$(cat ${VCF_DIR}/${seqID}.realign.mutect2.vcf | grep ^# | wc -l)
    newVcfVariants=$(sed 1,${newVcfNoFilterSkip}d ${VCF_DIR}/${seqID}.realign.mutect2.vcf | wc -l)

    # old vcf (no-local-realign, no-filtered)
    oldVcfFilterSkip=$(cat ${VCF_DIR}/${seqID}.mutect2.filtered.vcf | grep ^# | wc -l)
    oldVcfVariants2=$(sed 1,${oldVcfFilterSkip}d ${VCF_DIR}/${seqID}.mutect2.filtered.vcf | grep PASS | wc -l)

    # new vcf (yes-local-realign, no-filtered)
    newVcfFilterSkip=$(cat ${VCF_DIR}/${seqID}.realign.mutect2.filtered.vcf | grep ^# | wc -l)
    newVcfVariant2s=$(sed 1,${newVcfFilterSkip}d ${VCF_DIR}/${seqID}.realign.mutect2.filtered.vcf | grep PASS | wc -l)

    echo -e "${seqID}\t${oldVcfVariants}\t${newVcfVariants}\t${oldVcfVariants2}\t${newVcfVariant2s}" >> /data/wes/23WES02/meta/VCF-Compare2.txt

done





