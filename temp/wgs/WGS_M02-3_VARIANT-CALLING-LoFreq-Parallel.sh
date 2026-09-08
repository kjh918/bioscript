
    TumorSeqID="${TumorSeqID}"
    TumorBamFile="${TumorBamDir}/${TumorSeqID}.analysisReady.bam"
    NormalBamFile="${NormalBamDir}/${NormalSeqID}.analysisReady.bam"
    LofreqRunDir="${LofreqDir_NT}"
    GenomeAssembly="hg19"
    TargetBed="${WgsCallRegionBed}"
    Threads=16
        
    RunLofreq="singularity exec -B /storage,/data /storage/images/lofreq-v2.1.sif lofreq"
    GenomeFasta="/storage/references_and_index/${GenomeAssembly}/fasta/${GenomeAssembly}.primary_decoy_HLA_no.ALT.fa"
    ChromTmpDir="${LofreqRunDir}/chr_tmp"
    mkdir -p "${ChromTmpDir}/final_vcfs"

    TargetChrom="1"


    # TARGET CHROM BED    
    TargetChromBed="${ChromTmpDir}/target.chrom.${TargetChrom}.bed"
    awk -v chr="${TargetChrom}" '$1==chr' "${TargetBed}" > "${TargetChromBed}"

    # TARGET-CHROMOSOME BAM
    samtools view -bh -@ "${Threads}" "${TumorBamFile}" "${TargetChrom}" > "${ChromTmpDir}/tumor.${TargetChrom}.bam"
    samtools view -bh -@ "${Threads}" "${NormalBamFile}" "${TargetChrom}" > "${ChromTmpDir}/normal.${TargetChrom}.bam"
    samtools index -@ "${Threads}" "${ChromTmpDir}/tumor.${TargetChrom}.bam"
    samtools index -@ "${Threads}" "${ChromTmpDir}/normal.${TargetChrom}.bam"
 
    # INDELQUAL
    ${RunLofreq} indelqual --dindel -f "${GenomeFasta}" "${ChromTmpDir}/tumor.${TargetChrom}.bam" -o "${ChromTmpDir}/tumor.${TargetChrom}.indelqual.bam"
    ${RunLofreq} indelqual --dindel -f "${GenomeFasta}" "${ChromTmpDir}/normal.${TargetChrom}.bam" -o "${ChromTmpDir}/normal.${TargetChrom}.indelqual.bam"
    samtools index -@ "${Threads}" "${ChromTmpDir}/tumor.${TargetChrom}.indelqual.bam"
    samtools index -@ "${Threads}" "${ChromTmpDir}/normal.${TargetChrom}.indelqual.bam"

    # CALL SOMATIC VARIANTS
    ${RunLofreq} somatic \
        -n "${ChromTmpDir}/normal.${TargetChrom}.indelqual.bam" \
        -t "${ChromTmpDir}/tumor.${TargetChrom}.indelqual.bam" \
        -f "${GenomeFasta}" -l "${TargetChromBed}" \
        --call-indels --threads "${Threads}" \
        -o "${ChromTmpDir}/${TumorSeqID}.target.chrom.${TargetChrom}."

    # COLLECT FINAL VCFs
    cp ${ChromTmpDir}/${TumorSeqID}.target.chrom.${TargetChrom}.somatic_final.indels.vcf.gz  ${ChromTmpDir}/final_vcfs/
    cp ${ChromTmpDir}/${TumorSeqID}.target.chrom.${TargetChrom}.somatic_final.indels.vcf.gz.tbi ${ChromTmpDir}/final_vcfs/
    cp ${ChromTmpDir}/${TumorSeqID}.target.chrom.${TargetChrom}.somatic_final.snvs.vcf.gz  ${ChromTmpDir}/final_vcfs/
    cp ${ChromTmpDir}/${TumorSeqID}.target.chrom.${TargetChrom}.somatic_final.snvs.vcf.gz.tbi  ${ChromTmpDir}/final_vcfs/

    # REMVOE TMP FILES
    rm ${ChromTmpDir}/${TumorSeqID}.target.chrom.${TargetChrom}.*
    rm ${ChromTmpDir}/tumor.${TargetChrom}.*
    rm ${ChromTmpDir}/normal.${TargetChrom}.*
    rm ${ChromTmpDir}/target.chrom.${TargetChrom}.bed



