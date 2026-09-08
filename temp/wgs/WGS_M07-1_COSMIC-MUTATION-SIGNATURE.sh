


#---/ Manual Params /-------------------------------------------------------------------------------
    # BaseDir="/data/wgs"
    # BatchID="seqc2_wgs"
    # TumorSeqID="HCC1395"
    # NormalSeqID="HCC1395BL"
    # Threads=24
    # GenomeAssembly="hg38"
    # WithRnaseqData="true"
    # RnaseqHmftoolsDataDir="/data/wgs/seqc2_wgs/HCC1395_RNA/hmftools"
#---------------------------------------------------------------------------------------------------

#---| DATA and FOLDERS |----------------------------------------------------------------------------

    CosMutSigDir=${BaseDir}/${BatchID}/${TumorSeqID}/cosmic_mutsig
    mkdir -p ${CosMutSigDir}/input_data/

#---| RESOURCES |

    GenomeFasta=



    PurpleSomaticVCF=${BaseDir}/${BatchID}/${TumorSeqID}/hmftools/purple/${TumorSeqID}.purple.somatic.vcf.gz   
    
    if [ -f "${PurpleSomaticVCF}" ]; then
        bcftools view -f "PASS" ${PurpleSomaticVCF} -Oz -o ${CosMutSigDir}/input_data/${TumorSeqID}.purple.somatic.pass.vcf.gz
        tabix -p vcf ${CosMutSigDir}/input_data/${TumorSeqID}.purple.somatic.pass.vcf.gz
        bcftools norm -f ${GenomeFasta} -m -any ${CosMutSigDir}/input_data/${TumorSeqID}.purple.somatic.pass.vcf.gz \
            -Ov -o ${CosMutSigDir}/input_data/${TumorSeqID}.purple.somatic.pass.norm.vcf
        
    fi



    bcftools norm -f ${GenomeFasta} -m -any ${OutputVCF_Merged}.gz -Ov -o ${OutputVCF_Normalized}
