

# 
    BaseDir=/storage/home/kangsm/analysis/gcx.rnaseq.maqc
    BatchID=no_rep_subsampling_test

#
    ln -s /data/wts/WTS_MAQC/MAQC_II_A_rep1/fastq/MAQC_II_A_rep1.trimmed_R1.fastq.gz ${BaseDir}/${BatchID}/fastq/MAQC_II_NoRep_A_R1.fastq.gz
    ln -s /data/wts/WTS_MAQC/MAQC_II_A_rep1/fastq/MAQC_II_A_rep1.trimmed_R2.fastq.gz ${BaseDir}/${BatchID}/fastq/MAQC_II_NoRep_A_R2.fastq.gz
    ln -s /data/wts/WTS_MAQC/MAQC_II_B_rep1/fastq/MAQC_II_B_rep1.trimmed_R1.fastq.gz ${BaseDir}/${BatchID}/fastq/MAQC_II_NoRep_B_R1.fastq.gz
    ln -s /data/wts/WTS_MAQC/MAQC_II_B_rep1/fastq/MAQC_II_B_rep1.trimmed_R2.fastq.gz ${BaseDir}/${BatchID}/fastq/MAQC_II_NoRep_B_R2.fastq.gz
    

### 75% sub-sampling
    SampleID_List="MAQC_II_NoRep_A MAQC_II_NoRep_B"

    parallel -j 10 -k "/storage/apps/seqkit-2.8.2/seqkit sample -p 0.75 -s 1357 ${BaseDir}/${BatchID}/fastq/{}_R1.fastq.gz -o ${BaseDir}/${BatchID}/fastq/{}_S75_rep1_R1.fastq.gz" ::: $SampleID_List
    parallel -j 10 -k "/storage/apps/seqkit-2.8.2/seqkit sample -p 0.75 -s 1357 ${BaseDir}/${BatchID}/fastq/{}_R2.fastq.gz -o ${BaseDir}/${BatchID}/fastq/{}_S75_rep1_R2.fastq.gz" ::: $SampleID_List
    parallel -j 10 -k "/storage/apps/seqkit-2.8.2/seqkit sample -p 0.75 -s 2468 ${BaseDir}/${BatchID}/fastq/{}_R1.fastq.gz -o ${BaseDir}/${BatchID}/fastq/{}_S75_rep2_R1.fastq.gz" ::: $SampleID_List
    parallel -j 10 -k "/storage/apps/seqkit-2.8.2/seqkit sample -p 0.75 -s 2468 ${BaseDir}/${BatchID}/fastq/{}_R2.fastq.gz -o ${BaseDir}/${BatchID}/fastq/{}_S75_rep2_R2.fastq.gz" ::: $SampleID_List


### 
    qsub -q all.q@ngsnode1 /storage/home/kangsm/myScripts/WTS/WTS_MAQC/MAQC_II_SGE_Subsample_Processing.sh


