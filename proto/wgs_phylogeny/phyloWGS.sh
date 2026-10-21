




    WorkDir="/data/wgs/WGS-26-01/analysis/phylowgs"
    enomeFasta="/storage/references_and_index/hg19/hmftools/Homo_sapiens.GRCh37.GATK.illumina.fasta"

    # vcf
    TumorList=($(cat /data/wgs/WGS-26-01/meta/analysis.list | grep -v "normal" | awk '{print $1}'))
    for TID in ${TumorList[@]}
    do
        hmftools_purple_vcf=${BaseDir}/${BatchID}/${TID}/hmftools/purple/${TID}.purple.somatic.vcf.gz

        bcftools norm -m -any -f "${GenomeFasta}" "${hmftools_purple_vcf}" \
        | bcftools view -f PASS -s ${TID} \
            -Oz -o ${WorkDir}/${TID}.purple.somatic.filtered.vcf.gz
        tabix -p vcf ${WorkDir}/${TID}.purple.somatic.filtered.vcf.gz
    done    

    # cnv
    #/storage/home/kangsm/runScripts/00_Pipeline_Development/wgs_phylogeny/create-phylowgs-cnv-data.R

    /storage/home/kangsm/runScripts/00_Pipeline_Development/wgs_phylogeny/select-phylowgs-input-variants.R
    
    # create phylowgs input

    RunInputCreator="singularity exec -B /storage,/data /storage/images/phylowgs.sif python /opt/conda/phylowgs/parser/create_phylowgs_inputs.py" 

    ${RunInputCreator} \
        --vcf-type T1=mutect_smchet \
        --vcf-type T2=mutect_smchet \
        --cnvs T1=/data/wgs/WGS-26-01/analysis/phylowgs/WGS-26-01_01.cnv.data.tsv \
        --cnvs T2=/data/wgs/WGS-26-01/analysis/phylowgs/WGS-26-01_02.cnv.data.tsv \
        T1=/data/wgs/WGS-26-01/analysis/phylowgs/WGS-26-01_01.phylowgs.input.filtered.vcf.gz \
        T2=/data/wgs/WGS-26-01/analysis/phylowgs/WGS-26-01_02.phylowgs.input.filtered.vcf.gz 

    # run phylowgs

    RunPhylowgs="singularity exec -B /storage,/data /storage/images/phylowgs.sif python /opt/conda/phylowgs/multievolve.py" 
    ${RunPhylowgs} \
        --num-chains 4 \
        
        --output-dir ${WorkDir}/results \
        --ssms ${WorkDir}/ssm_data.txt \
        --cnvs ${WorkDir}/cnv_data.txt 

    WritePhylowgs="singularity exec -B /storage,/data /storage/images/phylowgs.sif python /opt/conda/phylowgs/write_results.py" 
    ${WritePhylowgs} \
        "WGS-26-01" \
        ${WorkDir}/results/trees.zip \
        ${WorkDir}/WGS-26-01.summ.json.gz \
        ${WorkDir}/WGS-26-01.muts.json.gz \
        ${WorkDir}/WGS-26-01.mutass.zip

    cp -r /storage/home/kangsm/swtools/phylowgs/witness ${WorkDir}/

    gzip -dc ${WorkDir}/WGS-26-01.summ.json.gz > ${WorkDir}/witness/data/WGS-26-01.summ.json
    gzip -dc ${WorkDir}/WGS-26-01.muts.json.gz > ${WorkDir}/witness/data/WGS-26-01.muts.json

    