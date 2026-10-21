

#
    git clone https://github.com/AmpliconSuite/AmpliconSuite-pipeline
    sudo cp -r AmpliconSuite-pipeline /storage/apps/
    cd /storage/apps/AmpliconSuite-pipeline/
    ./install.sh --finalize_only --data_repo_loc /storage/references_and_index/hg38/amplicon_architect/GRCh38

#---| Manual Params --------------------------------------------------------------------------------
    BaseDir="/data/wgs"
    BatchID="cancer_cell_lines"
    SeqID="SRR8652105_MCF7"
    Threads=16
    GenomeAssembly="hg38"

#---| Input Files and Folders ----------------------------------------------------------------------
    AmpArchResDir=${BaseDir}/${BatchID}/${SeqID}/ecDNA/ampliconarchitect
    mkdir -p ${AmpArchResDir}
    cd ${AmpArchResDir}
    InputBamFile=${BaseDir}/${BatchID}/${SeqID}/bam/${SeqID}.analysisReady.bam

    FastqR1=${BaseDir}/${BatchID}/${SeqID}/fastq/${SeqID}.trimmed_R1.fastq.gz
    FastqR2=${BaseDir}/${BatchID}/${SeqID}/fastq/${SeqID}.trimmed_R2.fastq.gz

    

#---| Swtools Run ----------------------------------------------------------------------------------
    #RunAmpArch="singularity exec -B /storage,/data /storage/images/ampliconsuite-pipeline.sif /home/programs/AmpliconSuite-pipeline-master/AmpliconSuite-pipeline.py"
    export MOSEKLM_LICENSE_FILE=/storage/home/kangsm/myDB/Mosek_Lisence/Trial_20260309/mosek.lic
   
    if [ "${GenomeAssembly}" = "hg19" ]; then 
        AA_ReferenceVersion="GRCh37"
    else
        AA_ReferenceVersion="GRCh38"
    fi
    export AA_DATA_REPO=/storage/references_and_index/${GenomeAssembly}/amplicon_architect

    AA_SIF=/storage/images/ampliconsuite-pipeline_v1.5.0.sif

    python3 /storage/apps/AmpliconSuite-pipeline/singularity/run_paa_singularity.py \
        --sif ${AA_SIF} \
        --output_directory ${AmpArchResDir} \
        --sample_name ${SeqID} \
        --nthreads ${Threads} \
        --ref ${AA_ReferenceVersion} \
        --bam ${InputBamFile} \
        --run_AA \
        --run_AC 

