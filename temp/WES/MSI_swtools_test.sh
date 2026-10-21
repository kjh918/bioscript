
#! /usr/bin/bash

#---/ 
    JOBS=1
    DATA_DIR=/BiO/DATA/WES/23WES02
    SAMPLE_LIST=$(cat ${DATA_DIR}/sample.list)
    SAMPLE_NUMBERS=$(echo ${SAMPLE_LIST} | wc -w)
    if [[ "${JOBS}" -eq "1"  &&  "${SAMPLE_NUMBERS}" -eq "1" ]]; then
        JOBS=${JOBS}
    else
        JOBS=${SAMPLE_NUMBERS}
    fi

#---| msisensor2 test |
    INPUT_BAM_DIR=${DATA_DIR}/bam_files
    MSISENSOR2_OUTPUT_DIR=${DATA_DIR}/extra_analysis/MSI/msisensor2_test

    /BiO/Install/msisensor2-0.1/msisensor2 scan \
        -d /BiO/Reference_and_Index/gcx.hg19/gcx.hg19.fasta \
        -o ${MSISENSOR2_OUTPUT_DIR}/reference.list

    parallel -j ${JOBS} -k "cp ${INPUT_BAM_DIR}/{}.sorted.dup.marked.recal.bai ${INPUT_BAM_DIR}/{}.sorted.dup.marked.recal.bam.bai" ::: ${SAMPLE_LIST}
    
    parallel -j ${JOBS} -k /BiO/Install/msisensor2-0.1/msisensor2 msi -l 1 -q 1 -b 2 \
        -M /BiO/Install/msisensor2-0.1/models_b37_HumanG1Kv37 \
        -t ${DATA_DIR}/bam_files/{}.sorted.dup.marked.recal.bam \
        -o ${MSISENSOR2_OUTPUT_DIR}/{}.msisensor2.model.only ::: ${SAMPLE_LIST}

    parallel -j ${JOBS} -k /BiO/Install/msisensor2-0.1/msisensor2 msi -l 1 -q 1 -b 2 \
        -d ${MSISENSOR2_OUTPUT_DIR}/reference.list \
        -t ${DATA_DIR}/bam_files/{}.sorted.dup.marked.recal.bam \
        -o ${MSISENSOR2_OUTPUT_DIR}/{}.msisensor2.ref.only ::: ${SAMPLE_LIST}

    parallel -j ${JOBS} -k /BiO/Install/msisensor2-0.1/msisensor2 msi -l 1 -q 1 -b 2 \
        -M /BiO/Install/msisensor2-0.1/models_b37_HumanG1Kv37 \
        -d ${MSISENSOR2_OUTPUT_DIR}/reference.list \
        -t ${DATA_DIR}/bam_files/{}.sorted.dup.marked.recal.bam \
        -o ${MSISENSOR2_OUTPUT_DIR}/{}.msisensor2.ref.model ::: ${SAMPLE_LIST}


/BiO/Install/msisensor2-0.1/msisensor2 msi \
    -l 1 -q 1 -b 2 \
    -d ${MSISENSOR2_OUTPUT_DIR}/reference.list \
    -t ${DATA_DIR}/bam_files/3-100X.sorted.dup.marked.recal.bam \
    -o ${MSISENSOR2_OUTPUT_DIR}/3-100X.msisensor2.score.test

echo -e "SampleID \t Total_Number_of_Sites \t Number_of_Somatic_Sites \t %" >> ./summary.txt
echo -e "3-100X \t $(sed -n '2p' 3-100X.msisensor2.ref.only)" >> ./summary.txt
echo -e "3-200X \t $(sed -n '2p' 3-200X.msisensor2.ref.only)" >> ./summary.txt
echo -e "3-400X \t $(sed -n '2p' 3-400X.msisensor2.ref.only)" >> ./summary.txt
echo -e "17-100X \t $(sed -n '2p' 17-100X.msisensor2.ref.only)" >> ./summary.txt
echo -e "17-200X \t $(sed -n '2p' 17-200X.msisensor2.ref.only)" >> ./summary.txt
echo -e "17-400X \t $(sed -n '2p' 17-400X.msisensor2.ref.only)" >> ./summary.txt
echo -e "tissue-200X \t $(sed -n '2p' tissue-200X.msisensor2.ref.only)" >> ./summary.txt


#---/

#---| 





