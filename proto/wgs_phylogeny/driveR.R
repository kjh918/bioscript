

./prepare_annovar_user.pl -dbtype cosmic /storage/DB/cosmic/v103/GRCh37/


/storage/apps/annovar/convert2annovar.pl


bcftools view -s WGS-26-01_01 /data/wgs/WGS-26-01/WGS-26-01_01/vcf/WGS-26-01_01.NT.VARIANTS.VEP.ANNOTATED.vcf -Ou \
    | bcftools view -Ov | /storage/apps/annovar/convert2annovar.pl \
    -format vcf4  - > /data/wgs/WGS-26-01/WGS-26-01_01/vcf/WGS-26-01_01.NT.annovar.avinput

/storage/apps/annovar/table_annovar.pl \
    -buildver hg19 \
    -out /data/wgs/WGS-26-01/WGS-26-01_01/vcf/WGS-26-01_01.NT.annovar.output.csv \
    -remove \
    -protocol refGene,cytoBand,exac03,avsnp150,dbnsfp30a,avsift,cosmic91 \
    -operation gx,r,f,f,f,f \
    -nastring . \
    -csvout \
    -polish \
    /data/wgs/WGS-26-01/WGS-26-01_01/vcf/WGS-26-01_01.NT.annovar.avinput \
    /storage/apps/annovar/humandb/ 

    cd /storage/apps/CancerVar
    python3 /storage/apps/CancerVar/CancerVar.py \
        -i /data/wgs/WGS-26-01/WGS-26-01_01/vcf/WGS-26-01_01.NT.annovar.avinput \
        -c /storage/apps/CancerVar/config.ini \
        --input_type AVinput \
        --cancer_type Breast \
        -o /data/wgs/WGS-26-01/WGS-26-01_01/vcf/WGS-26-01_01.NT.cancervar.output \
        -t /storage/apps/CancerVar/cancervardb

    python3 /storage/apps/CancerVar/OPAI/scripts/feature_preprocess.py \
        -a /data/wgs/WGS-26-01/WGS-26-01_01/vcf/WGS-26-01_01.NT.cancervar.output.hg19_multianno.txt.grl_p \
        -c /data/wgs/WGS-26-01/WGS-26-01_01/vcf/WGS-26-01_01.NT.cancervar.output.hg19_multianno.txt.cancervar \
        -m ensemble \
        -n 5 \
        -d /storage/apps/CancerVar/OPAI/saves/nonmissing_db.npy \
        -o /data/wgs/WGS-26-01/WGS-26-01_01/vcf/WGS-26-01_01.NT.cancervar.output.hg19_multianno.cancervar.ensemble.csv

    python3 /storage/apps/CancerVar/OPAI/scripts//opai_predictor.py \
        -i /data/wgs/WGS-26-01/WGS-26-01_01/vcf/WGS-26-01_01.NT.cancervar.output.hg19_multianno.cancervar.ensemble.csv \
        -v /data/wgs/WGS-26-01/WGS-26-01_01/vcf/WGS-26-01_01.NT.cancervar.output.hg19_multianno.txt.cancervar \
        -m ensemble -c /storage/apps/CancerVar/OPAI/saves/ensemble.pt \
        -d cpu \
        -o /data/wgs/WGS-26-01/WGS-26-01_01/vcf/WGS-26-01_01.NT.cancervar.output.hg19_multianno.txt.cancervar.ensemble.pred


