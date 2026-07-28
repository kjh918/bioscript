import os, sys 
from glob import glob 



file_list = glob('/storage/home/jhkim/Projects/NIPT/GCX-NIPT-260121/Results/WGS.cfDNA.ControlSamples/*/bam/*.recal.filtered.bam')


for path in file_list:
    print(path)
    sample_id = path.split('/')[-3]
    cmd = f'python3.11 check_fragment_count.py --bam {path} --bed /storage/home/jhkim/Projects/NIPT/GCX-NIPT-260121/Resources/mrd_ic.bed --out /storage/home/jhkim/Projects/NIPT/GCX-NIPT-260121/Results/{sample_id}.tsv --short-max 150 --mono-min 151 --mono-max  220 --min-mapq 20'
    os.system(cmd)