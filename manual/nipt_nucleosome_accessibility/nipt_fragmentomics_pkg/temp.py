import os 
from glob import glob 

path = '/storage/home/jhkim/workspace/test/projects/NIPT/260729-GCX-NIPT-PieplineSetup/Results/PRJNA837410'

#path = '/storage/home/jhkim/workspace/test/projects/NIPT/260729-GCX-NIPT-PieplineSetup/Results/PRJNA385180'

bam_file_list = glob(f'{path}/*/bam/*.recal.filtered*.bam')

for bam_path in bam_file_list:
    print(bam_path)
    bam_id = os.path.basename(bam_path).split('.')[0]
    cmd = f'''/storage/home/jhkim/Apps/Python-3.11.13/python -m run \
    --bam {bam_path} \
    --bin-bed  /storage/home/jhkim/Projects/NIPT/GCX-NIPT-260121/Resources/reference/Binning/hg38.fa.bin_100.0K.bed.gz \
    --out-dir  /storage/home/jhkim/workspace/test/projects/NIPT/260729-GCX-NIPT-PieplineSetup/Analysis/WPS/{bam_id} \
    --marker-bed /storage/home/jhkim/workspace/test/projects/NIPT/260729-GCX-NIPT-PieplineSetup/Resources/local/refs/merged_marker.bed \
    --fasta    /storage/references_and_index/hg38/fasta/cbNIPT/hg38.fa \
    --vcf    /storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Resources/reference/KOVA_v7/kova_sites_vcf/KOVA_v7_merged.vcf.gz \
    --bw       /storage/home/jhkim/Projects/NIPT/GCX-NIPT-260121/Resources/reference/hg38.100mer.bw \
    --wps-win 1000 --jobs 8 \
    --min-mappability 0.9 \
    --sample-id {bam_id} --resume'''

    os.system(cmd)
    #exit()