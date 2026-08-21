import os 
from glob import glob 

path = '/storage/home/jhkim/workspace/projects/NIPT/260729-GCX-NIPT-PieplineSetup/Results/PRJNA837410'

#path = '/storage/home/jhkim/workspace/projects/NIPT/260729-GCX-NIPT-PieplineSetup/Results/PRJNA385180'

bam_file_list = glob(f'{path}/*/bam/*.recal.filtered*.bam')

for bam_path in bam_file_list:

    bam_id = os.path.basename(bam_path).split('.')[0]
    #cmd = f'''/storage/home/jhkim/Apps/Python-3.11.13/python run.py cnv --bam {bam_path} --bw /storage/home/jhkim/Projects/NIPT/GCX-NIPT-260121/Resources/reference/hg38.100mer.bw \
    #    --fasta /storage/references_and_index/hg38/fasta/cbNIPT/hg38.fa \
    #    --outdir /storage/home/jhkim/workspace/projects/NIPT/260729-GCX-NIPT-PieplineSetup/Analysis/WPS/{bam_id} \
    #    --sample {bam_id} \
    #    --mode L S \
    #    --marker-bed /storage/home/jhkim/workspace/projects/NIPT/260729-GCX-NIPT-PieplineSetup/Resources/local/refs/merged_marker.bed \
    #    --jobs 8'''
    #os.system(cmd)
    cmd = f'''/storage/home/jhkim/Apps/Python-3.11.13/python run.py cnv --bam {bam_path} --bw /storage/home/jhkim/Projects/NIPT/GCX-NIPT-260121/Resources/reference/hg38.100mer.bw \
        --fasta /storage/references_and_index/hg38/fasta/cbNIPT/hg38.fa \
        --outdir /storage/home/jhkim/workspace/projects/NIPT/260729-GCX-NIPT-PieplineSetup/Analysis/WPS/{bam_id} \
        --sample {bam_id}
        '''
    os.system(cmd)
    #cmd = f'''/storage/home/jhkim/Apps/Python-3.11.13/python run.py wps --bam {bam_path} \
    #    --outdir /storage/home/jhkim/workspace/projects/NIPT/260729-GCX-NIPT-PieplineSetup/Analysis/WPS/{bam_id} \
    #    --sample {bam_id} \
    #    --mode L S \
    #    --marker-bed /storage/home/jhkim/workspace/projects/NIPT/260729-GCX-NIPT-PieplineSetup/Resources/local/refs/merged_marker.bed \
    #    --jobs 8'''
    #os.system(cmd)
    #print(cmd)
