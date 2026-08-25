import pandas as pd
import os

from glob import glob

result_dict = {
    "sample_id":[],
    "group":[],
    "wps_dir":[],
}

wps_dir = '/storage/home/jhkim/workspace/projects/NIPT/260729-GCX-NIPT-PieplineSetup/Analysis/WPS.cnv'

sample_list = glob('/storage/home/jhkim/workspace/projects/NIPT/260729-GCX-NIPT-PieplineSetup/Results/PRJNA837410/S*')
for sample in sample_list:

    sample_id = sample.split('/')[-1]

    _wps_dir = f'{wps_dir}/{sample_id}/wps'

    result_dict['sample_id'].append(sample_id)
    result_dict['group'].append('Pragent_Abnormal')
    result_dict['wps_dir'].append(_wps_dir)


sample_list = glob('/storage/home/jhkim/workspace/projects/NIPT/260729-GCX-NIPT-PieplineSetup/Results/PRJNA385180/AnalysisUse/S*')
for sample in sample_list:

    sample_id = sample.split('/')[-1]

    _wps_dir = f'{wps_dir}/{sample_id}/wps'

    result_dict['sample_id'].append(sample_id)
    result_dict['group'].append('Normal')
    result_dict['wps_dir'].append(_wps_dir)


df = pd.DataFrame(result_dict)# .items())
print(df)

df.to_csv('sample.tsv',index=False,sep='\t')