import os
import pysam
import numpy as np
import pandas as pd
import sys
sys.path.append(os.path.dirname(__file__))  # 현재 디렉토리를 모듈 검색 경로에 추가

from visualization import plot_gc_correction, plot_genome_wide_cnv_with_annotations

sample_id = "cbNIPT_24_04_03_DS"

path = f'/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Results/manual/{sample_id}/data/{sample_id}.segments.tsv'

df_global = pd.read_csv(f"/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Results/manual/{sample_id}/data/{sample_id}.normalized.tsv", sep="\t")

df = pd.read_csv(path,sep='\t')
plot_genome_wide_cnv_with_annotations(
            bins_df=df_global, 
            segments_df=df, 
            output_path=os.getcwd() + 'temp.png',
            run_id=f'{sample_id}'
        )
print(df)