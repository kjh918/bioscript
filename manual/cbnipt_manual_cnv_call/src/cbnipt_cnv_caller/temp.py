import pandas as pd 
from glob import glob 
import os
import sys


path = '/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Analysis/cbNIPT_24_04_04_2X_s1235/data'

df = pd.read_csv(f'{path}/cbNIPT_24_04_04_2X_s1235.normalized.new.tsv', sep='\t')

df = df[['chrom','start','end','bin_id','bin_BAF','copy_number_signal']]

print(df)