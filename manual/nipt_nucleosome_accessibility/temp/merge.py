import pandas as pd
from glob import glob 



path_list = glob('/storage/home/jhkim/Projects/NIPT/GCX-NIPT-260121/Results/cfDNA.MRD.IC.FragmentCount/*tsv')


total_list = []
index_list = []
for path in path_list:

    df = pd.read_csv(path, sep='\t')
    total_list.append(df)
    name = path.split('/')[-1].replace('.tsv','')
    #df.index= [name]
    #print(df)
    index_list.append(name)
    total_list.append(df)


total_df = pd.concat(total_list)
total_df.index = index_list
print(total_df)

df.to_csv('temp.csv',sep=',')