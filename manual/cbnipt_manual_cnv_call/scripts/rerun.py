
import os 
import sys 
from glob import glob 


bash_script_path = glob('/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Results/temp/*/logs/02.cnv_call/*sh')


for bash in bash_script_path:
    print(f"bash: {bash}")
    os.system(f"bash {bash}")

