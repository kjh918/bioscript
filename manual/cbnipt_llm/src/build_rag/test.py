import pandas as pd 
import os 

gcx_cbnipt_target_disease_path = '/storage/home/jhkim/scripts/bioscript/manual/cbnipt_llm/src/resources/gcx_nipt_target.tsv'

gcx_cbnipt_target_disease = pd.read_csv(gcx_cbnipt_target_disease_path, sep='\t')

target_disease = list(set(gcx_cbnipt_target_disease['SYNDROME']))

#target_disease = [
#    'Double Y syndrome',
#    'Trisomy 9'
#]
for disease in target_disease:
    
    output = disease.replace(' ','_') + '.json'

    cmd = f'python3.11 medgen_parser.py --disease "{disease}" --max-results 1 --email you@example.com --output {output}'

    cmd = f"python3.11 pubmed_parser.py --input /storage/home/jhkim/scripts/bioscript/manual/cbnipt_llm/database/medgen/{disease.replace(' ','_')}.json --output /storage/home/jhkim/scripts/bioscript/manual/cbnipt_llm/database/pubmed/{disease.replace(' ','_')}.pubmed.json"
    #print(cmd)
#
    os.system(cmd)








