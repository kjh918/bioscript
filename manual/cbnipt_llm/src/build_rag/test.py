import pandas as pd 
import os 

gcx_cbnipt_target_disease_path = '/storage/home/jhkim/scripts/bioscript/manual/cbnipt_llm/src/resources/gcx_nipt_target.tsv'

gcx_cbnipt_target_disease = pd.read_csv(gcx_cbnipt_target_disease_path, sep='\t')

target_disease = list(set(gcx_cbnipt_target_disease['SYNDROME']))

#target_disease = [
#    'Double Y syndrome',
#    'Trisomy 9'
#]
target_disease = ['Williams syndrome']
for disease in target_disease:
    
    output = disease.replace(' ','_') + '.json'

    cmd = f'python3.11 medgen_parser.py --disease "{disease}" --max-results 1 --email you@example.com --output {output}'

    cmd = f"python3.11 pubmed_parser.py --input /storage/home/jhkim/scripts/bioscript/manual/cbnipt_llm/database/medgen/{disease.replace(' ','_')}.json --output /storage/home/jhkim/scripts/bioscript/manual/cbnipt_llm/database/pubmed/{disease.replace(' ','_')}.pubmed.json"
    #print(cmd)
#
    
    #cmd = f"python3.11 build_rag.py /storage/home/jhkim/scripts/bioscript/manual/cbnipt_llm/database/medgen/{disease.replace(' ','_')}.json "
    #cmd += f" /storage/home/jhkim/scripts/bioscript/manual/cbnipt_llm/database/pubmed/{disease.replace(' ','_')}.pubmed.json --output {disease.replace(' ','_')}.jsonl"
    #os.system(cmd)

    cmd = f"python3.11 chunk_builder.py --medgen /storage/home/jhkim/scripts/bioscript/manual/cbnipt_llm/database/medgen/{disease.replace(' ','_')}.json "
    cmd += f" --pubmed /storage/home/jhkim/scripts/bioscript/manual/cbnipt_llm/database/pubmed/{disease.replace(' ','_')}.pubmed.json --output {disease.replace(' ','_')}.json"
    
    resources_path = '/storage/home/jhkim/Projects/cbNIPT/260427-GCX-cbNIPT-LLM/Resources/temp'

    cmd = f'''python3.11 {resources_path}/syndrome_discovery.py --syndrome "{disease}" \
      --email you@example.com \
      --orphanet {resources_path}/en_product6.xml \
      --gencc {resources_path}/gencc_submissions.tsv \
      --clingen-gene {resources_path}/ClinGen_gene_curation_list_GRCh38.tsv \
      --clingen-region {resources_path}/ClinGen_region_curation_list_GRCh38.tsv \
      --hpo {resources_path}/phenotype_to_genes.txt \
      --tier Tier1_High,Tier2_Medium \
      --output-dir /storage/home/jhkim/Projects/cbNIPT/260427-GCX-cbNIPT-LLM/Resources/temp/output
      '''
    print(cmd)
    os.system(cmd)







