


mkdir -p viral_ref
cd viral_ref

# HPV16
curl -L "https://www.ncbi.nlm.nih.gov/search/api/download-sequence/?db=nuccore&id=NC_001526.4" \
  | sed '1s/.*/>HPV16|NC_001526.4/' > HPV16.fa

# HPV18
curl -L "https://www.ncbi.nlm.nih.gov/search/api/download-sequence/?db=nuccore&id=NC_001357.1" \
  | sed '1s/.*/>HPV18|NC_001357.1/' > HPV18.fa

# HBV
curl -L "https://www.ncbi.nlm.nih.gov/search/api/download-sequence/?db=nuccore&id=NC_003977.2" \
  | sed '1s/.*/>HBV|NC_003977.2/' > HBV.fa

# EBV
curl -L "https://www.ncbi.nlm.nih.gov/search/api/download-sequence/?db=nuccore&id=NC_007605.1" \
  | sed '1s/.*/>EBV|NC_007605.1/' > EBV.fa

# MCPyV
curl -L "https://www.ncbi.nlm.nih.gov/search/api/download-sequence/?db=nuccore&id=NC_010277.2" \
  | sed '1s/.*/>MCPyV|NC_010277.2/' > MCPyV.fa

# HHV-8 / KSHV
curl -L "https://www.ncbi.nlm.nih.gov/search/api/download-sequence/?db=nuccore&id=NC_009333.1" \
  | sed '1s/.*/>HHV8|NC_009333.1/' > HHV8.fa

# 하나로 합치기
cat HPV16.fa HPV18.fa HBV.fa EBV.fa MCPyV.fa HHV8.fa > viral.fa