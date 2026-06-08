import sys
import time
import json
import argparse
import urllib.parse
import xml.etree.ElementTree as ET
from typing import List, Dict, Any
from Bio import Entrez
import xmltodict
# [MODIFIED] Imported libraries for HTTP requests and HTML parsing
import requests
from bs4 import BeautifulSoup
from pubmed import fetch_pubmed_details
import sys
import time
import json
import argparse
import urllib.parse
import xml.etree.ElementTree as ET
from typing import List, Dict, Any
from Bio import Entrez
import xmltodict
import requests
from bs4 import BeautifulSoup

def convert_element_to_dict(doc: ET.Element) -> dict:
    if doc is None:
        raise ValueError("The 'doc' element must be explicitly provided. None is not allowed.")
        
    try:
        doc_xml_str = ET.tostring(doc, encoding='unicode')
    except Exception as e:
        raise RuntimeError(f"Failed to convert XML element to string: {e}")

    if not doc_xml_str.strip():
        raise ValueError("The converted XML string is empty.")

    try:
        doc_dict = xmltodict.parse(doc_xml_str)
        return doc_dict
    except Exception as e:
        raise RuntimeError(f"Failed to parse XML string into dictionary using xmltodict: {e}")
    
def _decode_xml(payload: bytes) -> str:
    try:
        return payload.decode("utf-8")
    except UnicodeDecodeError:
        return payload.decode("latin-1", errors="replace")

def scrape_clinvar_page(var_id: str) -> Dict[str, Any]:
    url = f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{var_id}/"
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36"
    }
    
    try:
        response = requests.get(url, headers=headers, timeout=10)
        response.raise_for_status()
    except Exception as e:
        print(f"  -> [Warning] Failed to fetch HTML for ID {var_id}: {e}")
        return {"scraped_pubmed_ids": [], "error": str(e)}

    soup = BeautifulSoup(response.text, "html.parser")
    scraped_pubmed_ids = []
    
    for a_tag in soup.find_all("a", href=True):
        href = a_tag["href"]
        if "pubmed.ncbi.nlm.nih.gov/" in href:
            parts = href.strip("/").split("/")
            if parts and parts[-1].isdigit():
                pid = parts[-1]
                if pid not in scraped_pubmed_ids:
                    scraped_pubmed_ids.append(pid)

    return {
        "scraped_pubmed_ids": scraped_pubmed_ids
    }

# [MODIFIED] New function to calculate priority score based on severity and review status
def calculate_evidence_score(classification_dict: dict) -> int:
    """
    Calculates a numeric score to rank ClinVar variations.
    Higher score = more relevant/reliable evidence.
    """
    if not classification_dict:
        return 0
        
    score = 0
    description = classification_dict.get("description", "").lower()
    review_status = classification_dict.get("review_status", "").lower()

    # 1. Score by Severity (Pathogenicity)
    if "pathogenic" in description and "likely" not in description:
        score += 50
    elif "likely pathogenic" in description:
        score += 40
    elif "drug response" in description or "risk factor" in description:
        score += 30
    elif "benign" in description or "likely benign" in description:
        score += 10
    else: # Uncertain significance (VUS) or other
        score += 5

    # 2. Score by Review Status (ClinVar Stars)
    if "practice guideline" in review_status: # 4 stars
        score += 100
    elif "reviewed by expert panel" in review_status: # 3 stars
        score += 80
    elif "criteria provided, multiple submitters, no conflicts" in review_status: # 2 stars
        score += 60
    elif "criteria provided, single submitter" in review_status or "conflicting interpretations" in review_status: # 1 star
        score += 40
    else: # 0 stars (no assertion criteria provided)
        score += 0
        
    return score

def fetch_clinvar_evidence_by_disease(disease_name: str, max_results: int) -> List[Dict[str, Any]]:
    if not disease_name: raise ValueError("disease_name is missing.")
    if max_results is None or max_results <= 0: raise ValueError("max_results must be 1 or greater.")

    clean_disease_name = urllib.parse.unquote(disease_name).strip()
    if not clean_disease_name:
        raise ValueError("Decoded disease_name is empty.")

    search_term = f'"{clean_disease_name}"[Condition] OR "{clean_disease_name}"'
    
    # [MODIFIED] Oversample the ESearch. We grab more results than requested so we have a large pool to sort through.
    fetch_pool_size = max(max_results * 5, 50) 
    print(f"  -> [Debug] NCBI ESearch query: {search_term}")
    print(f"  -> [Debug] Oversampling pool size set to {fetch_pool_size} for optimal sorting.")

    try:
        with Entrez.esearch(db="clinvar", term=search_term, retmax=fetch_pool_size) as handle:
            record = Entrez.read(handle)
        id_list = record.get("IdList", [])
    except Exception as e:
        raise RuntimeError(f"ClinVar esearch failed: {e}")

    if not id_list:
        raise RuntimeError(f"0 results found. Query: '{search_term}'")

    time.sleep(0.35)

    try:
        with Entrez.esummary(db="clinvar", id=",".join(id_list)) as handle:
            raw_xml_bytes = handle.read()
            xml_text = _decode_xml(raw_xml_bytes)
            root = ET.fromstring(xml_text)
    except Exception as e:
        raise RuntimeError(f"ClinVar esummary call & XML parsing failed: {e}")
    
    raw_evidences = []
    
    for doc in root.findall(".//DocumentSummary"):
        var_id = doc.get("uid")
        if not var_id:
            continue
            
        dict_data = convert_element_to_dict(doc)
        doc_summary = dict_data.get("DocumentSummary")
        if not doc_summary:
            continue

        title_elem = doc.find("title")
        title = title_elem.text if title_elem is not None else "Unknown Title"

        accession = doc_summary.get("accession")
        
        # Safe extraction of nested elements
        variation_set = doc_summary.get("variation_set", {})
        variation = variation_set.get("variation") if isinstance(variation_set, dict) else None
        if not variation:
            continue

        variation_data = variation[0] if isinstance(variation, list) else variation
        variation_loc = variation_data.get("variation_loc", {})
        assembly_set_raw = variation_loc.get("assembly_set", [])
        assembly_set = assembly_set_raw[0] if isinstance(assembly_set_raw, list) else assembly_set_raw

        germline_classification = doc_summary.get("germline_classification")
        somatic_classification = doc_summary.get("somatic_classification")

        if not germline_classification and not somatic_classification:
            continue
            
        # [MODIFIED] Calculate the priority score based on the classification metadata
        target_classification = germline_classification if germline_classification else somatic_classification
        priority_score = calculate_evidence_score(target_classification)

        raw_evidences.append({
            "variation_id": var_id,
            "title": title,
            "accession": accession,
            "assembly_set": assembly_set,
            "germline_classification": germline_classification,
            "url": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{var_id}/",
            "priority_score": priority_score # Store score temporarily for sorting
        })

    # [MODIFIED] Sort the list descending by priority score and slice to get the best 'max_results'
    raw_evidences.sort(key=lambda x: x["priority_score"], reverse=True)
    top_evidences = raw_evidences[:max_results]
    
    print(f"\n  -> [Info] Successfully sorted and isolated the top {len(top_evidences)} high-quality variations.")

    final_evidences = []
    total_pmid_list = []

    # [MODIFIED] Now, loop ONLY through the top selected evidences to perform the slow web scraping
    for evidence in top_evidences:
        var_id = evidence["variation_id"]
        
        print(f"  -> [Scraping] Fetching HTML page details for Variation ID: {var_id} (Score: {evidence['priority_score']})...")
        scraped_data = scrape_clinvar_page(var_id)
        
        found_pmids = scraped_data.get("scraped_pubmed_ids", [])
        total_pmid_list.extend(found_pmids)
        
        evidence["scraped_pubmed_ids"] = found_pmids
        
        # Remove the temporary score key to keep JSON output clean
        del evidence["priority_score"] 
        
        final_evidences.append(evidence)
        time.sleep(1.0) # Polite delay for NCBI web servers

    # Deduplicate the total PMID list
    total_pmid_list = list(set(total_pmid_list))
    print(f"\n  -> [Summary] A total of {len(total_pmid_list)} unique PMIDs were gathered from the top results.")

    return final_evidences, total_pmid_list

def main():
    parser = argparse.ArgumentParser(description="Fetch ClinVar evidence by Disease Name (Raw XML Parsing & Web Scraping)")
    parser.add_argument("--email", required=True, help="NCBI API 이메일")
    parser.add_argument("--disease", required=True, help="질환명 (예: Achondroplasia)")
    parser.add_argument("--max-results", type=int, required=True, help="최대 검색 건수")

    args = parser.parse_args()

    if args.max_results <= 0:
        print("[Error] max-results는 1 이상이어야 합니다.")
        sys.exit(1)

    try:
        print(f"▶ '{args.disease}' 질환에 대한 ClinVar 근거 검색 시작...")
        results, pmid_list = fetch_clinvar_evidence_by_disease(
            #email=args.email, 
            disease_name=args.disease, 
            max_results=args.max_results
        )
        
        
        print(f"\n▶ 검색 완료: 총 {len(results)}건 확보")
        print(json.dumps(results, ensure_ascii=False, indent=2))
        
        x = fetch_pubmed_details(pmid_list)
        print(json.dumps(x, ensure_ascii=False, indent=2))


        with open("pubmed_articles.json", "w", encoding="utf-8") as f:
            for paper in x:
                f.write(json.dumps(paper, ensure_ascii=False) + "\n")

        with open("evidence.json", "w", encoding="utf-8") as f:
            for paper in results:
                f.write(json.dumps(paper, ensure_ascii=False) + "\n")


    except Exception as e:
        print(f"\n[Fatal Error] {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()