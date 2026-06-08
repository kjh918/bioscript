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
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from bs4 import BeautifulSoup

# [MODIFIED] Added robust session generator to prevent connection aborts during scraping
def get_robust_session() -> requests.Session:
    session = requests.Session()
    retries = Retry(
        total=5,
        backoff_factor=1, 
        status_forcelist=[500, 502, 503, 504],
        allowed_methods=["HEAD", "GET", "OPTIONS"]
    )
    adapter = HTTPAdapter(max_retries=retries)
    session.mount("http://", adapter)
    session.mount("https://", adapter)
    return session

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
        raise RuntimeError(f"Failed to parse XML string into dictionary using xmltodict : {e}")
    
def _decode_xml(payload: bytes) -> str:
    try:
        return payload.decode("utf-8")
    except UnicodeDecodeError:
        return payload.decode("latin-1", errors="replace")

def fetch_related_literature(uid: str, max_papers: int = 10) -> List[Dict[str, str]]:
    """ELink를 사용하여 MedGen ID에 연결된 PubMed 논문 정보 추출"""
    literature_list = []
    try:
        with Entrez.elink(dbfrom="medgen", db="pubmed", id=uid) as handle:
            link_record = Entrez.read(handle)
            
        pmids = []
        for dblist in link_record[0].get("LinkSetDb", []):
            if dblist.get("DbTo") == "pubmed":
                for link in dblist.get("Link", []):
                    pmids.append(link["Id"])
                    
        if not pmids:
            return []
            
        target_pmids = pmids[:max_papers]
        
        with Entrez.esummary(db="pubmed", id=",".join(target_pmids)) as handle:
            pub_records = Entrez.read(handle)
            
        for article in pub_records:
            literature_list.append({
                "pmid": str(article["Id"]),
                "title": article.get("Title", "No title available")
            })
    except Exception as e:
        print(f"  -> [Warning] UID {uid}의 논문 정보 추출 실패: {e}")
    return literature_list

def extract_medgen_section(soup: BeautifulSoup, keywords: List[str]) -> str:
    """
    HTML 문서 내에서 특정 키워드가 포함된 헤딩(h2, h3 등)을 찾아 
    그 아래에 있는 텍스트 콘텐츠를 병합하여 반환합니다.
    """
    for heading in soup.find_all(["h2", "h3", "h4", "dt"]):
        if heading.text and any(k.lower() in heading.text.lower() for k in keywords):
            content = []
            current_node = heading.find_next_sibling()
            
            # 다음 헤딩이 나오기 전까지의 모든 텍스트 추출
            while current_node and current_node.name not in ["h1", "h2", "h3", "h4", "dt"]:
                text = current_node.get_text(separator=" ", strip=True)
                if text and text not in content:
                    content.append(text)
                current_node = current_node.find_next_sibling()
                
            return " | ".join(content) if content else "No content found"
            
    return "Not provided"

def extract_term_hierarchy(soup: BeautifulSoup) -> str:
    """Term Hierarchy 섹션에서 뱃지(C,R,O,G,V)를 제외하고 계층 구조만 추출"""
    hierarchy = []
    header = soup.find(lambda tag: tag.name in ["h2", "h3", "h4", "dt"] and "term hierarchy" in tag.text.lower())
    
    if header:
        current_node = header.find_next_sibling()
        while current_node and current_node.name not in ["h1", "h2", "h3", "h4", "dt"]:
            tl_lines = current_node.find_all("span", class_="TLline")
            direct_text = current_node.find(string=True, recursive=False)
            
            if direct_text and direct_text.strip():
                hierarchy.append(direct_text.strip())
                
            for line in tl_lines:
                a_tag = line.find("a")
                if a_tag:
                    hierarchy.append(a_tag.text.strip())
                else:
                    hierarchy.append(line.text.strip())
            current_node = current_node.find_next_sibling()

    clean_hierarchy = []
    for item in hierarchy:
        clean_item = " ".join(item.split())
        if clean_item and clean_item not in clean_hierarchy:
            clean_hierarchy.append(clean_item)
            
    return " > ".join(clean_hierarchy) if clean_hierarchy else "Not provided"

def extract_ds_tree_to_dict(html_content: str) -> dict:
    """
    Parses the nested 'ds_tree' HTML structure and extracts the labels, 
    attributes, and hierarchical relationships into a dictionary.
    """
    soup = BeautifulSoup(html_content, "html.parser")
    
    # Find the main container
    tree_container = soup.find("div", class_="ds_tree")
    if not tree_container:
        return {"error": "Could not find 'ds_tree' class in the HTML."}
        
    # Find the root unordered list
    root_ul = tree_container.find("ul")
    if not root_ul:
        return {"error": "Could not find root 'ul' inside 'ds_tree'."}

    # Recursive function to parse nested lists
    def _parse_ul(ul_element):
        nodes = []
        
        # recursive=False ensures we only process direct children, preventing deep overlap
        for li in ul_element.find_all("li", recursive=False):
            node_data = {}
            
            # Find the span containing the data
            tl_line = li.find("span", class_="TLline")
            
            if tl_line:
                a_tag = tl_line.find("a")
                if a_tag:
                    # Extract text and all attributes (href, data-ga-label, etc.)
                    node_data["label"] = a_tag.get_text(strip=True)
                    node_data["attributes"] = a_tag.attrs
                else:
                    # Fallback if there is no <a> tag (just raw text)
                    node_data["label"] = tl_line.get_text(strip=True)
                    node_data["attributes"] = {}
            
            # Check for nested <ul> representing children
            child_ul = li.find("ul", recursive=False)
            if child_ul:
                node_data["children"] = _parse_ul(child_ul)
            else:
                node_data["children"] = []
                
            nodes.append(node_data)
            
        return nodes

    # Start the recursive parsing from the root
    parsed_tree = {
        "tree_hierarchy": _parse_ul(root_ul)
    }
    
    return parsed_tree

def fetch_medgen_details(disease_name: str, max_results: int) -> List[Dict[str, Any]]:
    """
    Searches the MedGen database by disease name to get UIDs, 
    and extracts full text definitions by scraping the target URLs.
    """
    if not disease_name: 
        raise ValueError("disease_name parameter is missing.")
    if max_results is None or max_results <= 0: 
        raise ValueError("max_results must be 1 or greater.")

    clean_disease_name = urllib.parse.unquote(disease_name).strip()
    if not clean_disease_name:
        raise ValueError("Decoded disease_name is empty.")

    search_term = f'"{clean_disease_name}"'
    
    print(f"  -> [Debug] NCBI MedGen ESearch query: {search_term}")

    # 1. ESearch (Get MedGen UIDs)
    try:
        with Entrez.esearch(db="medgen", term=search_term, retmax=max_results, sort="relevance") as handle:
            record = Entrez.read(handle)
        id_list = record.get("IdList", [])
    except Exception as e:
        raise RuntimeError(f"MedGen esearch failed: {e}")

    print(f"  -> [Debug] ESearch retrieved MedGen UIDs: {len(id_list)}")

    if not id_list:
        print("  -> [Warning] 0 results found. Returning empty list.")
        return []

    time.sleep(0.35)

    medgen_data = []
    session = get_robust_session()
    
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36",
        "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8",
        "Accept-Language": "en-US,en;q=0.5",
        "Connection": "keep-alive",
    }
    
    for uid in id_list:
        url = f"https://www.ncbi.nlm.nih.gov/medgen/{uid}"
        print(f"  -> [Scraping] Fetching full contents for UID: {uid}...")
        
        try:
            response = session.get(url, headers=headers, timeout=15)
            response.raise_for_status()
        except Exception as e:
            print(f"  -> [Error] Failed to fetch {url}: {e}")
            continue

        soup = BeautifulSoup(response.text, "html.parser")
        
        # [HTML 구조 확인용 - 파싱 테스트가 끝나면 아래 4줄을 삭제/주석처리 하세요]
        #print("\n" + "="*60)
        print(f"[DEBUG] UID {uid}의 전체 HTML 데이터 구조:")
        print("="*60)
        print(soup.prettify()) 
        print("="*60 + "\n")
        exit()
        target_div = soup.find("div", class_='portlet_content')
        #print(target_div)
        #exit()
        # 4. 텍스트 추출 및 정제
        # separator=" " : <br>이나 <p> 같은 태그를 공백 하나로 치환하여 단어가 붙어버리는 것을 방지합니다.
        # 줄바꿈으로 나누고 싶다면 separator="\n" 을 사용하세요.
        # strip=True : 텍스트 앞뒤의 불필요한 공백과 줄바꿈을 모두 제거합니다.
        definition_text = target_div.get_text(separator=" ", strip=True)
        target_div = soup.find("h1", class_='nl',id='Term_Hierarchy')
        definition_text = target_div.get_text(separator=" ", strip=True)
        print(definition_text)
        exit() # HTML 구조 확인 후 파싱을 진행하려면 이 줄을 삭제하세요.
        
        '''
        <div class="portlet_content ln">


        '''
        # --- 새로운 파싱 로직 적용 ---
        
        # 1. 제목 (MedGenTitleText 우선 탐색)
        title_div = soup.find("div", class_="MedGenTitleText")
        if not title_div:
            title_div = soup.find("h1")
        disease_title = title_div.text.strip() if title_div else "Unknown Disease"
        
        # 2. 텍스트 섹션 
        definition = extract_medgen_section(soup, ["definition"])
        synonyms = extract_medgen_section(soup, ["synonym"])
        
        # 3. Term Hierarchy
        term_hierarchy = extract_term_hierarchy(soup)
        
        # 4. 연관 문헌
        literature = fetch_related_literature(uid, max_papers=10)

        # 5. 결과 저장
        parsed_record = {
            "medgen_uid": uid,
            "disease_name": disease_title,
            "definition": definition,
            "synonyms": synonyms,
            "term_hierarchy": term_hierarchy,
            "literature": literature,
            "url": url
        }
        
        medgen_data.append(parsed_record)
        time.sleep(1.5)

    return medgen_data


def main():
    parser = argparse.ArgumentParser(description="Fetch Full Disease Information and Texts from NCBI MedGen")
    parser.add_argument("--disease", required=True, help="Disease Name (e.g., Achondroplasia)")
    parser.add_argument("--max-results", type=int, default=5, help="Number of TOP results to fetch (default: 5)")
    parser.add_argument("--output", type=str, default="medgen_results.json", help="Output JSON file name")

    args = parser.parse_args()

    # Configure NCBI API email globally
    Entrez.email = "your.email@example.com"  # Remember to replace this with your actual email

    try:
        print(f"▶ Starting MedGen information search for '{args.disease}'...")
        results = fetch_medgen_details(
            disease_name=args.disease, 
            max_results=args.max_results
        )
        
        print(f"\n▶ Search complete: {len(results)} items found")
        
        # Save output to JSON
        with open(args.output, "w", encoding="utf-8") as json_file:
            json.dump(results, json_file, ensure_ascii=False, indent=4)
            
        print(f"▶ Results successfully saved to: {args.output}")
        
    except Exception as e:
        print(f"\n[Fatal Error] {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()