import sys
import json
import argparse
from typing import List, Dict, Any
from Bio import Entrez

# [MODIFIED] Changed the return type hint from Dict to List of Dicts
def fetch_pubmed_details(pmids: List[str]) -> List[Dict[str, Any]]:
    """
    Takes a list of PMIDs, queries the PubMed database, and parses the XML 
    to extract Journal, Article Title, Abstract text, and Publication Date.
    Returns a list of objects containing the pmid and its details.
    """
    if not pmids:
        # [MODIFIED] Return an empty list instead of an empty dictionary
        return []

    pmid_str = ",".join(pmids)
    
    # [MODIFIED] Initialize pubmed_data as a list
    pubmed_data = []

    try:
        with Entrez.efetch(db="pubmed", id=pmid_str, retmode="xml") as handle:
            records = Entrez.read(handle)
            
        for article in records.get('PubmedArticle', []):
            medline_citation = article.get('MedlineCitation', {})
            pmid = str(medline_citation.get('PMID', ''))
            
            article_info = medline_citation.get('Article', {})
            article_title = article_info.get('ArticleTitle', 'No title available')
            
            journal = article_info.get('Journal', {})
            journal_title = journal.get('Title', 'No journal available')
            
            journal_issue = journal.get('JournalIssue', {})
            pub_date_node = journal_issue.get('PubDate', {})
            
            year = pub_date_node.get('Year', '')
            month = pub_date_node.get('Month', '')
            day = pub_date_node.get('Day', '')
            medline_date = pub_date_node.get('MedlineDate', '')
            
            if medline_date:
                publish_date = medline_date
            else:
                date_parts = [part for part in (year, month, day) if part]
                publish_date = "-".join(date_parts) if date_parts else "No date available"
            
            abstract = article_info.get('Abstract', {}).get('AbstractText', [])
            abstract_text = " ".join([str(text) for text in abstract]) if abstract else "No abstract available."
            
            # [MODIFIED] Append the data as a dictionary object with "pmid" as a key
            pubmed_data.append({
                "pmid": pmid,
                "title": article_title,
                "journal": journal_title,
                "publish_date": publish_date,
                "abstract": abstract_text
            })
    except Exception as e:
        print(f"  -> [Warning] Failed to fetch PubMed details for {pmids}: {e}")

    return pubmed_data


def main():
    parser = argparse.ArgumentParser(description="Fetch PubMed details by PMID")
    parser.add_argument(
        "--pmid", 
        nargs="+", 
        required=True, 
        help="One or more PubMed IDs (e.g., --pmid 123456 789012)"
    )
    parser.add_argument(
        "--output",
        type=str,
        default="pubmed_results.json",
        help="Path to save the output JSON file (default: pubmed_results.json)"
    )

    args = parser.parse_args()

    # Be sure to replace this with your actual email address
    Entrez.email = "your.email@example.com"

    try:
        print(f"▶ Fetching details for {len(args.pmid)} PubMed ID(s)...")
        
        results = fetch_pubmed_details(pmids=args.pmid)
        
        print("\n▶ Search complete")
        print(json.dumps(results, ensure_ascii=False, indent=2))
        
        with open(args.output, "w", encoding="utf-8") as json_file:
            json.dump(results, json_file, ensure_ascii=False, indent=4)
        
        print(f"\n▶ Results successfully saved to: {args.output}")
        
    except Exception as e:
        print(f"\n[Fatal Error] {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()