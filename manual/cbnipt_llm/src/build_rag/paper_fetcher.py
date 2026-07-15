import sys
import time
import json
import argparse
import xml.etree.ElementTree as ET
from typing import List, Dict, Any
from Bio import Entrez

# ---------------------------------------------------------
# 1. ClinVar 파싱 로직
# ---------------------------------------------------------

def get_review_status_stars(status_text: str) -> int:
    # [Rule 0] Default 배제 및 엄격한 예외 처리
    if not status_text:
        raise ValueError("status_text 파라미터가 엄격히 요구됩니다.")
    
    mapping = {
        "practice guideline": 4,
        "reviewed by expert panel": 3,
        "criteria provided, multiple submitters, no conflicts": 2,
        "criteria provided, single submitter": 1,
        "criteria provided, conflicting classifications": 1,
        "no assertion criteria provided": 0,
        "no assertion provided": 0,
        "no classifications from unflagged records": 0
    }
    return mapping.get(status_text.lower(), 0)

def extract_clinvar_evidence(email: str, variation_id: str) -> Dict[str, Any]:
    # [Rule 0]
    if not email: raise ValueError("email 파라미터가 누락되었습니다.")
    if not variation_id: raise ValueError("variation_id 파라미터가 누락되었습니다.")

    Entrez.email = email
    time.sleep(0.35)

    try:
        with Entrez.efetch(db="clinvar", id=variation_id, rettype="vcv", retmode="xml") as handle:
            raw_xml = handle.read()
            
        xml_text = raw_xml.decode("utf-8") if isinstance(raw_xml, bytes) else raw_xml
        
        if not xml_text.strip():
            raise RuntimeError("NCBI로부터 빈 응답이 반환되었습니다.")
            
        root = ET.fromstring(xml_text)
    except Exception as e:
        raise RuntimeError(f"ClinVar VCV XML 다운로드 및 파싱 실패 (Variation ID: {variation_id}): {e}")

    # Review Status 파싱
    rcv_list = root.find(".//RCVList")
    review_status = "no assertion provided"
    if rcv_list is not None:
        rcv_accession = rcv_list.find(".//RCVAccession")
        if rcv_accession is not None:
            review_status = rcv_accession.get("ReviewStatus")
            if not review_status:
                raise ValueError("RCVAccession 태그는 존재하나 ReviewStatus 속성이 없습니다.")

    stars = get_review_status_stars(status_text=review_status)

    # Citation (PMID) 파싱
    pmids = set()
    for citation in root.findall(".//Citation"):
        for id_node in citation.findall(".//ID"):
            print(citation)
            if id_node.get("Source") == "PubMed" and id_node.text:
                pmids.add(id_node.text.strip())

    return {
        "variation_id": variation_id,
        "review_status": review_status,
        "star_rating": stars,
        "pmids": list(pmids)
    }

# ---------------------------------------------------------
# 2. PubMed 파싱 로직
# ---------------------------------------------------------

def fetch_pubmed_details(email: str, pmids: List[str]) -> List[Dict[str, str]]:
    # [Rule 0]
    if not email: raise ValueError("email 파라미터가 누락되었습니다.")
    if not pmids or len(pmids) == 0: raise ValueError("pmids 리스트가 비어있습니다.")

    Entrez.email = email
    time.sleep(0.35)

    try:
        with Entrez.esummary(db="pubmed", id=",".join(pmids)) as handle:
            summaries = Entrez.read(handle)
    except Exception as e:
        raise RuntimeError(f"PubMed esummary API 호출 실패: {e}")

    papers = []
    for summary in summaries:
        pmid = str(summary.get("Id"))
        title = summary.get("Title")
        journal = summary.get("Source")
        
        if not pmid or not title:
            raise ValueError(f"PubMed 메타데이터 누락 (PMID: {pmid}, Title: {title})")
            
        papers.append({
            "pmid": pmid,
            "title": title,
            "journal": journal,
            "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
        })

    return papers

def search_disease_papers(email: str, disease_name: str, max_results: int) -> List[Dict[str, str]]:
    # [Rule 0]
    if not email: raise ValueError("email 파라미터가 누락되었습니다.")
    if not disease_name: raise ValueError("disease_name 파라미터가 누락되었습니다.")
    if not isinstance(max_results, int) or max_results <= 0:
        raise ValueError("max_results는 1 이상의 정수여야 합니다.")

    Entrez.email = email
    time.sleep(0.35)

    try:
        # 논문 제목/초록에 질병명이 포함되며, Publication Type이 Review인 논문 검색
        query = f'"{disease_name}"[Title/Abstract] AND Review[Publication Type]'
        with Entrez.esearch(db="pubmed", term=query, retmax=max_results) as handle:
            record = Entrez.read(handle)
        pmids = record.get("IdList", [])
    except Exception as e:
        raise RuntimeError(f"PubMed 질병 논문 검색 실패: {e}")

    if not pmids:
        return []

    return fetch_pubmed_details(email=email, pmids=pmids)

# ---------------------------------------------------------
# 3. 통합 실행 파이프라인
# ---------------------------------------------------------

def run_evidence_pipeline(email: str, variation_id: str, disease_name: str, max_disease_papers: int) -> Dict[str, Any]:
    # [Rule 0] 필수 파라미터 검증
    if not email: raise ValueError("email은 필수입니다.")
    if not variation_id: raise ValueError("variation_id는 필수입니다.")
    if not disease_name: raise ValueError("disease_name은 필수입니다.")
    if not max_disease_papers: raise ValueError("max_disease_papers는 필수입니다.")

    print(f"[*] ClinVar에서 변이 {variation_id}의 증거 메타데이터를 추출합니다...")
    clinvar_data = extract_clinvar_evidence(email=email, variation_id=variation_id)
    
    evidence_papers = []
    if clinvar_data["pmids"]:
        print(f"[*] 추출된 {len(clinvar_data['pmids'])}개의 증거 논문(PMID) 메타데이터를 가져옵니다...")
        evidence_papers = fetch_pubmed_details(email=email, pmids=clinvar_data["pmids"])
    print(evidence_papers)
    disease_papers = []
    print(f"[*] '{disease_name}'에 대한 리뷰 논문을 최대 {max_disease_papers}개 검색합니다...")
    disease_papers = search_disease_papers(email=email, disease_name=disease_name, max_results=max_disease_papers)

    return {
        "variation_id": clinvar_data["variation_id"],
        "review_status": clinvar_data["review_status"],
        "star_rating": clinvar_data["star_rating"],
        "evidence_papers": evidence_papers,
        "disease_review_papers": disease_papers
    }

# ---------------------------------------------------------
# CLI 실행부
# ---------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="ClinVar Variation 증거 및 질병 리뷰 논문 단일 추출 스크립트")
    
    # [Rule 0] default 미설정, required=True 강제
    parser.add_argument("--email", required=True, help="NCBI API 인증용 이메일")
    parser.add_argument("--vid", required=True, help="ClinVar Variation ID (예: 16209)")
    parser.add_argument("--disease", required=True, help="백그라운드 검색용 질병명 (예: Achondroplasia)")
    parser.add_argument("--max-papers", required=True, type=int, help="수집할 질병 리뷰 논문의 최대 개수")

    args = parser.parse_args()

    try:
        result = run_evidence_pipeline(
            email=args.email,
            variation_id=args.vid,
            disease_name=args.disease,
            max_disease_papers=args.max_papers
        )
        
        print("\n=======================================================")
        print(f"   [결과] Variation ID: {args.vid} / Disease: {args.disease}   ")
        print("=======================================================")
        print(json.dumps(result, ensure_ascii=False, indent=2))
        
    except Exception as e:
        print(f"\n[Fatal Error] 파이프라인 실행 중 치명적 오류 발생: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()