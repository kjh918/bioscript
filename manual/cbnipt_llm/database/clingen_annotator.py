"""
syndrome_annotator.py
─────────────────────
Syndrome 이름을 input으로 받아 아래 DB에서 데이터를 통합 수집하여
JSON + TSV 형태로 저장합니다.

수집 소스:
  1. 마커셋 TSV      → 해당 syndrome의 feature 목록 (CoreGene, CoreRegion 등)
  2. MedGen          → 질환 정의, synonyms, cross-refs, 카테고리별 논문
  3. PubMed          → abstract + PDF 링크
  4. NCBI Gene       → CoreGene별 gene_id, full_name, summary, map_location
  5. UniProt         → 단백질 기능, 도메인, 질환 연관
  6. ClinGen Dosage  → HI/TS score (FTP 다운로드 또는 API)
  7. GenCC           → 유전자-질환 분류 강도, MOI, 제출기관
  8. HPO             → 유전자-표현형 연관 HPO term

출력:
  {syndrome_id}_annotated.json   ← 전체 구조 (중첩 JSON)
  {syndrome_id}_genes.tsv        ← CoreGene 행 flat TSV
  {syndrome_id}_regions.tsv      ← CoreRegion/PartialChromosome 행 flat TSV

Usage:
  python syndrome_annotator.py \\
      --syndrome "DiGeorge syndrome" \\
      --markerset markerset.tsv \\
      --email your@email.com \\
      --output-dir ./output

  python syndrome_annotator.py \\
      --nipt-id NIPT_DGS \\
      --markerset markerset.tsv \\
      --email your@email.com

사전 다운로드 필요 (대용량 파일):
  ClinGen Dosage:
    wget https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv
    wget https://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh38.tsv
  GenCC:
    wget https://search.thegencc.org/download/action/submissions-export-tsv -O gencc_submissions.tsv
  HPO:
    wget https://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt
"""

import os
import re
import csv
import sys
import json
import time
import gzip
import argparse
import urllib.parse
from io import StringIO
from pathlib import Path
from typing import Any

import requests
from requests.adapters import HTTPAdapter
from clingen_parser import parse_clingen_gene, parse_clingen_region, get_region_overlaps
from urllib3.util.retry import Retry

# ── 상수 ──────────────────────────────────────────────────────
EUTILS   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UNIPROT  = "https://rest.uniprot.org/uniprotkb"
PUBMED   = "https://pubmed.ncbi.nlm.nih.gov"

HI_SCORE_LABEL = {
    "0":  "No Evidence",
    "1":  "Little Evidence",
    "2":  "Emerging Evidence",
    "3":  "Sufficient Evidence",
    "30": "Associated with Autosomal Recessive",
    "40": "Dosage Sensitivity Unlikely",
}

# ── HTTP 세션 ─────────────────────────────────────────────────

def _session(email: str) -> requests.Session:
    s = requests.Session()
    retry = Retry(total=5, backoff_factor=1.5,
                  status_forcelist=[429, 500, 502, 503, 504],
                  allowed_methods=["GET"])
    s.mount("https://", HTTPAdapter(max_retries=retry))
    s.mount("http://",  HTTPAdapter(max_retries=retry))
    s.headers.update({"User-Agent": f"syndrome-annotator/1.0 (mailto:{email})"})
    return s

SESSION: requests.Session | None = None

def _get(url: str, params: dict | None = None, delay: float = 0.35) -> requests.Response:
    time.sleep(delay)
    r = SESSION.get(url, params=params, timeout=30)
    r.raise_for_status()
    return r


# ════════════════════════════════════════════════════════════════
# 1. 마커셋 TSV 파서
# ════════════════════════════════════════════════════════════════

def load_markerset(tsv_path: str) -> list[dict]:
    rows = []
    with open(tsv_path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append({k.strip(): v.strip() for k, v in row.items()})
    return rows


def find_syndrome_rows(rows: list[dict],
                       syndrome: str | None,
                       nipt_id: str | None) -> list[dict]:
    results = []
    for row in rows:
        if nipt_id and row.get("NIPT_ID", "").upper() == nipt_id.upper():
            results.append(row)
        elif syndrome and syndrome.lower() in row.get("SYNDROME", "").lower():
            results.append(row)
    return results


def parse_markerset_features(rows: list[dict]) -> dict:
    """
    마커셋 행들을 feature_type별로 분류.
    반환:
      {
        nipt_id, syndrome, nipt_group,
        target_chromosomes: [...],
        partial_chromosomes: [...],
        core_regions: [...],
        core_genes: [{ feature_name, chrom, start, end, size_mb }]
      }
    """
    if not rows:
        return {}

    first = rows[0]
    result = {
        "nipt_id":           first.get("NIPT_ID", ""),
        "syndrome":          first.get("SYNDROME", ""),
        "nipt_group":        first.get("NIPT_GROUP", ""),
        "target_chromosomes": [],
        "partial_chromosomes": [],
        "core_regions":      [],
        "core_genes":        [],
    }

    for row in rows:
        ftype = row.get("FEATURE_TYPE", "")
        entry = {
            "feature_name":     row.get("FEATURE_NAME", ""),
            "chromosome":       row.get("CHROMOSOME", ""),
            "start":            int(row.get("GENOMIC_POS_START", 0) or 0),
            "end":              int(row.get("GENOMIC_POS_END", 0) or 0),
            "size_mb":          float(row.get("SIZE_Mb", 0) or 0),
        }
        if ftype == "TargetChromosome":
            result["target_chromosomes"].append(entry)
        elif ftype == "PartialChromosome":
            result["partial_chromosomes"].append(entry)
        elif ftype in ("CoreRegion", "PrimaryTargetRegion"):
            entry["feature_type"] = ftype   # 원래 타입 보존
            result["core_regions"].append(entry)
        elif ftype == "CoreGene":
            result["core_genes"].append(entry)

    return result


# ════════════════════════════════════════════════════════════════
# 2. MedGen 파서 (질환 정의 + 논문)
# ════════════════════════════════════════════════════════════════

def fetch_medgen(syndrome_name: str, email: str) -> dict:
    print(f"  [MedGen] '{syndrome_name}' 검색 중...")
    try:
        r = _get(f"{EUTILS}/esearch.fcgi",
                 params={"db": "medgen", "term": f'"{syndrome_name}"',
                         "retmax": 1, "sort": "relevance",
                         "retmode": "json", "email": email})
        id_list = r.json().get("esearchresult", {}).get("idlist", [])
        if not id_list:
            print(f"  [MedGen] 결과 없음")
            return {}

        uid = id_list[0]
        url = f"https://www.ncbi.nlm.nih.gov/medgen/{uid}"

        time.sleep(1.5)
        resp = SESSION.get(url, timeout=20,
                           headers={"User-Agent": "Mozilla/5.0 (Macintosh) Chrome/124"})
        resp.raise_for_status()

        from bs4 import BeautifulSoup
        soup = BeautifulSoup(resp.text, "lxml")

        # definition
        definition = ""
        meta = soup.find("meta", attrs={"name": "description"})
        if meta and meta.get("content"):
            definition = meta["content"].strip()
        portlet = soup.find("div", id="ID_100")
        if portlet:
            content = portlet.find("div", class_="portlet_content")
            if content:
                definition = re.sub(r"\s+", " ",
                                    content.get_text(separator=" ", strip=True))

        # synonyms
        synonyms = []
        table = soup.find("table", class_="medgenTable")
        if table:
            for tr in table.find_all("tr"):
                tds = tr.find_all("td")
                if len(tds) >= 2 and "synonym" in tds[0].get_text().lower():
                    raw = tds[1].get_text(separator=";", strip=True)
                    synonyms = [s.strip() for s in raw.split(";") if s.strip()]

        # cross-refs
        xrefs: dict[str, str] = {}
        if table:
            for tr in table.find_all("tr"):
                tds = tr.find_all("td")
                if len(tds) < 2:
                    continue
                label = tds[0].get_text().strip().rstrip(":").lower()
                a = tds[1].find("a")
                val = (a.get_text() if a else tds[1].get_text()).strip()
                if "monarch" in label or "mondo" in label:
                    xrefs["MONDO"] = val
                elif "omim" in label:
                    xrefs["OMIM"] = val
                elif "orphanet" in label:
                    xrefs["Orphanet"] = re.sub(r"^ORPHA", "", val)
                elif "snomed" in label:
                    m = re.search(r"\((\d+)\)", val)
                    xrefs["SNOMED-CT"] = m.group(1) if m else val

        # 논문 (카테고리별)
        literature = _parse_medgen_literature(soup)

        print(f"  [MedGen] UID={uid} | def={bool(definition)} | "
              f"xrefs={list(xrefs.keys())} | "
              f"lit_cats={list(literature.keys())}")

        return {
            "medgen_uid":       uid,
            "medgen_url":       url,
            "definition":       definition,
            "synonyms":         synonyms,
            "cross_references": xrefs,
            "literature":       literature,
        }

    except Exception as e:
        print(f"  [MedGen] 오류: {e}")
        return {}


def _parse_medgen_literature(soup) -> dict:
    """MedGen HTML에서 카테고리별 논문 파싱."""
    SUBHEAD_MAP = {
        "etiology":              "etiology",
        "diagnosis":             "diagnosis",
        "therapy":               "therapy",
        "prognosis":             "prognosis",
        "clinical prediction":   "clinical_prediction_guides",
    }
    literature: dict[str, dict] = {}

    # Professional guidelines (ID_105)
    portlet105 = soup.find("div", id="ID_105")
    if portlet105:
        content = portlet105.find("div", class_="portlet_content")
        if content:
            papers, see_all = _extract_papers_from_section(content)
            literature["professional_guidelines"] = {
                "papers": papers, "see_all": see_all
            }

    # Recent clinical studies (ID_103)
    portlet103 = soup.find("div", id="ID_103")
    if portlet103:
        content = portlet103.find("div", class_="portlet_content")
        if content:
            current_cat = "_default"
            for child in list(content.children):
                from bs4 import Tag
                if not isinstance(child, Tag):
                    continue
                if child.name == "h3" and "subhead" in child.get("class", []):
                    h3t = child.get_text(strip=True).lower()
                    current_cat = next(
                        (v for k, v in SUBHEAD_MAP.items() if k in h3t), h3t
                    )
                    if current_cat not in literature:
                        literature[current_cat] = {"papers": [], "see_all": {}}

            # 재파싱: 헤딩 기준으로 논문 분리
            current_cat = "_default"
            children = list(content.children)
            i = 0
            while i < len(children):
                node = children[i]
                from bs4 import Tag
                if not isinstance(node, Tag):
                    i += 1; continue

                if node.name == "h3" and "subhead" in node.get("class", []):
                    h3t = node.get_text(strip=True).lower()
                    current_cat = next(
                        (v for k, v in SUBHEAD_MAP.items() if k in h3t), h3t
                    )
                    if current_cat not in literature:
                        literature[current_cat] = {"papers": [], "see_all": {}}
                    i += 1; continue

                if node.name == "div" and "nl" in node.get("class", []):
                    # See all 체크
                    see_a = node.find("a", attrs={"data-ga-label":
                                                   re.compile(r"See all", re.I)})
                    if see_a:
                        label = see_a.get("data-ga-label", "")
                        m = re.search(r"See all \((\d+)\)", label)
                        if m and current_cat in literature:
                            href = see_a.get("href", "")
                            if href.startswith("/"):
                                href = "https://www.ncbi.nlm.nih.gov" + href
                            literature[current_cat]["see_all"] = {
                                "count": int(m.group(1)), "url": href
                            }
                        i += 1; continue

                    # 논문 제목 a 태그
                    a = node.find("a")
                    if a:
                        title = re.sub(r"\s+", " ", a.get_text(strip=True))
                        pmid = ""
                        m = re.search(r"/pubmed/(\d+)", a.get("href", ""))
                        if m:
                            pmid = m.group(1)

                        # 다음 detail div
                        detail = None
                        j = i + 1
                        while j < len(children):
                            nxt = children[j]
                            if not isinstance(nxt, Tag):
                                j += 1; continue
                            if ("portlet_content" in nxt.get("class", []) and
                                    "ln" in nxt.get("class", [])):
                                detail = nxt
                                i = j
                            break

                        authors, journal, year = [], "", ""
                        if detail:
                            asp = detail.find("span", class_="medgenPMauthor")
                            if asp:
                                authors = [a2.strip() for a2 in
                                           asp.get_text(strip=True).split(",")
                                           if a2.strip()]
                            jsp = detail.find("span", class_="medgenPMjournal")
                            if jsp:
                                journal = jsp.get_text(strip=True)
                            ym = re.search(r"\b(19|20)\d{2}\b",
                                           detail.get_text())
                            if ym:
                                year = ym.group(0)
                            if not pmid:
                                pa = detail.find("a", attrs={
                                    "data-ga-action": "PMID"})
                                if pa:
                                    pmid = pa.get_text(strip=True)

                        pmc_free = bool(detail and
                                        detail.find("a", class_="PubMedFree"))

                        paper = {
                            "pmid": pmid, "title": title,
                            "authors": authors, "journal": journal,
                            "year": year, "pmc_free": pmc_free,
                            "url": (f"{PUBMED}/{pmid}/" if pmid else ""),
                        }
                        if current_cat in literature:
                            literature[current_cat]["papers"].append(paper)
                i += 1

    # Systematic reviews (ID_104)
    portlet104 = soup.find("div", id="ID_104")
    if portlet104:
        content = portlet104.find("div", class_="portlet_content")
        if content:
            papers, see_all = _extract_papers_from_section(content)
            literature["systematic_reviews"] = {
                "papers": papers, "see_all": see_all
            }

    return literature


def _extract_papers_from_section(content) -> tuple[list, dict]:
    """단일 섹션에서 논문 목록 + see_all 추출."""
    from bs4 import Tag
    papers, see_all = [], {}
    children = list(content.children)
    i = 0
    while i < len(children):
        node = children[i]
        if not isinstance(node, Tag):
            i += 1; continue

        if node.name == "div" and "nl" in node.get("class", []):
            see_a = node.find("a", attrs={"data-ga-label":
                                           re.compile(r"See all", re.I)})
            if see_a:
                label = see_a.get("data-ga-label", "")
                m = re.search(r"See all \((\d+)\)", label)
                if m:
                    href = see_a.get("href", "")
                    if href.startswith("/"):
                        href = "https://www.ncbi.nlm.nih.gov" + href
                    see_all = {"count": int(m.group(1)), "url": href}
                i += 1; continue

            a = node.find("a")
            if a:
                title = re.sub(r"\s+", " ", a.get_text(strip=True))
                pmid = ""
                m = re.search(r"/pubmed/(\d+)", a.get("href", ""))
                if m:
                    pmid = m.group(1)
                papers.append({"pmid": pmid, "title": title,
                                "url": f"{PUBMED}/{pmid}/" if pmid else ""})
        i += 1
    return papers, see_all


# ════════════════════════════════════════════════════════════════
# 3. PubMed abstract 수집
# ════════════════════════════════════════════════════════════════

def fetch_pubmed_abstracts(pmids: list[str], email: str,
                           batch_size: int = 100) -> dict[str, dict]:
    """PMID 목록 → {pmid: {title, abstract, abstract_sections, authors, ...}}"""
    if not pmids:
        return {}

    import xml.etree.ElementTree as ET
    lookup: dict[str, dict] = {}
    pmids = [p for p in pmids if p]

    print(f"  [PubMed] {len(pmids)}개 abstract 수집 중...")
    for i in range(0, len(pmids), batch_size):
        batch = pmids[i:i + batch_size]
        try:
            r = _get(f"{EUTILS}/efetch.fcgi",
                     params={"db": "pubmed", "id": ",".join(batch),
                             "rettype": "abstract", "retmode": "xml",
                             "email": email})
            root = ET.fromstring(r.text.encode("utf-8"))
        except Exception as e:
            print(f"  [PubMed] 배치 오류: {e}")
            continue

        for pa in root.findall(".//PubmedArticle"):
            mc = pa.find("MedlineCitation")
            if mc is None:
                continue
            pmid_el = mc.find("PMID")
            pmid = pmid_el.text.strip() if pmid_el is not None else ""
            art = mc.find("Article")
            if art is None:
                continue

            title_el = art.find("ArticleTitle")
            title = "".join(title_el.itertext()).strip() if title_el is not None else ""

            # Abstract
            abstract_el = art.find(".//Abstract")
            abstract, sections = "", {}
            if abstract_el:
                parts = []
                for at in abstract_el.findall("AbstractText"):
                    label = at.get("Label") or at.get("NlmCategory") or ""
                    text = re.sub(r"\s+", " ", "".join(at.itertext()).strip())
                    if label:
                        sections[label] = text
                        parts.append(f"{label}: {text}")
                    else:
                        parts.append(text)
                abstract = " ".join(parts)

            # Authors
            authors = []
            for au in art.findall(".//AuthorList/Author"):
                last = (au.findtext("LastName") or "").strip()
                fore = (au.findtext("ForeName") or
                        au.findtext("Initials") or "").strip()
                cname = (au.findtext("CollectiveName") or "").strip()
                if cname:
                    authors.append(cname)
                elif last:
                    authors.append(f"{last} {fore}".strip())

            # Journal + year
            journal = (art.findtext(".//Journal/Title") or
                       art.findtext(".//Journal/ISOAbbreviation") or "")
            pub_date = art.find(".//PubDate")
            year = ""
            if pub_date is not None:
                year = pub_date.findtext("Year") or ""
                if not year:
                    med = pub_date.findtext("MedlineDate") or ""
                    m = re.match(r"(\d{4})", med)
                    year = m.group(1) if m else ""

            # DOI + PMC
            doi, pmc_id = "", ""
            pd = pa.find("PubmedData")
            if pd is not None:
                for aid in pd.findall(".//ArticleIdList/ArticleId"):
                    t = aid.get("IdType", "")
                    v = (aid.text or "").strip()
                    if t == "doi":
                        doi = v
                    elif t == "pmc":
                        pmc_id = v

            # PDF links
            pdf_links = []
            if pmc_id:
                pmc_num = pmc_id.replace("PMC", "")
                pdf_links.append({
                    "source": "PMC", "type": "pdf",
                    "url": f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmc_id}/pdf/",
                    "note": "PMC free full text PDF",
                })
                pdf_links.append({
                    "source": "PMC", "type": "html",
                    "url": f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmc_id}/",
                    "note": "PMC full text HTML",
                })
            if doi:
                pdf_links.append({
                    "source": "DOI", "type": "publisher",
                    "url": f"https://doi.org/{doi}",
                    "note": "Publisher page",
                })
            pdf_links.append({
                "source": "PubMed", "type": "abstract",
                "url": f"{PUBMED}/{pmid}/",
                "note": "PubMed abstract",
            })

            # MeSH
            mesh = [mh.findtext("DescriptorName") or ""
                    for mh in mc.findall(".//MeshHeadingList/MeshHeading")]

            lookup[pmid] = {
                "pmid": pmid, "doi": doi, "pmc_id": pmc_id,
                "title": title, "authors": authors,
                "journal": journal, "year": year,
                "abstract": abstract,
                "abstract_sections": sections,
                "mesh_terms": [m for m in mesh if m],
                "full_text_available": bool(pmc_id),
                "pdf_links": pdf_links,
                "pubmed_url": f"{PUBMED}/{pmid}/",
            }

    print(f"  [PubMed] {len(lookup)}개 abstract 수집 완료")
    return lookup


# ════════════════════════════════════════════════════════════════
# 4. NCBI Gene
# ════════════════════════════════════════════════════════════════

def fetch_ncbi_gene(gene_symbol: str, email: str,
                    hgnc_id: str = "", cytoband: str = "") -> dict:
    """
    NCBI Gene에서 Human(Homo sapiens, taxid=9606) 유전자만 정확히 조회.

    [방어 전략 3단계]
    1순위: HGNC ID → Gene Accession 직접 검색 (마우스 ortholog 혼입 방지)
    2순위: Gene Name + Homo sapiens + 염색체 번호 필터
    3순위: Gene Name + Homo sapiens → taxid=9606 검증
    """
    def _esummary(gene_id: str) -> dict:
        time.sleep(0.35)
        r = _get(f"{EUTILS}/esummary.fcgi",
                 params={"db": "gene", "id": gene_id,
                         "retmode": "json", "email": email})
        return r.json().get("result", {}).get(gene_id, {})

    def _is_human(data: dict) -> bool:
        org = data.get("organism", {})
        return (org.get("taxid", 0) == 9606 or
                "homo sapiens" in org.get("scientificname", "").lower())

    def _build(gene_id: str, data: dict) -> dict:
        ginfo = (data.get("genomicinfo", [{}])[0]
                 if data.get("genomicinfo") else {})
        return {
            "ncbi_gene_id":    gene_id,
            "gene_full_name":  data.get("description", ""),
            "gene_summary":    data.get("summary", ""),
            "hgnc_id":         data.get("nomenclatureauthorityid", ""),
            "map_location":    data.get("maplocation", ""),
            "chromosome":      data.get("chromosome", ""),
            "chr_start":       ginfo.get("chrstart", ""),
            "chr_stop":        ginfo.get("chrstop", ""),
            "chr_accver":      ginfo.get("chraccver", ""),
            "ensembl_gene_id": next(
                (xr.get("value", "") for xr in data.get("locationhist", [])
                 if "ENSG" in xr.get("value", "")), ""
            ),
            "ncbi_gene_url":   f"https://www.ncbi.nlm.nih.gov/gene/{gene_id}",
        }

    try:
        # 1순위: HGNC ID 직접 검색
        hgnc_num = re.sub(r"^HGNC:", "", hgnc_id).strip() if hgnc_id else ""
        if hgnc_num:
            r = _get(f"{EUTILS}/esearch.fcgi",
                     params={"db": "gene",
                             "term": f"HGNC:{hgnc_num}[Gene Accession]",
                             "retmax": 1, "retmode": "json", "email": email})
            ids = r.json().get("esearchresult", {}).get("idlist", [])
            if ids:
                data = _esummary(ids[0])
                if _is_human(data):
                    print(f"  [NCBI Gene] {gene_symbol} → HGNC:{hgnc_num} "
                          f"gene_id={ids[0]} chr={data.get('chromosome','')} "
                          f"loc={data.get('maplocation','')}")
                    return _build(ids[0], data)

        # 2순위: Gene Name + Organism + CHR
        chrom_num = ""
        if cytoband:
            m = re.match(r"^(\d+|X|Y)", str(cytoband).replace("chr", ""))
            if m:
                chrom_num = m.group(1)

        if chrom_num:
            query = (f"{gene_symbol}[Gene Name] AND "
                     f"Homo sapiens[Organism] AND {chrom_num}[CHR]")
            r = _get(f"{EUTILS}/esearch.fcgi",
                     params={"db": "gene", "term": query,
                             "retmax": 5, "sort": "relevance",
                             "retmode": "json", "email": email})
            for gene_id in r.json().get("esearchresult", {}).get("idlist", []):
                data = _esummary(gene_id)
                if _is_human(data):
                    print(f"  [NCBI Gene] {gene_symbol} → chr{chrom_num} filter "
                          f"gene_id={gene_id} loc={data.get('maplocation','')}")
                    return _build(gene_id, data)

        # 3순위: Gene Name + Organism → taxid 검증
        r = _get(f"{EUTILS}/esearch.fcgi",
                 params={"db": "gene",
                         "term": f"{gene_symbol}[Gene Name] AND Homo sapiens[Organism]",
                         "retmax": 5, "sort": "relevance",
                         "retmode": "json", "email": email})
        for gene_id in r.json().get("esearchresult", {}).get("idlist", []):
            data = _esummary(gene_id)
            if _is_human(data):
                print(f"  [NCBI Gene] {gene_symbol} → organism filter "
                      f"gene_id={gene_id} loc={data.get('maplocation','')}")
                return _build(gene_id, data)

        print(f"  [NCBI Gene] {gene_symbol} → Human 결과 없음")
        return {}

    except Exception as e:
        print(f"  [NCBI Gene] {gene_symbol} 오류: {e}")
        return {}


# ════════════════════════════════════════════════════════════════
# 5. UniProt
# ════════════════════════════════════════════════════════════════

def fetch_uniprot(gene_symbol: str) -> dict:
    """유전자 심볼 → UniProt reviewed entry (Swiss-Prot 우선)"""
    try:
        query = (f"gene_exact:{gene_symbol} AND organism_id:9606 AND reviewed:true")
        r = _get(UNIPROT,
                 params={"query": query, "format": "json", "size": 1})
        results = r.json().get("results", [])
        if not results:
            return {}

        entry = results[0]
        accession = entry.get("primaryAccession", "")

        # 기능 설명
        function = ""
        for comment in entry.get("comments", []):
            if comment.get("commentType") == "FUNCTION":
                texts = comment.get("texts", [])
                if texts:
                    function = texts[0].get("value", "")
                    break

        # 도메인 목록
        domains = [
            {
                "name":  f.get("description", ""),
                "start": f.get("location", {}).get("start", {}).get("value", ""),
                "end":   f.get("location", {}).get("end", {}).get("value", ""),
            }
            for f in entry.get("features", [])
            if f.get("type") == "Domain"
        ]

        # 질환 연관
        diseases = []
        for comment in entry.get("comments", []):
            if comment.get("commentType") == "DISEASE":
                d = comment.get("disease", {})
                diseases.append({
                    "name":        d.get("diseaseId", ""),
                    "description": next(
                        (t.get("value", "") for t in comment.get("texts", [])), ""
                    ),
                    "omim":        d.get("diseaseCrossReference", {}).get("id", ""),
                })

        # 키워드
        keywords = [kw.get("name", "") for kw in entry.get("keywords", [])]

        # Subcellular location
        subcell = []
        for comment in entry.get("comments", []):
            if comment.get("commentType") == "SUBCELLULAR LOCATION":
                for loc in comment.get("subcellularLocations", []):
                    loc_name = loc.get("location", {}).get("value", "")
                    if loc_name:
                        subcell.append(loc_name)

        return {
            "uniprot_accession":   accession,
            "uniprot_id":          entry.get("uniProtkbId", ""),
            "protein_name":        (entry.get("proteinDescription", {})
                                    .get("recommendedName", {})
                                    .get("fullName", {}).get("value", "")),
            "protein_function":    function,
            "protein_domains":     domains,
            "disease_associations":diseases,
            "subcellular_location":subcell,
            "keywords":            keywords[:10],
            "uniprot_url":         f"https://www.uniprot.org/uniprotkb/{accession}",
        }
    except Exception as e:
        print(f"  [UniProt] {gene_symbol} 오류: {e}")
        return {}


# ════════════════════════════════════════════════════════════════
# 6. ClinGen Dosage Sensitivity (로컬 파일)
# ════════════════════════════════════════════════════════════════

_clingen_gene_cache: dict[str, dict] | None = None
_clingen_region_cache: list[dict] | None = None

def load_clingen_gene(path: str) -> dict[str, dict]:
    """ClinGen_gene_curation_list_GRCh38.tsv → {gene_symbol: {...}}"""
    global _clingen_gene_cache
    if _clingen_gene_cache is not None:
        return _clingen_gene_cache

    lookup: dict[str, dict] = {}
    if not path or not Path(path).exists():
        print(f"  [ClinGen] gene 파일 없음: {path}")
        return lookup

    with open(path, encoding="utf-8") as f:
        # 헤더 라인이 #으로 시작하는 경우 스킵
        lines = [l for l in f if not l.startswith("#")]

    reader = csv.DictReader(StringIO("".join(lines)), delimiter="\t")
    for row in reader:
        sym = row.get("Gene Symbol", "").strip()
        if sym:
            lookup[sym] = {
                "ncbi_gene_id_clingen": row.get("Gene ID", "").strip(),
                "cytoband":             row.get("cytoBand", "").strip(),
                "hi_score":             row.get("Haploinsufficiency Score", "").strip(),
                "hi_score_label":       HI_SCORE_LABEL.get(
                    row.get("Haploinsufficiency Score", "").strip(), ""),
                "hi_disease_id":        row.get("Haploinsufficiency Disease ID", "").strip(),
                "hi_disease_name":      row.get("Haploinsufficiency Disease Name", "").strip(),
                "ts_score":             row.get("Triplosensitivity Score", "").strip(),
                "ts_score_label":       HI_SCORE_LABEL.get(
                    row.get("Triplosensitivity Score", "").strip(), ""),
                "ts_disease_id":        row.get("Triplosensitivity Disease ID", "").strip(),
                "ts_disease_name":      row.get("Triplosensitivity Disease Name", "").strip(),
            }
    _clingen_gene_cache = lookup
    print(f"  [ClinGen] gene {len(lookup)}개 로드")
    return lookup


def load_clingen_region(path: str) -> list[dict]:
    """ClinGen_region_curation_list_GRCh38.tsv → 목록"""
    global _clingen_region_cache
    if _clingen_region_cache is not None:
        return _clingen_region_cache

    regions = []
    if not path or not Path(path).exists():
        print(f"  [ClinGen] region 파일 없음: {path}")
        return regions

    with open(path, encoding="utf-8") as f:
        lines = [l for l in f if not l.startswith("#")]

    reader = csv.DictReader(StringIO("".join(lines)), delimiter="\t")
    
    for row in reader:
        print(row)
        exit()
        loc = row.get("Genomic Location GRCh38", "")
        m = re.match(r"chr([\dXY]+):(\d+)-(\d+)", loc)
        if m:
            regions.append({
                "chrom":          "chr" + m.group(1),
                "start":          int(m.group(2)),
                "end":            int(m.group(3)),
                "region_name":    row.get("ISCA Region Name", "").strip(),
                "hi_score":       row.get("Haploinsufficiency Score", "").strip(),
                "hi_score_label": HI_SCORE_LABEL.get(
                    row.get("Haploinsufficiency Score", "").strip(), ""),
                "hi_disease_id":  row.get("Haploinsufficiency Disease ID", "").strip(),
                "ts_score":       row.get("Triplosensitivity Score", "").strip(),
                "ts_score_label": HI_SCORE_LABEL.get(
                    row.get("Triplosensitivity Score", "").strip(), ""),
            })
    _clingen_region_cache = regions
    print(f"  [ClinGen] region {len(regions)}개 로드")
    return regions


def get_region_overlaps(chrom: str, start: int, end: int,
                                regions: list[dict]) -> list[dict]:
    """좌표 overlap으로 ClinGen region 매칭"""
    return [
        r for r in regions
        if r["chrom"] == chrom and r["start"] <= end and r["end"] >= start
    ]


# ════════════════════════════════════════════════════════════════
# 7. GenCC (로컬 TSV)
# ════════════════════════════════════════════════════════════════

_gencc_cache: dict[str, list[dict]] | None = None

def load_gencc(path: str) -> dict[str, list[dict]]:
    """GenCC submissions TSV → {gene_symbol: [entries]}"""
    global _gencc_cache
    if _gencc_cache is not None:
        return _gencc_cache

    lookup: dict[str, list[dict]] = {}
    if not path or not Path(path).exists():
        print(f"  [GenCC] 파일 없음: {path}")
        return lookup

    with open(path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sym = row.get("gene_symbol", "").strip()
            if not sym:
                continue
            entry = {
                "disease_curie":       row.get("disease_curie", "").strip(),
                "disease_title":       row.get("disease_title", "").strip(),
                "classification":      row.get("classification_title", "").strip(),
                "moi":                 row.get("moi_title", "").strip(),
                "submitter":           row.get("submitter_title", "").strip(),
                "pmids":               row.get("submitted_as_pmids", "").strip(),
                "report_url":          row.get("submitted_as_public_report_url", "").strip(),
                "submitted_date":      row.get("submitted_as_date", "").strip(),
            }
            lookup.setdefault(sym, []).append(entry)

    _gencc_cache = lookup
    print(f"  [GenCC] {len(lookup)}개 유전자 로드")
    return lookup


# ════════════════════════════════════════════════════════════════
# 8. HPO (로컬 파일)
# ════════════════════════════════════════════════════════════════

_hpo_cache: dict[str, list[dict]] | None = None

def load_hpo(path: str) -> dict[str, list[dict]]:
    """phenotype_to_genes.txt → {gene_symbol: [{hpo_id, hpo_name}]}"""
    global _hpo_cache
    if _hpo_cache is not None:
        return _hpo_cache

    lookup: dict[str, list[dict]] = {}
    if not path or not Path(path).exists():
        print(f"  [HPO] 파일 없음: {path}")
        return lookup

    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            hpo_id, hpo_name = parts[0], parts[1]
            gene_sym = parts[3] if len(parts) > 3 else ""
            if gene_sym:
                lookup.setdefault(gene_sym, []).append(
                    {"hpo_id": hpo_id, "hpo_name": hpo_name}
                )

    _hpo_cache = lookup
    print(f"  [HPO] {len(lookup)}개 유전자-표현형 로드")
    return lookup


# ════════════════════════════════════════════════════════════════
# 통합 수집
# ════════════════════════════════════════════════════════════════

def annotate_syndrome(
    syndrome_name: str | None,
    nipt_id: str | None,
    markerset_path: str,
    email: str,
    clingen_gene_path: str = "",
    clingen_region_path: str = "",
    gencc_path: str = "",
    hpo_path: str = "",
) -> dict[str, Any]:

    # 1. 마커셋 로드
    print(f"\n[1/8] 마커셋 로드: {markerset_path}")
    all_rows = load_markerset(markerset_path)
    rows = find_syndrome_rows(all_rows, syndrome_name, nipt_id)
    if not rows:
        raise ValueError(f"마커셋에서 '{syndrome_name or nipt_id}' 를 찾을 수 없음")

    features = parse_markerset_features(rows)
    actual_syndrome = features["syndrome"]
    actual_nipt_id  = features["nipt_id"]
    core_genes      = [g["feature_name"] for g in features["core_genes"]]
    print(f"  → {actual_syndrome} ({actual_nipt_id}) | CoreGene {len(core_genes)}개: {core_genes}")

    # 2. MedGen
    print(f"\n[2/8] MedGen: '{actual_syndrome}'")
    medgen = fetch_medgen(actual_syndrome, email)

    # 3. PubMed abstracts
    print(f"\n[3/8] PubMed abstracts 수집")
    all_pmids: list[str] = []
    for cat_data in medgen.get("literature", {}).values():
        if isinstance(cat_data, dict):
            for paper in cat_data.get("papers", []):
                if paper.get("pmid"):
                    all_pmids.append(paper["pmid"])
    all_pmids = list(dict.fromkeys(all_pmids))  # 중복 제거 (순서 유지)
    pubmed_lookup = fetch_pubmed_abstracts(all_pmids, email)

    # literature에 abstract 병합
    for cat_key, cat_data in medgen.get("literature", {}).items():
        if isinstance(cat_data, dict):
            for paper in cat_data.get("papers", []):
                pmid = paper.get("pmid", "")
                if pmid and pmid in pubmed_lookup:
                    abs_data = pubmed_lookup[pmid]
                    paper["abstract"]           = abs_data.get("abstract", "")
                    paper["abstract_sections"]  = abs_data.get("abstract_sections", {})
                    paper["mesh_terms"]         = abs_data.get("mesh_terms", [])
                    paper["pdf_links"]          = abs_data.get("pdf_links", [])
                    paper["full_text_available"]= abs_data.get("full_text_available", False)
                    paper["doi"]                = abs_data.get("doi", "")
                    paper["pmc_id"]             = abs_data.get("pmc_id", "")

    # 4~8. CoreGene annotation (외부 DB)
    print(f"\n[4-8/8] CoreGene annotation ({len(core_genes)}개)")
    clingen_gene_db  = parse_clingen_gene(clingen_gene_path)
    clingen_region_db = parse_clingen_region(clingen_region_path)
    gencc_db         = load_gencc(gencc_path)
    hpo_db           = load_hpo(hpo_path)

    for gene_entry in features["core_genes"]:
        sym = gene_entry["feature_name"]
        # ClinGen에서 cytoband 미리 조회 (human 필터용)
        clingen_pre = clingen_gene_db.get(sym, {})
        cytoband_hint = clingen_pre.get("cytoband", "")
        print(f"  → {sym}  cytoband_hint={cytoband_hint}")

        # NCBI Gene (human 검증 포함)
        print(f"    [NCBI Gene]")
        gene_entry["ncbi_gene"] = fetch_ncbi_gene(
            sym, email, cytoband=cytoband_hint
        )

        # UniProt
        print(f"    [UniProt]")
        gene_entry["uniprot"] = fetch_uniprot(sym)

        # ClinGen Dosage
        gene_entry["clingen_dosage"] = clingen_gene_db.get(sym, {})

        # GenCC
        gencc_entries = gencc_db.get(sym, [])
        # 해당 질환과 관련된 항목 우선, 없으면 전체
        disease_lower = actual_syndrome.lower()
        relevant = [e for e in gencc_entries
                    if disease_lower in e["disease_title"].lower()]
        gene_entry["gencc"] = relevant if relevant else gencc_entries[:5]

        # HPO
        gene_entry["hpo_terms"] = hpo_db.get(sym, [])[:20]

        time.sleep(0.5)

    # CoreRegion에 ClinGen region overlap 추가
    for region_entry in features["core_regions"] + features["partial_chromosomes"]:
        overlaps = get_region_overlaps(
            region_entry["chromosome"],
            region_entry["start"],
            region_entry["end"],
            clingen_region_db,
        )
        region_entry["clingen_region"] = overlaps

    # 9. feature별 critical region 논문
    print(f"\n[9/9] Feature별 논문 수집 (critical region / CoreGene)")
    feature_literature = fetch_all_feature_literature(
        features, actual_syndrome, email, max_per_query=10
    )

    # feature 결과에 문헌 병합
    for gene_entry in features["core_genes"]:
        sym = gene_entry["feature_name"]
        if sym in feature_literature:
            gene_entry["feature_literature"] = feature_literature[sym]["literature"]

    for region_entry in features["core_regions"]:
        fname = region_entry["feature_name"]
        if fname in feature_literature:
            region_entry["feature_literature"] = feature_literature[fname]["literature"]

    for region_entry in features["partial_chromosomes"]:
        fname = region_entry["feature_name"]
        if fname in feature_literature:
            region_entry["feature_literature"] = feature_literature[fname]["literature"]

    # 최종 결과
    result = {
        "nipt_id":   actual_nipt_id,
        "syndrome":  actual_syndrome,
        "nipt_group":features["nipt_group"],
        "markerset": features,
        "medgen":    medgen,
        "annotation_sources": {
            "ncbi_gene":       "NCBI Gene ESearch + ESummary (공개)",
            "uniprot":         "UniProt REST API reviewed Swiss-Prot (CC BY 4.0)",
            "clingen_dosage":  f"ClinGen Dosage Sensitivity GRCh38 ({clingen_gene_path or 'N/A'})",
            "gencc":           f"GenCC submissions TSV ({gencc_path or 'N/A'}) (CC0)",
            "hpo":             f"HPO phenotype_to_genes ({hpo_path or 'N/A'}) (CC BY 4.0)",
            "medgen":          "NCBI MedGen HTML scraping (공개)",
            "pubmed":          "NCBI PubMed EFetch (공개)",
        },
        "collected_at": __import__("datetime").datetime.utcnow().isoformat() + "Z",
    }
    return result


# ════════════════════════════════════════════════════════════════
# TSV 출력
# ════════════════════════════════════════════════════════════════

def export_genes_tsv(result: dict, output_path: str):
    """CoreGene 행 → flat TSV"""
    rows = []
    base = {
        "nipt_id":   result["nipt_id"],
        "syndrome":  result["syndrome"],
        "nipt_group":result["nipt_group"],
    }
    med = result.get("medgen", {})
    xrefs = med.get("cross_references", {})

    for gene in result["markerset"]["core_genes"]:
        sym = gene["feature_name"]
        ng  = gene.get("ncbi_gene", {})
        up  = gene.get("uniprot", {})
        cd  = gene.get("clingen_dosage", {})
        gc  = gene.get("gencc", [{}])
        gc0 = gc[0] if gc else {}

        row = {
            **base,
            "chromosome":         gene["chromosome"],
            "feature_name":       sym,
            "feature_type":       "CoreGene",
            "genomic_pos_start":  gene["start"],
            "genomic_pos_end":    gene["end"],
            "size_mb":            gene["size_mb"],
            # MedGen cross-refs
            "mondo_id":           xrefs.get("MONDO", ""),
            "omim_id":            xrefs.get("OMIM", ""),
            "orphanet_id":        xrefs.get("Orphanet", ""),
            "snomed_ct":          xrefs.get("SNOMED-CT", ""),
            # NCBI Gene
            "ncbi_gene_id":       ng.get("ncbi_gene_id", ""),
            "gene_full_name":     ng.get("gene_full_name", ""),
            "hgnc_id":            ng.get("hgnc_id", ""),
            "ensembl_gene_id":    ng.get("ensembl_gene_id", ""),
            "map_location":       ng.get("map_location", ""),
            "gene_summary":       ng.get("gene_summary", "")[:500],
            # UniProt
            "uniprot_accession":  up.get("uniprot_accession", ""),
            "protein_name":       up.get("protein_name", ""),
            "protein_function":   up.get("protein_function", "")[:500],
            "protein_domains":    "; ".join(
                f"{d['name']}({d['start']}-{d['end']})"
                for d in up.get("protein_domains", [])
            ),
            # ClinGen Dosage
            "hi_score":           cd.get("hi_score", ""),
            "hi_score_label":     cd.get("hi_score_label", ""),
            "hi_disease":         cd.get("hi_disease_name", ""),
            "ts_score":           cd.get("ts_score", ""),
            "ts_score_label":     cd.get("ts_score_label", ""),
            # GenCC
            "gencc_classification":gc0.get("classification", ""),
            "gencc_moi":           gc0.get("moi", ""),
            "gencc_submitters":    "; ".join(
                set(e.get("submitter", "") for e in gc if e.get("submitter"))
            ),
            "gencc_pmids":         gc0.get("pmids", ""),
            # HPO
            "hpo_terms":           "; ".join(
                f"{h['hpo_id']}:{h['hpo_name']}"
                for h in gene.get("hpo_terms", [])[:10]
            ),
        }
        rows.append(row)

    if not rows:
        return
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()),
                                delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"  → {output_path} ({len(rows)}행)")


def export_regions_tsv(result: dict, output_path: str):
    """CoreRegion/PartialChromosome/TargetChromosome 행 → TSV"""
    rows = []
    base = {
        "nipt_id":    result["nipt_id"],
        "syndrome":   result["syndrome"],
        "nipt_group": result["nipt_group"],
    }
    xrefs = result.get("medgen", {}).get("cross_references", {})

    all_regions = (
        result["markerset"]["target_chromosomes"] +
        result["markerset"]["partial_chromosomes"] +
        result["markerset"]["core_regions"]
    )

    for reg in all_regions:
        clingen_hits = reg.get("clingen_region", [])
        ch = clingen_hits[0] if clingen_hits else {}
        row = {
            **base,
            "chromosome":        reg["chromosome"],
            "feature_name":      reg["feature_name"],
            "feature_type":      "Region",
            "genomic_pos_start": reg["start"],
            "genomic_pos_end":   reg["end"],
            "size_mb":           reg["size_mb"],
            "mondo_id":          xrefs.get("MONDO", ""),
            "omim_id":           xrefs.get("OMIM", ""),
            "orphanet_id":       xrefs.get("Orphanet", ""),
            # ClinGen region overlap
            "clingen_region_name": ch.get("region_name", ""),
            "region_hi_score":     ch.get("hi_score", ""),
            "region_hi_label":     ch.get("hi_score_label", ""),
            "region_hi_disease":   ch.get("hi_disease_id", ""),
            "region_ts_score":     ch.get("ts_score", ""),
            "region_ts_label":     ch.get("ts_score_label", ""),
            "clingen_overlaps":    len(clingen_hits),
        }
        rows.append(row)

    if not rows:
        return
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()),
                                delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"  → {output_path} ({len(rows)}행)")


# ════════════════════════════════════════════════════════════════
# CLI
# ════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Syndrome → multi-DB annotation (JSON + TSV)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
사전 다운로드 (선택, 없으면 해당 소스 skip):
  # ClinGen Dosage Sensitivity
  wget https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv
  wget https://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh38.tsv

  # GenCC (CC0, OMIM 제외)
  wget "https://search.thegencc.org/download/action/submissions-export-tsv" -O gencc_submissions.tsv

  # HPO
  wget https://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt

예시:
  python syndrome_annotator.py \\
      --syndrome "DiGeorge syndrome" \\
      --markerset markerset.tsv \\
      --email you@example.com \\
      --clingen-gene ClinGen_gene_curation_list_GRCh38.tsv \\
      --clingen-region ClinGen_region_curation_list_GRCh38.tsv \\
      --gencc gencc_submissions.tsv \\
      --hpo phenotype_to_genes.txt \\
      --output-dir ./output

  python syndrome_annotator.py --nipt-id NIPT_PWS --markerset markerset.tsv --email you@example.com
        """
    )
    parser.add_argument("--syndrome",        type=str, help="Syndrome 이름 (부분 매칭)")
    parser.add_argument("--nipt-id",         type=str, help="NIPT_ID 직접 지정 (예: NIPT_DGS)")
    parser.add_argument("--markerset",       required=True, help="마커셋 TSV 경로")
    parser.add_argument("--email",           default="your@email.com")
    parser.add_argument("--clingen-gene",    default="", help="ClinGen gene dosage TSV")
    parser.add_argument("--clingen-region",  default="", help="ClinGen region dosage TSV")
    parser.add_argument("--gencc",           default="", help="GenCC submissions TSV")
    parser.add_argument("--hpo",             default="", help="HPO phenotype_to_genes.txt")
    parser.add_argument("--output-dir",      default=".", help="출력 디렉토리")
    args = parser.parse_args()

    if not args.syndrome and not args.nipt_id:
        parser.error("--syndrome 또는 --nipt-id 중 하나는 필수")

    global SESSION
    SESSION = _session(args.email)

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    result = annotate_syndrome(
        syndrome_name     = args.syndrome,
        nipt_id           = args.nipt_id,
        markerset_path    = args.markerset,
        email             = args.email,
        clingen_gene_path = args.clingen_gene,
        clingen_region_path = args.clingen_region,
        gencc_path        = args.gencc,
        hpo_path          = args.hpo,
    )

    sid = result["nipt_id"].lower()
    out = Path(args.output_dir)

    # JSON
    json_path = out / f"{sid}_annotated.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)

    # TSV
    genes_tsv   = out / f"{sid}_genes.tsv"
    regions_tsv = out / f"{sid}_regions.tsv"
    export_genes_tsv(result, str(genes_tsv))
    export_regions_tsv(result, str(regions_tsv))

    print(f"\n▶ 완료")
    print(f"  JSON   : {json_path}")
    print(f"  genes  : {genes_tsv}")
    print(f"  regions: {regions_tsv}")

    # 요약
    mg  = result.get("medgen", {})
    lit = mg.get("literature", {})
    n_papers = sum(
        len(v.get("papers", [])) for v in lit.values()
        if isinstance(v, dict)
    )
    n_abs = sum(
        1 for v in lit.values() if isinstance(v, dict)
        for p in v.get("papers", []) if p.get("abstract")
    )
    print(f"\n  질환     : {result['syndrome']} ({result['nipt_id']})")
    print(f"  CoreGene : {len(result['markerset']['core_genes'])}개")
    print(f"  논문     : {n_papers}편 (abstract {n_abs}편)")
    print(f"  lit 카테고리: {list(lit.keys())}")


if __name__ == "__main__":
    main()


# ════════════════════════════════════════════════════════════════
# 9. Critical Region / CoreGene 논문 수집 (PubMed ESearch 기반)
# ════════════════════════════════════════════════════════════════

def _build_feature_queries(
    syndrome: str,
    feature_name: str,
    feature_type: str,
    chrom: str,
) -> dict[str, str]:
    """
    feature_type별 PubMed ESearch 쿼리 생성.

    CoreGene       → 유전자-질환 연관 + 유전자 기능 쿼리
    CoreRegion /
    PrimaryTargetRegion → critical region + genotype-phenotype 쿼리
    PartialChromosome   → partial trisomy/deletion + 표현형 쿼리
    TargetChromosome    → 전체 염색체 review 쿼리
    """
    queries: dict[str, str] = {}
    chrom_num = chrom.replace("chr", "").upper()

    if feature_type == "CoreGene":
        queries["gene_disease"] = (
            f'("{feature_name}"[Gene Name] OR "{feature_name}"[tiab]) '
            f'AND ("{syndrome}"[MeSH Terms] OR "{syndrome}"[tiab]) '
            f'AND ("haploinsufficiency"[tiab] OR "deletion"[tiab] OR '
            f'"critical"[tiab] OR "dosage sensitivity"[tiab])'
        )
        queries["gene_function"] = (
            f'"{feature_name}"[Gene Name] '
            f'AND ("function"[tiab] OR "mechanism"[tiab] OR '
            f'"pathway"[tiab] OR "expression"[tiab]) '
            f'AND "Homo sapiens"[Organism]'
        )

    elif feature_type in ("CoreRegion", "PrimaryTargetRegion"):
        queries["critical_region"] = (
            f'("{feature_name}"[tiab] OR "chromosome {chrom_num}"[tiab]) '
            f'AND ("{syndrome}"[MeSH Terms] OR "{syndrome}"[tiab]) '
            f'AND ("critical region"[tiab] OR "minimal region"[tiab] OR '
            f'"smallest region"[tiab] OR "genotype"[tiab] OR "deletion"[tiab])'
        )
        queries["genotype_phenotype"] = (
            f'("{feature_name}"[tiab]) '
            f'AND ("phenotype"[tiab] OR "genotype-phenotype"[tiab] OR '
            f'"clinical features"[tiab] OR "clinical manifestation"[tiab]) '
            f'AND ("{syndrome}"[MeSH Terms] OR "{syndrome}"[tiab])'
        )

    elif feature_type == "PartialChromosome":
        queries["partial_region"] = (
            f'("{feature_name}"[tiab] OR '
            f'"partial trisomy {chrom_num}"[tiab] OR '
            f'"partial monosomy {chrom_num}"[tiab]) '
            f'AND ("{syndrome}"[MeSH Terms] OR "{syndrome}"[tiab]) '
            f'AND ("phenotype"[tiab] OR "critical region"[tiab])'
        )

    elif feature_type == "TargetChromosome":
        queries["chromosome_review"] = (
            f'"trisomy {chrom_num}"[tiab] '
            f'AND ("{syndrome}"[MeSH Terms] OR "{syndrome}"[tiab]) '
            f'AND "review"[Publication Type]'
        )

    return queries


def fetch_feature_literature(
    syndrome: str,
    feature_name: str,
    feature_type: str,
    chrom: str,
    email: str,
    max_per_query: int = 10,
) -> dict[str, list[dict]]:
    """
    feature 하나에 대한 PubMed 논문 수집.

    반환:
      { query_type: [{ pmid, title, authors, journal, year, abstract,
                       abstract_sections, mesh_terms, pdf_links,
                       full_text_available, doi, pmc_id }] }
    """
    queries = _build_feature_queries(syndrome, feature_name, feature_type, chrom)
    results: dict[str, list[dict]] = {}

    for qtype, query in queries.items():
        try:
            # ESearch
            r = _get(
                f"{EUTILS}/esearch.fcgi",
                params={
                    "db": "pubmed", "term": query,
                    "retmax": max_per_query,
                    "sort": "relevance",
                    "retmode": "json", "email": email,
                },
            )
            pmids = r.json().get("esearchresult", {}).get("idlist", [])
            if not pmids:
                results[qtype] = []
                continue

            # EFetch → abstract
            abstracts = fetch_pubmed_abstracts(pmids, email)
            results[qtype] = [
                abstracts[pmid] for pmid in pmids if pmid in abstracts
            ]
            print(f"    [{qtype}] {len(results[qtype])}편")

        except Exception as e:
            print(f"    [{qtype}] 오류: {e}")
            results[qtype] = []

        time.sleep(0.35)

    return results


def fetch_all_feature_literature(
    features: dict,
    syndrome: str,
    email: str,
    max_per_query: int = 10,
) -> dict[str, dict]:
    """
    마커셋 전체 feature에 대해 critical region 논문 수집.

    반환:
      {
        feature_name: {
          feature_type: str,
          literature: { query_type: [articles] }
        }
      }
    """
    all_lit: dict[str, dict] = {}

    # CoreGene
    for gene in features.get("core_genes", []):
        sym = gene["feature_name"]
        print(f"  [PubMed] CoreGene '{sym}' 논문 수집")
        lit = fetch_feature_literature(
            syndrome, sym, "CoreGene",
            gene["chromosome"], email, max_per_query,
        )
        all_lit[sym] = {"feature_type": "CoreGene", "literature": lit}

    # CoreRegion + PrimaryTargetRegion
    for region in features.get("core_regions", []):
        fname = region["feature_name"]
        ftype = region.get("feature_type", "CoreRegion")
        print(f"  [PubMed] {ftype} '{fname}' 논문 수집")
        lit = fetch_feature_literature(
            syndrome, fname, ftype,
            region["chromosome"], email, max_per_query,
        )
        all_lit[fname] = {"feature_type": ftype, "literature": lit}

    # PartialChromosome
    for region in features.get("partial_chromosomes", []):
        fname = region["feature_name"]
        print(f"  [PubMed] PartialChromosome '{fname}' 논문 수집")
        lit = fetch_feature_literature(
            syndrome, fname, "PartialChromosome",
            region["chromosome"], email, max_per_query,
        )
        all_lit[fname] = {"feature_type": "PartialChromosome", "literature": lit}

    return all_lit