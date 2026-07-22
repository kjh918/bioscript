"""
syndrome_discovery.py
─────────────────────
Syndrome 이름만으로 연관 유전자/영역을 발견하고
각 유전자의 임상적 근거를 multi-DB 교차 검증하여
마커셋 형태의 JSON + TSV로 출력합니다.

수집 소스 (우선순위 순):
  1. Orphanet en_product6.xml  → 유전자 + 역할 타입 + 관계 상태
  2. GenCC TSV                 → 유전자-질환 분류 강도 (Definitive/Strong/...)
  3. ClinGen Dosage FTP        → HI score / TS score
  4. HPO phenotype_to_genes    → 표현형 연관 유전자
  5. NCBI Gene ESummary        → 유전자 기능 요약, map_location
  6. UniProt REST              → 단백질 기능, 도메인
  7. PubMed ESearch            → 유전자별 근거 논문

출력:
  {syndrome_id}_discovered.json   ← 전체 구조
  {syndrome_id}_discovered_genes.tsv  ← flat TSV (마커셋과 동일 구조 + evidence 컬럼)

Usage:
  python syndrome_discovery.py \\
      --syndrome "DiGeorge syndrome" \\
      --email you@example.com \\
      --orphanet en_product6.xml \\
      --gencc gencc_submissions.tsv \\
      --clingen-gene ClinGen_gene_curation_list_GRCh38.tsv \\
      --hpo phenotype_to_genes.txt \\
      --output-dir ./output
"""

import re
import os
import csv
import sys
import json
import time
import argparse
import xml.etree.ElementTree as ET
from io import StringIO
from pathlib import Path
from typing import Any
from collections import defaultdict

import requests
from requests.adapters import HTTPAdapter
from clingen_parser import parse_clingen_gene, parse_clingen_region, get_region_overlaps
from urllib3.util.retry import Retry

# ── 상수 ──────────────────────────────────────────────────────
EUTILS  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UNIPROT = "https://rest.uniprot.org/uniprotkb"
PUBMED  = "https://pubmed.ncbi.nlm.nih.gov"

# GenCC 분류 강도 → 숫자 점수
GENCC_SCORE = {
    "Definitive":           4,
    "Strong":               3,
    "Moderate":             2,
    "Limited":              1,
    "Animal Model Only":    1,
    "Disputed Evidence":    0,
    "Refuted Evidence":    -1,
    "No Known Disease Relationship": 0,
}

# ClinGen HI score → 숫자 점수
HI_SCORE_MAP = {
    "3": 3,   # Sufficient Evidence
    "2": 2,   # Emerging Evidence
    "1": 1,   # Little Evidence
    "0": 0,   # No Evidence
    "30": 1,  # AR gene
    "40": 0,  # Unlikely
}
HI_LABEL = {
    "3": "Sufficient Evidence",
    "2": "Emerging Evidence",
    "1": "Little Evidence",
    "0": "No Evidence",
    "30": "Associated with AR Phenotype",
    "40": "Dosage Sensitivity Unlikely",
}

# Orphanet 역할 타입 → 점수
ORPHANET_ROLE_SCORE = {
    "Disease-causing germline mutation(s) in": 4,
    "Disease-causing germline mutation(s) (loss of function) in": 4,
    "Disease-causing germline mutation(s) (gain of function) in": 4,
    "Disease-causing somatic mutation(s) in": 3,
    "Modifying germline mutation in": 2,
    "Major susceptibility factor in": 2,
    "Role in the phenotype of": 1,
    "Candidate gene tested in": 1,
    "Part of a fusion gene in": 1,
}

SESSION: requests.Session | None = None


def _make_session(email: str) -> requests.Session:
    s = requests.Session()
    retry = Retry(total=5, backoff_factor=1.5,
                  status_forcelist=[429, 500, 502, 503, 504],
                  allowed_methods=["GET"])
    s.mount("https://", HTTPAdapter(max_retries=retry))
    s.mount("http://",  HTTPAdapter(max_retries=retry))
    s.headers.update({"User-Agent": f"syndrome-discovery/1.0 (mailto:{email})"})
    return s


def _get(url: str, params: dict | None = None, delay: float = 0.35,
         **kw) -> requests.Response:
    time.sleep(delay)
    r = SESSION.get(url, params=params, timeout=30, **kw)
    r.raise_for_status()
    return r


# ════════════════════════════════════════════════════════════════
# Step 1. MedGen → 질환 ID 획득
# ════════════════════════════════════════════════════════════════

def resolve_disease_ids(syndrome_name: str, email: str) -> dict:
    """
    syndrome 이름 → MedGen UID, MONDO, OMIM, Orphanet ID,
                    gene_locations (Gene(location) 행에서 직접 추출)
    """
    print(f"  [MedGen] '{syndrome_name}' ID 조회")
    try:
        r = _get(f"{EUTILS}/esearch.fcgi",
                 params={"db": "medgen", "term": f'"{syndrome_name}"',
                         "retmax": 1, "sort": "relevance",
                         "retmode": "json", "email": email})
        ids = r.json().get("esearchresult", {}).get("idlist", [])
        if not ids:
            return {}

        uid = ids[0]
        time.sleep(1.0)
        resp = SESSION.get(f"https://www.ncbi.nlm.nih.gov/medgen/{uid}",
                           timeout=20, headers={
                               "User-Agent": "Mozilla/5.0 Chrome/124"})
        resp.raise_for_status()

        from bs4 import BeautifulSoup
        soup = BeautifulSoup(resp.text, "lxml")

        xrefs: dict[str, str] = {}
        gene_locations: list[dict] = []

        table = soup.find("table", class_="medgenTable")
        if table:
            for tr in table.find_all("tr"):
                tds = tr.find_all("td")
                if len(tds) < 2:
                    continue
                label = tds[0].get_text().strip().rstrip(":").lower()
                a_tags = tds[1].find_all("a")
                first_a = a_tags[0] if a_tags else None
                val = (first_a.get_text() if first_a
                       else tds[1].get_text()).strip()

                # ── Gene (location) 행 ──────────────────────
                if "gene" in label and "location" in label:
                    for a in a_tags:
                        sym = a.get_text(strip=True)
                        if not sym:
                            continue
                        href = a.get("href", "")
                        gid_m = re.search(r"/gene/(\d+)", href)
                        gene_id = gid_m.group(1) if gid_m else ""

                        cytoband = ""
                        for sibling in a.next_siblings:
                            cb_m = re.search(
                                r"\(([0-9XY]+[pq][\d.]+(?:[- ][pq]?[\d.]+)?)\)",
                                str(sibling)
                            )
                            if cb_m:
                                cytoband = cb_m.group(1).strip()
                                break

                        gene_locations.append({
                            "gene_symbol":  sym,
                            "ncbi_gene_id": gene_id,
                            "cytoband":     cytoband,
                            "source":       "MedGen",
                        })

                # ── Cross-references ─────────────────────────
                elif "monarch" in label or "mondo" in label:
                    xrefs["MONDO"] = val
                elif "omim" in label:
                    xrefs["OMIM"] = val
                elif "orphanet" in label:
                    xrefs["Orphanet"] = re.sub(r"^ORPHA", "", val)
                elif "snomed" in label:
                    m = re.search(r"\((\d+)\)", val)
                    xrefs["SNOMED-CT"] = m.group(1) if m else val

        print(f"  [MedGen] UID={uid} xrefs={xrefs} "
              f"gene_loc={[g['gene_symbol'] for g in gene_locations]}")
        return {"medgen_uid": uid, "gene_locations": gene_locations, **xrefs}

    except Exception as e:
        print(f"  [MedGen] 오류: {e}")
        return {}


# ════════════════════════════════════════════════════════════════
# Step 2. 소스별 유전자 수집
# ════════════════════════════════════════════════════════════════

def load_orphanet_genes(xml_path: str,
                        orpha_code: str) -> list[dict]:
    """
    Orphanet en_product6.xml → 특정 질환의 연관 유전자 목록.
    OrphaCode로 질환을 찾고 연관 유전자 + 역할 타입 반환.

    반환:
      [{ gene_symbol, hgnc_id, role_type, status }]
    """
    if not xml_path or not Path(xml_path).exists():
        print(f"  [Orphanet] 파일 없음: {xml_path}")
        return []

    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
    except Exception as e:
        print(f"  [Orphanet] XML 파싱 오류: {e}")
        return []

    genes = []
    # OrphaCode 정규화 (앞의 "ORPHA:" 제거)
    target_code = str(orpha_code).replace("ORPHA:", "").strip()

    for disorder in root.iter("Disorder"):
        code_el = disorder.find("OrphaCode")
        if code_el is None or code_el.text != target_code:
            continue

        for assoc in disorder.iter("DisorderGeneAssociation"):
            gene_el    = assoc.find("Gene")
            role_el    = assoc.find("DisorderGeneAssociationType/Name")
            status_el  = assoc.find("DisorderGeneAssociationStatus/Name")

            if gene_el is None:
                continue

            sym_el   = gene_el.find("Symbol")
            # ET의 XPath는 속성 조건에 @ 하나만 사용
            hgnc_el = None
            for xref in gene_el.findall(".//ExternalReference"):
                if xref.findtext("Source") == "HGNC":
                    hgnc_el = xref.find("Reference")
                    break
            if hgnc_el is None:
                hgnc_el = gene_el.find(".//ExternalReference/Reference")

            genes.append({
                "gene_symbol": sym_el.text.strip() if sym_el is not None else "",
                "hgnc_id":     (hgnc_el.text.strip()
                                if hgnc_el is not None else ""),
                "role_type":   (role_el.text.strip()
                                if role_el is not None else ""),
                "status":      (status_el.text.strip()
                                if status_el is not None else ""),
            })

    print(f"  [Orphanet] OrphaCode={target_code} → {len(genes)}개 유전자")
    return genes


def load_gencc_by_disease(tsv_path: str,
                           mondo_id: str = "",
                           syndrome_name: str = "") -> list[dict]:
    """
    GenCC TSV → 특정 질환(MONDO ID 또는 질환명)의 유전자-분류 목록.

    반환:
      [{ gene_symbol, hgnc_id, classification, moi, submitter, pmids }]
    """
    if not tsv_path or not Path(tsv_path).exists():
        print(f"  [GenCC] 파일 없음: {tsv_path}")
        return []

    results = []
    syndrome_lower = syndrome_name.lower()
    # MONDO ID 정규화: "MONDO:0008608" → "MONDO:0008608"
    target_mondo = mondo_id.upper().replace(" ", "")

    with open(tsv_path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            disease_curie = row.get("disease_curie", "").strip().upper()
            disease_title = row.get("disease_title", "").strip().lower()
            orig_curie    = row.get("disease_original_curie", "").strip().upper()

            # MONDO ID 매칭 또는 질환명 부분 매칭
            matched = (
                (target_mondo and (disease_curie == target_mondo or
                                   orig_curie == target_mondo)) or
                (syndrome_lower and syndrome_lower in disease_title)
            )
            if not matched:
                continue

            sym = row.get("gene_symbol", "").strip()
            if not sym:
                continue

            classification = row.get("classification_title", "").strip()
            if GENCC_SCORE.get(classification, 0) < 0:
                # Refuted는 제외
                continue

            results.append({
                "gene_symbol":    sym,
                "hgnc_id":        row.get("gene_curie", "").strip(),
                "disease_curie":  row.get("disease_curie", "").strip(),
                "disease_title":  row.get("disease_title", "").strip(),
                "classification": classification,
                "moi":            row.get("moi_title", "").strip(),
                "submitter":      row.get("submitter_title", "").strip(),
                "pmids":          row.get("submitted_as_pmids", "").strip(),
                "report_url":     row.get("submitted_as_public_report_url", "").strip(),
            })

    print(f"  [GenCC] {len(results)}개 유전자-분류 (MONDO={target_mondo})")
    return results


def load_clingen_dosage_by_disease(
    gene_tsv: str,
    region_tsv: str,
    syndrome_name: str = "",
    mondo_id: str = "",
    omim_id: str = "",
) -> tuple[list[dict], list[dict]]:
    """
    ClinGen Dosage FTP → 질환과 연관된 유전자 + 영역 목록.
    clingen_parser.parse_clingen_gene/region 기반으로 인덱스 파싱.
    """
    def _match(text: str) -> bool:
        t = (text or "").lower()
        return (
            (syndrome_name and syndrome_name.lower() in t) or
            (omim_id and omim_id in t) or
            (mondo_id and mondo_id.lower() in t)
        )

    # ── 유전자 ───────────────────────────────────────────────
    all_genes = parse_clingen_gene(gene_tsv)
    matched_genes = [
        g for g in all_genes.values()
        if _match(g["hi_disease_name"]) or
           _match(g["ts_disease_name"]) or
           _match(g["hi_disease_id"])
    ]

    # ── 영역 ────────────────────────────────────────────────
    all_regions = parse_clingen_region(region_tsv)
    matched_regions = [
        r for r in all_regions
        if _match(r["hi_disease_name"]) or
           _match(r["ts_disease_name"]) or
           _match(r["hi_disease_id"])
    ]

    print(f"  [ClinGen Dosage] 유전자={len(matched_genes)}, 영역={len(matched_regions)}")
    return matched_genes, matched_regions


def load_hpo_by_disease(hpo_path: str,
                         syndrome_name: str = "",
                         omim_id: str = "") -> list[dict]:
    """
    HPO phenotype_to_genes.txt → 질환 연관 유전자.
    disease_id(OMIM) 또는 질환명으로 필터.

    반환: [{ gene_symbol, hpo_id, hpo_name, disease_id }]
    """
    if not hpo_path or not Path(hpo_path).exists():
        print(f"  [HPO] 파일 없음: {hpo_path}")
        return []

    # phenotype_to_genes.txt 포맷:
    # HPO-id  HPO label  entrez-gene-id  entrez-gene-symbol  ...  disease-ID
    # 또는 genes_to_disease.txt 포맷:
    # entrez-gene-id  entrez-gene-symbol  Association  disease-ID  disease-name
    results = []
    omim_target = f"OMIM:{omim_id}" if omim_id and not omim_id.startswith("OMIM:") else omim_id
    syn_lower = syndrome_name.lower()

    with open(hpo_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue

            # 포맷 자동 감지
            if parts[0].startswith("HP:"):
                # phenotype_to_genes.txt
                hpo_id, hpo_name = parts[0], parts[1]
                gene_sym = parts[3] if len(parts) > 3 else ""
                disease_id = parts[5] if len(parts) > 5 else ""
            else:
                # genes_to_disease.txt (gene_id, symbol, assoc, disease_id, disease_name)
                gene_sym   = parts[1] if len(parts) > 1 else ""
                disease_id = parts[3] if len(parts) > 3 else ""
                dis_name   = parts[4] if len(parts) > 4 else ""
                hpo_id, hpo_name = "", dis_name

            matched = (
                (omim_target and omim_target in disease_id) or
                (syn_lower and syn_lower in hpo_name.lower())
            )
            if matched and gene_sym:
                results.append({
                    "gene_symbol": gene_sym,
                    "hpo_id":      hpo_id,
                    "hpo_name":    hpo_name,
                    "disease_id":  disease_id,
                })

    # 유전자 기준 집계
    gene_hpo: dict[str, list] = defaultdict(list)
    for r in results:
        gene_hpo[r["gene_symbol"]].append({
            "hpo_id":   r["hpo_id"],
            "hpo_name": r["hpo_name"],
        })

    aggregated = [
        {"gene_symbol": sym, "hpo_terms": terms}
        for sym, terms in gene_hpo.items()
    ]
    print(f"  [HPO] {len(aggregated)}개 유전자")
    return aggregated


# ════════════════════════════════════════════════════════════════
# Step 3. 교차 검증 + Evidence Tier 산정
# ════════════════════════════════════════════════════════════════

def cross_validate_genes(
    orphanet_genes:   list[dict],
    gencc_genes:      list[dict],
    clingen_genes:    list[dict],
    hpo_genes:        list[dict],
    medgen_gene_locs: list[dict] | None = None,
) -> list[dict]:
    """
    5개 소스 병합 + evidence_tier 산정.
    MedGen Gene(location) 직접 명시 → 가장 높은 신뢰도.
    """
    orpha_map: dict[str, dict] = {
        g["gene_symbol"]: g for g in orphanet_genes if g.get("gene_symbol")
    }
    gencc_map: dict[str, list[dict]] = defaultdict(list)
    for g in gencc_genes:
        gencc_map[g["gene_symbol"]].append(g)
    clingen_map: dict[str, dict] = {
        g["gene_symbol"]: g for g in clingen_genes if g.get("gene_symbol")
    }
    hpo_map: dict[str, dict] = {
        g["gene_symbol"]: g for g in hpo_genes if g.get("gene_symbol")
    }
    medgen_map: dict[str, dict] = {
        g["gene_symbol"]: g
        for g in (medgen_gene_locs or [])
        if g.get("gene_symbol")
    }

    all_symbols = (set(orpha_map) | set(gencc_map) |
                   set(clingen_map) | set(hpo_map) | set(medgen_map))

    merged: list[dict] = []
    for sym in all_symbols:
        orpha   = orpha_map.get(sym, {})
        gencc_l = gencc_map.get(sym, [])
        clingen = clingen_map.get(sym, {})
        hpo     = hpo_map.get(sym, {})
        medgen  = medgen_map.get(sym, {})

        orpha_score  = ORPHANET_ROLE_SCORE.get(orpha.get("role_type", ""), 0)
        gencc_score  = max(
            (GENCC_SCORE.get(g["classification"], 0) for g in gencc_l),
            default=0
        )
        hi_score_num = HI_SCORE_MAP.get(clingen.get("hi_score", ""), 0)
        hpo_score    = 1 if hpo else 0
        medgen_score = 3 if medgen else 0   # MedGen 직접 명시 가중치

        source_count = sum([
            bool(orpha), bool(gencc_l), bool(clingen),
            bool(hpo),   bool(medgen)
        ])

        evidence_score = (orpha_score * 2 + gencc_score * 3 +
                          hi_score_num * 2 + hpo_score +
                          medgen_score * 3 + source_count * 2)

        # MedGen 직접 명시 → 무조건 Tier1
        if medgen or gencc_score >= 3 or source_count >= 3:
            tier = "Tier1_High"
        elif source_count == 2 or gencc_score == 2:
            tier = "Tier2_Medium"
        else:
            tier = "Tier3_Low"

        # cytoband / ncbi_gene_id 우선순위
        cytoband        = (clingen.get("cytoband") or medgen.get("cytoband") or "")
        ncbi_gene_hint  = (clingen.get("ncbi_gene_id") or medgen.get("ncbi_gene_id") or "")

        best_gencc = max(gencc_l, key=lambda g: GENCC_SCORE.get(
            g["classification"], 0), default={})

        merged.append({
            "gene_symbol":         sym,
            "evidence_tier":       tier,
            "evidence_score":      evidence_score,
            "source_count":        source_count,
            "sources_confirmed":   sorted([
                s for s, m in [
                    ("MedGen",   medgen),
                    ("Orphanet", orpha),
                    ("GenCC",    gencc_l),
                    ("ClinGen",  clingen),
                    ("HPO",      hpo),
                ] if m
            ]),
            # MedGen
            "medgen_gene_loc":     medgen.get("cytoband", ""),
            "medgen_ncbi_gene_id": medgen.get("ncbi_gene_id", ""),
            # Orphanet
            "orphanet_role":        orpha.get("role_type", ""),
            "orphanet_status":      orpha.get("status", ""),
            # GenCC
            "gencc_classification": best_gencc.get("classification", ""),
            "gencc_moi":            best_gencc.get("moi", ""),
            "gencc_submitters":     "; ".join(
                {g["submitter"] for g in gencc_l if g.get("submitter")}
            ),
            "gencc_pmids":          best_gencc.get("pmids", ""),
            # ClinGen Dosage
            "hi_score":             clingen.get("hi_score", ""),
            "hi_score_label":       clingen.get("hi_score_label", ""),
            "hi_disease":           clingen.get("hi_disease_name", ""),
            "ts_score":             clingen.get("ts_score", ""),
            "ts_score_label":       clingen.get("ts_score_label", ""),
            "cytoband":             cytoband,
            "hgnc_id":              "",
            # HPO
            "hpo_terms":            hpo.get("hpo_terms", []),
        })

    merged.sort(key=lambda x: -x["evidence_score"])
    return merged


# ════════════════════════════════════════════════════════════════
# Step 4. NCBI Gene + UniProt annotation
# ════════════════════════════════════════════════════════════════

def fetch_ncbi_gene(gene_symbol: str, email: str,
                    hgnc_id: str = "", cytoband: str = "") -> dict:
    """
    NCBI Gene에서 Human(Homo sapiens) 유전자만 정확히 조회.

    [방어 전략 3단계]
    1순위: HGNC ID가 있으면 Gene Accession으로 직접 검색 (가장 정확)
    2순위: Gene Name + Homo sapiens[Organism] + 염색체 번호(cytoband 파싱)
    3순위: Gene Name + Homo sapiens[Organism] → ESummary에서 taxid=9606 검증

    마우스/기타 생물 ortholog가 잘못 잡히는 문제를 방지합니다.
    """
    def _esummary(gene_id: str) -> dict:
        time.sleep(0.35)
        r = _get(f"{EUTILS}/esummary.fcgi",
                 params={"db": "gene", "id": gene_id,
                         "retmode": "json", "email": email})
        return r.json().get("result", {}).get(gene_id, {})

    def _is_human(data: dict) -> bool:
        """ESummary 결과가 Human(taxid=9606)인지 검증"""
        org = data.get("organism", {})
        taxid = org.get("taxid", 0)
        sci_name = org.get("scientificname", "").lower()
        # taxid 또는 scientific name 중 하나로 확인
        return (taxid == 9606 or "homo sapiens" in sci_name)

    def _build_result(gene_id: str, data: dict) -> dict:
        ginfo = (data.get("genomicinfo", [{}])[0]
                 if data.get("genomicinfo") else {})
        # chraccver: GRCh38 = NC_000022.11 형태, 확인용
        accver = ginfo.get("chraccver", "")
        return {
            "ncbi_gene_id":   gene_id,
            "gene_full_name": data.get("description", ""),
            "gene_summary":   data.get("summary", ""),
            "hgnc_id":        data.get("nomenclatureauthorityid", ""),
            "map_location":   data.get("maplocation", ""),
            "chromosome":     data.get("chromosome", ""),
            "chr_start":      ginfo.get("chrstart", ""),
            "chr_stop":       ginfo.get("chrstop", ""),
            "chr_accver":     accver,   # GRCh38 확인용
        }

    try:
        # ── 1순위: HGNC ID 직접 검색 ──────────────────────────
        # GenCC gene_curie: "HGNC:11532" → "11532" 추출
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
                    print(f"    [NCBI Gene] {gene_symbol} → "
                          f"HGNC:{hgnc_num} → gene_id={ids[0]} "
                          f"chr={data.get('chromosome','')} "
                          f"loc={data.get('maplocation','')}")
                    return _build_result(ids[0], data)

        # ── 2순위: Gene Name + Organism + 염색체 번호 ──────────
        # cytoband "22q11.21" → chrom_num "22"
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
            ids = r.json().get("esearchresult", {}).get("idlist", [])
            for gene_id in ids:
                data = _esummary(gene_id)
                if _is_human(data):
                    print(f"    [NCBI Gene] {gene_symbol} → "
                          f"chr{chrom_num} filter → gene_id={gene_id} "
                          f"chr={data.get('chromosome','')} "
                          f"loc={data.get('maplocation','')}")
                    return _build_result(gene_id, data)

        # ── 3순위: Gene Name + Organism → taxid 검증 ───────────
        query = (f"{gene_symbol}[Gene Name] AND Homo sapiens[Organism]")
        r = _get(f"{EUTILS}/esearch.fcgi",
                 params={"db": "gene", "term": query,
                         "retmax": 5, "sort": "relevance",
                         "retmode": "json", "email": email})
        ids = r.json().get("esearchresult", {}).get("idlist", [])
        for gene_id in ids:
            data = _esummary(gene_id)
            if _is_human(data):
                print(f"    [NCBI Gene] {gene_symbol} → "
                      f"organism filter → gene_id={gene_id} "
                      f"chr={data.get('chromosome','')} "
                      f"loc={data.get('maplocation','')}")
                return _build_result(gene_id, data)

        print(f"    [NCBI Gene] {gene_symbol} → Human 결과 없음")
        return {}

    except Exception as e:
        print(f"    [NCBI Gene] {gene_symbol} 오류: {e}")
        return {}


def fetch_uniprot(gene_symbol: str) -> dict:
    try:
        r = _get(UNIPROT, params={
            "query": (f"gene_exact:{gene_symbol} AND "
                      f"organism_id:9606 AND reviewed:true"),
            "format": "json", "size": 1,
        })
        results = r.json().get("results", [])
        if not results:
            return {}
        entry = results[0]
        acc   = entry.get("primaryAccession", "")
        func  = next((
            c["texts"][0]["value"]
            for c in entry.get("comments", [])
            if c.get("commentType") == "FUNCTION" and c.get("texts")
        ), "")
        domains = [
            {"name": f.get("description", ""),
             "start": f.get("location", {}).get("start", {}).get("value", ""),
             "end":   f.get("location", {}).get("end",   {}).get("value", "")}
            for f in entry.get("features", [])
            if f.get("type") == "Domain"
        ]
        disease_assoc = [
            {"name": c.get("disease", {}).get("diseaseId", ""),
             "desc": (c.get("texts") or [{}])[0].get("value", "")}
            for c in entry.get("comments", [])
            if c.get("commentType") == "DISEASE"
        ]
        return {
            "uniprot_accession":    acc,
            "protein_name":         (entry.get("proteinDescription", {})
                                     .get("recommendedName", {})
                                     .get("fullName", {}).get("value", "")),
            "protein_function":     func,
            "protein_domains":      domains,
            "disease_associations": disease_assoc,
        }
    except Exception as e:
        print(f"    [UniProt] {gene_symbol} 오류: {e}")
        return {}


# ════════════════════════════════════════════════════════════════
# Step 5. 유전자별 PubMed 논문 수집
# ════════════════════════════════════════════════════════════════

def fetch_gene_literature(gene_symbol: str, syndrome_name: str,
                           email: str, max_results: int = 10) -> dict:
    """
    유전자 심볼 + syndrome → 근거 논문 수집.
    gene_disease / gene_function 두 쿼리 사용.
    """
    queries = {
        "gene_disease": (
            f'("{gene_symbol}"[Gene Name] OR "{gene_symbol}"[tiab]) '
            f'AND ("{syndrome_name}"[MeSH Terms] OR "{syndrome_name}"[tiab]) '
            f'AND ("haploinsufficiency"[tiab] OR "deletion"[tiab] OR '
            f'"pathogenic"[tiab] OR "critical"[tiab])'
        ),
        "gene_function": (
            f'"{gene_symbol}"[Gene Name] '
            f'AND ("function"[tiab] OR "mechanism"[tiab] OR "pathway"[tiab]) '
            f'AND "Homo sapiens"[Organism]'
        ),
    }

    import xml.etree.ElementTree as ET

    literature: dict[str, list[dict]] = {}

    for qtype, query in queries.items():
        try:
            r = _get(f"{EUTILS}/esearch.fcgi",
                     params={"db": "pubmed", "term": query,
                             "retmax": max_results,
                             "sort": "relevance",
                             "retmode": "json", "email": email})
            pmids = r.json().get("esearchresult", {}).get("idlist", [])
            if not pmids:
                literature[qtype] = []
                continue

            time.sleep(0.35)
            r2 = _get(f"{EUTILS}/efetch.fcgi",
                      params={"db": "pubmed", "id": ",".join(pmids),
                              "rettype": "abstract", "retmode": "xml",
                              "email": email})
            root = ET.fromstring(r2.text.encode("utf-8"))
            articles = []
            for pa in root.findall(".//PubmedArticle"):
                mc  = pa.find("MedlineCitation")
                art = mc.find("Article") if mc is not None else None
                if mc is None or art is None:
                    continue
                pmid_el = mc.find("PMID")
                pmid = pmid_el.text.strip() if pmid_el is not None else ""
                title_el = art.find("ArticleTitle")
                title = ("".join(title_el.itertext()).strip()
                         if title_el is not None else "")
                # abstract
                abs_el = art.find(".//Abstract")
                abstract, sections = "", {}
                if abs_el:
                    parts = []
                    for at in abs_el.findall("AbstractText"):
                        label = at.get("Label") or at.get("NlmCategory") or ""
                        text  = re.sub(r"\s+", " ",
                                       "".join(at.itertext()).strip())
                        if label:
                            sections[label] = text
                            parts.append(f"{label}: {text}")
                        else:
                            parts.append(text)
                    abstract = " ".join(parts)
                # PMC
                pmc_id = ""
                pd = pa.find("PubmedData")
                if pd is not None:
                    for aid in pd.findall(".//ArticleIdList/ArticleId"):
                        if aid.get("IdType") == "pmc":
                            pmc_id = aid.text.strip()
                pdf_links = []
                if pmc_id:
                    pdf_links.append({"source": "PMC", "type": "pdf",
                                      "url": f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmc_id}/pdf/"})
                pdf_links.append({"source": "PubMed", "type": "abstract",
                                   "url": f"{PUBMED}/{pmid}/"})
                articles.append({
                    "pmid": pmid, "title": title,
                    "abstract": abstract,
                    "abstract_sections": sections,
                    "full_text_available": bool(pmc_id),
                    "pmc_id": pmc_id,
                    "pdf_links": pdf_links,
                })
            literature[qtype] = articles
            print(f"    [{qtype}] {len(articles)}편")
        except Exception as e:
            print(f"    [{qtype}] 오류: {e}")
            literature[qtype] = []
        time.sleep(0.35)

    return literature


# ════════════════════════════════════════════════════════════════
# 메인 파이프라인
# ════════════════════════════════════════════════════════════════

def discover_syndrome_genes(
    syndrome_name:      str,
    email:              str,
    orphanet_xml:       str = "",
    gencc_tsv:          str = "",
    clingen_gene_tsv:   str = "",
    clingen_region_tsv: str = "",
    hpo_path:           str = "",
    tier_filter:        str = "Tier1_High,Tier2_Medium,Tier3_Low",
    max_lit:            int = 10,
) -> dict[str, Any]:

    tiers = [t.strip() for t in tier_filter.split(",")]

    # Step 1: 질환 ID + MedGen Gene(location) 직접 수집
    print(f"\n[1/5] 질환 ID 조회")
    disease_ids = resolve_disease_ids(syndrome_name, email)
    orphanet_id     = disease_ids.get("Orphanet", "")
    mondo_id        = disease_ids.get("MONDO", "")
    omim_id         = disease_ids.get("OMIM", "")
    medgen_gene_locs = disease_ids.get("gene_locations", [])

    if medgen_gene_locs:
        print(f"  [MedGen] Gene(location) 직접 수집: "
              f"{[g['gene_symbol'] for g in medgen_gene_locs]}")

    # Step 2: 소스별 유전자 수집
    print(f"\n[2/5] 소스별 유전자 수집")
    orphanet_genes = load_orphanet_genes(orphanet_xml, orphanet_id)
    gencc_genes    = load_gencc_by_disease(gencc_tsv, mondo_id, syndrome_name)
    clingen_genes, clingen_regions = load_clingen_dosage_by_disease(
        clingen_gene_tsv, clingen_region_tsv,
        syndrome_name, mondo_id, omim_id,
    )
    hpo_genes = load_hpo_by_disease(hpo_path, syndrome_name, omim_id)

    # Step 3: 교차 검증 (MedGen gene_locations 포함)
    print(f"\n[3/5] 교차 검증 + Evidence Tier 산정")
    merged = cross_validate_genes(
        orphanet_genes, gencc_genes, clingen_genes, hpo_genes,
        medgen_gene_locs=medgen_gene_locs,
    )
    filtered = [g for g in merged if g["evidence_tier"] in tiers]
    print(f"  전체 {len(merged)}개 → tier 필터 후 {len(filtered)}개")
    for tier in ["Tier1_High", "Tier2_Medium", "Tier3_Low"]:
        n = sum(1 for g in merged if g["evidence_tier"] == tier)
        print(f"  {tier}: {n}개")

    # Step 4: NCBI Gene + UniProt annotation
    print(f"\n[4/5] 유전자 상세 annotation ({len(filtered)}개)")
    for gene in filtered:
        sym      = gene["gene_symbol"]
        hgnc_id  = gene.get("hgnc_id", "")     # GenCC gene_curie
        cytoband = gene.get("cytoband", "")     # ClinGen Dosage cytoBand
        print(f"  → {sym} [{gene['evidence_tier']}]  hgnc={hgnc_id}  cytoband={cytoband}")
        gene["ncbi_gene"] = fetch_ncbi_gene(sym, email,
                                             hgnc_id=hgnc_id,
                                             cytoband=cytoband)
        gene["uniprot"]   = fetch_uniprot(sym)
        # cytoband 보완 (NCBI Gene map_location으로)
        if not gene.get("cytoband") and gene["ncbi_gene"].get("map_location"):
            gene["cytoband"] = gene["ncbi_gene"]["map_location"]
        time.sleep(0.5)

    # Step 5: 논문 수집
    print(f"\n[5/5] 유전자별 논문 수집")
    for gene in filtered:
        sym = gene["gene_symbol"]
        print(f"  → {sym}")
        gene["literature"] = fetch_gene_literature(
            sym, syndrome_name, email, max_results=max_lit
        )

    return {
        "syndrome":          syndrome_name,
        "disease_ids":       disease_ids,
        "clingen_regions":   clingen_regions,
        "genes":             filtered,
        "genes_all_tiers":   merged,
        "annotation_sources": {
            "orphanet":       f"en_product6.xml ({orphanet_xml or 'N/A'})",
            "gencc":          f"GenCC TSV ({gencc_tsv or 'N/A'}) (CC0)",
            "clingen_dosage": f"ClinGen FTP ({clingen_gene_tsv or 'N/A'})",
            "hpo":            f"HPO phenotype_to_genes ({hpo_path or 'N/A'})",
            "ncbi_gene":      "NCBI Gene ESearch+ESummary (공개)",
            "uniprot":        "UniProt REST reviewed (CC BY 4.0)",
            "pubmed":         "PubMed EFetch (공개)",
        },
        "collected_at": __import__("datetime").datetime.utcnow().isoformat() + "Z",
    }


# ════════════════════════════════════════════════════════════════
# TSV 출력
# ════════════════════════════════════════════════════════════════

def export_tsv(result: dict, output_path: str):
    """discovered genes → flat TSV (마커셋 형식 호환 + evidence 컬럼)"""
    rows = []
    syndrome = result["syndrome"]
    dids     = result["disease_ids"]

    for gene in result["genes"]:
        ng = gene.get("ncbi_gene", {})
        up = gene.get("uniprot", {})
        lit = gene.get("literature", {})
        n_papers = sum(len(v) for v in lit.values())
        n_abs    = sum(
            1 for v in lit.values() for a in v if a.get("abstract")
        )

        rows.append({
            # 마커셋 호환 컬럼
            "NIPT_ID":           f"DISCOVERY_{gene['gene_symbol']}",
            "SYNDROME":          syndrome,
            "NIPT_GROUP":        "Discovered",
            "CHROMOSOME":        f"chr{ng.get('chromosome', '')}",
            "FEATURE_NAME":      gene["gene_symbol"],
            "FEATURE_TYPE":      "CoreGene",
            "GENOMIC_POS_START": ng.get("chr_start", ""),
            "GENOMIC_POS_END":   ng.get("chr_stop", ""),
            "SIZE_Mb":           "",
            # 질환 ID
            "MONDO_ID":          dids.get("MONDO", ""),
            "OMIM_ID":           dids.get("OMIM", ""),
            "ORPHANET_ID":       dids.get("Orphanet", ""),
            # Evidence
            "EVIDENCE_TIER":     gene["evidence_tier"],
            "EVIDENCE_SCORE":    gene["evidence_score"],
            "SOURCE_COUNT":      gene["source_count"],
            "SOURCES_CONFIRMED": "; ".join(gene["sources_confirmed"]),
            # Orphanet
            "ORPHANET_ROLE":     gene["orphanet_role"],
            "ORPHANET_STATUS":   gene["orphanet_status"],
            # GenCC
            "GENCC_CLASSIFICATION": gene["gencc_classification"],
            "GENCC_MOI":         gene["gencc_moi"],
            "GENCC_SUBMITTERS":  gene["gencc_submitters"],
            "GENCC_PMIDS":       gene["gencc_pmids"],
            # ClinGen Dosage
            "HI_SCORE":          gene["hi_score"],
            "HI_SCORE_LABEL":    gene["hi_score_label"],
            "HI_DISEASE":        gene["hi_disease"],
            "TS_SCORE":          gene["ts_score"],
            "TS_SCORE_LABEL":    gene["ts_score_label"],
            "CYTOBAND":          gene.get("cytoband", ""),
            # NCBI Gene
            "NCBI_GENE_ID":      ng.get("ncbi_gene_id", ""),
            "GENE_FULL_NAME":    ng.get("gene_full_name", ""),
            "MAP_LOCATION":      ng.get("map_location", ""),
            "GENE_SUMMARY":      ng.get("gene_summary", "")[:400],
            # UniProt
            "UNIPROT_ACCESSION": up.get("uniprot_accession", ""),
            "PROTEIN_NAME":      up.get("protein_name", ""),
            "PROTEIN_FUNCTION":  up.get("protein_function", "")[:400],
            "PROTEIN_DOMAINS":   "; ".join(
                f"{d['name']}({d['start']}-{d['end']})"
                for d in up.get("protein_domains", [])
            ),
            # HPO
            "HPO_TERMS":         "; ".join(
                f"{h['hpo_id']}:{h['hpo_name']}"
                for h in gene.get("hpo_terms", [])[:10]
            ),
            # 논문
            "LITERATURE_COUNT":  n_papers,
            "ABSTRACT_COUNT":    n_abs,
            "GENE_DISEASE_PMIDS": "; ".join(
                a["pmid"] for a in lit.get("gene_disease", []) if a.get("pmid")
            ),
        })

    if not rows:
        print("  출력할 데이터 없음")
        return

    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()),
                                delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"  → {output_path} ({len(rows)}행, {len(rows[0])}컬럼)")


# ════════════════════════════════════════════════════════════════
# CLI
# ════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Syndrome → 유전자 발견 + 근거 annotation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
사전 다운로드:
  wget https://sciences.orphadata.com/genes/ -O en_product6.xml
  wget "https://search.thegencc.org/download/action/submissions-export-tsv" -O gencc_submissions.tsv
  wget https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv
  wget https://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh38.tsv
  wget https://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt

예시:
  python syndrome_discovery.py \\
      --syndrome "DiGeorge syndrome" \\
      --email you@example.com \\
      --orphanet en_product6.xml \\
      --gencc gencc_submissions.tsv \\
      --clingen-gene ClinGen_gene_curation_list_GRCh38.tsv \\
      --clingen-region ClinGen_region_curation_list_GRCh38.tsv \\
      --hpo phenotype_to_genes.txt \\
      --tier Tier1_High,Tier2_Medium \\
      --output-dir ./output
        """
    )
    parser.add_argument("--syndrome",        required=True)
    parser.add_argument("--email",           default="your@email.com")
    parser.add_argument("--orphanet",        default="", help="en_product6.xml")
    parser.add_argument("--gencc",           default="", help="gencc_submissions.tsv")
    parser.add_argument("--clingen-gene",    default="")
    parser.add_argument("--clingen-region",  default="")
    parser.add_argument("--hpo",             default="")
    parser.add_argument("--tier",            default="Tier1_High,Tier2_Medium",
                        help="포함할 tier (기본: Tier1+Tier2)")
    parser.add_argument("--max-lit",         type=int, default=10,
                        help="유전자별 논문 최대 수 (기본: 10)")
    parser.add_argument("--output-dir",      default=".")
    args = parser.parse_args()

    global SESSION
    SESSION = _make_session(args.email)
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    result = discover_syndrome_genes(
        syndrome_name      = args.syndrome,
        email              = args.email,
        orphanet_xml       = args.orphanet,
        gencc_tsv          = args.gencc,
        clingen_gene_tsv   = args.clingen_gene,
        clingen_region_tsv = args.clingen_region,
        hpo_path           = args.hpo,
        tier_filter        = args.tier,
        max_lit            = args.max_lit,
    )

    sid = re.sub(r"[^a-z0-9]", "_", args.syndrome.lower())
    out = Path(args.output_dir)

    json_path = out / f"{sid}_discovered.json"
    tsv_path  = out / f"{sid}_discovered_genes.tsv"

    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)

    export_tsv(result, str(tsv_path))

    print(f"\n▶ 완료")
    print(f"  JSON : {json_path}")
    print(f"  TSV  : {tsv_path}")
    print(f"\n  syndrome   : {result['syndrome']}")
    print(f"  disease_ids: {result['disease_ids']}")
    print(f"  발견 유전자: {len(result['genes'])}개")
    for tier in ["Tier1_High", "Tier2_Medium", "Tier3_Low"]:
        n = sum(1 for g in result["genes"] if g["evidence_tier"] == tier)
        if n:
            print(f"    {tier}: {n}개")
    n_regions = len(result.get("clingen_regions", []))
    if n_regions:
        print(f"  ClinGen 영역: {n_regions}개")


if __name__ == "__main__":
    main()