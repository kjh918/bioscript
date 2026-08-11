"""
syndrome_annotator.py
──────────────────────
질환명 기반 syndrome annotation orchestrator.

필수 입력: 질환명(--syndrome)
선택 입력: 기존 마커셋 TSV(--markerset)

마커셋에 질환이 없거나 --markerset을 생략해도 MedGen/PubMed annotation은 수행한다.
마커가 있는 경우에만 chromosome/region/gene annotation을 추가 수행한다.

출력 ID: cbNIPT_<disease_name>
Feature ID: chromosome → region → gene 순서로 001, 002, 003 ... 부여

Usage:
  python -m resource_parse_pipeline.syndrome_annotator \
      --syndrome "DiGeorge syndrome" \
      --email you@example.com \
      --output-dir ./output

  python -m resource_parse_pipeline.syndrome_annotator \
      --syndrome "DiGeorge syndrome" \
      --markerset markerset.tsv \
      --email you@example.com \
      --output-dir ./output
"""
import os
import sys 
import csv
import json
import re
import time
import argparse
from pathlib import Path
from typing import Any

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from resource_parse_pipeline import utils
from resource_parse_pipeline.parsers.medgen_parser   import fetch_medgen
from resource_parse_pipeline.parsers.pubmed_parser   import (fetch_pubmed_abstracts,
                                                    collect_pmids_from_literature,
                                                    discover_and_fetch_disease_literature)
from resource_parse_pipeline.parsers.clingen_parser import (
    parse_clingen_gene,
    parse_clingen_region,
    get_region_overlaps,
    find_clingen_genes,
    find_clingen_regions,
)
from resource_parse_pipeline.annotators.gene_annotator        import annotate_core_genes
from resource_parse_pipeline.annotators.region_annotator      import annotate_regions
from resource_parse_pipeline.annotators.chromosome_annotator  import annotate_all_chromosomal


# ════════════════════════════════════════════════════════════════
# 질환/Feature ID
# ════════════════════════════════════════════════════════════════

def make_cbnipt_disease_name(syndrome: str) -> str:
    """질환명을 파일/ID에서 안전하게 사용할 수 있는 cbNIPT 명으로 변환."""
    cleaned = re.sub(r"[^0-9A-Za-z가-힣]+", "_", syndrome.strip())
    cleaned = re.sub(r"_+", "_", cleaned).strip("_")
    if not cleaned:
        raise ValueError("유효한 질환명이 필요합니다.")
    return f"cbNIPT_{cleaned}"


def empty_features(syndrome: str, cbnipt_disease: str, nipt_group: str = "") -> dict[str, Any]:
    """마커셋이 없는 신규 질환용 빈 feature 구조."""
    return {
        "cbnipt_disease":      cbnipt_disease,
        "syndrome":            syndrome,
        "nipt_group":          nipt_group,
        "target_chromosomes":  [],
        "partial_chromosomes": [],
        "core_regions":        [],
        "core_genes":          [],
    }


def assign_feature_ids(features: dict[str, Any]) -> None:
    """
    chromosome → region → gene 순서로 feature_id를 001부터 부여한다.

    chromosome: TargetChromosome, PartialChromosome
    region: CoreRegion, PrimaryTargetRegion
    gene: CoreGene
    """
    ordered = (
        features.get("target_chromosomes", [])
        + features.get("partial_chromosomes", [])
        + features.get("core_regions", [])
        + features.get("core_genes", [])
    )
    for idx, feature in enumerate(ordered, start=1):
        feature["feature_id"] = f"{idx:03d}"


# ════════════════════════════════════════════════════════════════
# 마커셋 파서
# ════════════════════════════════════════════════════════════════

def load_markerset(tsv_path: str) -> list[dict]:
    with open(tsv_path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return [{k.strip(): v.strip() for k, v in row.items()}
                for row in reader]


def find_syndrome_rows(rows: list[dict],
                        syndrome: str | None,
                        nipt_id: str | None) -> list[dict]:
    out = []
    for row in rows:
        if nipt_id and row.get("NIPT_ID", "").upper() == nipt_id.upper():
            out.append(row)
        elif syndrome and syndrome.lower() in row.get("SYNDROME", "").lower():
            out.append(row)
    return out


def parse_markerset_features(rows: list[dict],
                             fallback_syndrome: str = "") -> dict[str, Any]:
    """마커셋 행을 feature_type별로 분류하고 cbNIPT 질환명을 생성."""
    if not rows:
        syndrome = fallback_syndrome.strip()
        return empty_features(syndrome, make_cbnipt_disease_name(syndrome))

    first = rows[0]
    syndrome = first.get("SYNDROME", "").strip() or fallback_syndrome.strip()
    cbnipt_disease = make_cbnipt_disease_name(syndrome)
    result: dict[str, Any] = {
        "cbnipt_disease":      cbnipt_disease,
        "syndrome":            syndrome,
        "nipt_group":          first.get("NIPT_GROUP", ""),
        "target_chromosomes":  [],
        "partial_chromosomes": [],
        "core_regions":        [],
        "core_genes":          [],
    }

    for row in rows:
        ftype = row.get("FEATURE_TYPE", "")
        entry = {
            "feature_name": row.get("FEATURE_NAME", ""),
            "chromosome":   row.get("CHROMOSOME", ""),
            "start":        int(row.get("GENOMIC_POS_START", 0) or 0),
            "end":          int(row.get("GENOMIC_POS_END", 0) or 0),
            "size_mb":      float(row.get("SIZE_Mb", 0) or 0),
            "feature_type": ftype,
        }
        if ftype == "TargetChromosome":
            result["target_chromosomes"].append(entry)
        elif ftype == "PartialChromosome":
            result["partial_chromosomes"].append(entry)
        elif ftype in ("CoreRegion", "PrimaryTargetRegion"):
            result["core_regions"].append(entry)
        elif ftype == "CoreGene":
            result["core_genes"].append(entry)

    assign_feature_ids(result)
    return result


def _extract_allowed_chromosomes(syndrome: str) -> set[str]:
    """질환명에 직접 명시된 chromosome 집합을 반환한다."""
    allowed: set[str] = set()
    for match in re.finditer(
        r"(?<![A-Za-z0-9])(?:chr)?(\d{1,2}|X|Y)"
        r"[pq]\d+(?:\.\d+)?(?:[-–][pq]?\d+(?:\.\d+)?)?",
        str(syndrome or ""),
        re.I,
    ):
        allowed.add(f"chr{match.group(1).upper()}")

    if not allowed:
        match = re.search(
            r"\b(?:trisomy|monosomy|chromosome)\s+(\d{1,2}|X|Y)\b",
            str(syndrome or ""),
            re.I,
        )
        if match:
            allowed.add(f"chr{match.group(1).upper()}")
    return allowed


def filter_discovered_features(
    discovered: dict[str, list[dict]],
    syndrome: str,
) -> dict[str, list[dict]]:
    """
    parser 결과에 대한 최종 방어 필터.

    TargetChromosome/CoreRegion은 질환명에 명시된 chromosome만 허용한다.
    primary chromosome이 명확할 때 CoreGene도 같은 chromosome만 유지한다.
    """
    allowed = _extract_allowed_chromosomes(syndrome)
    result = {
        "target_chromosomes": list(discovered.get("target_chromosomes", [])),
        "core_regions": list(discovered.get("core_regions", [])),
        "core_genes": list(discovered.get("core_genes", [])),
    }
    if not allowed:
        return result

    result["target_chromosomes"] = [
        item for item in result["target_chromosomes"]
        if item.get("chromosome") in allowed
    ]
    result["core_regions"] = [
        item for item in result["core_regions"]
        if item.get("chromosome") in allowed
    ]
    result["core_genes"] = [
        item for item in result["core_genes"]
        if not item.get("chromosome") or item.get("chromosome") in allowed
    ]
    return result


def merge_discovered_features(
    features: dict[str, Any],
    discovered: dict[str, list[dict]],
) -> None:
    """마커셋 feature와 MedGen feature를 중복 없이 병합한다."""
    for bucket in ("target_chromosomes", "core_regions", "core_genes"):
        existing = features.setdefault(bucket, [])
        seen = {
            (
                str(item.get("feature_name", "")).upper(),
                str(item.get("chromosome", "")).upper(),
            )
            for item in existing
        }
        for item in discovered.get(bucket, []):
            key = (
                str(item.get("feature_name", "")).upper(),
                str(item.get("chromosome", "")).upper(),
            )
            if key not in seen:
                existing.append(item)
                seen.add(key)



def _infer_cnv_event(*texts: str) -> str:
    """질환명/동의어/정의에서 loss/gain 방향을 추론한다."""
    text = " ".join(str(x or "") for x in texts).lower()
    if re.search(r"\b(deletion|deleted|monosomy|loss|minus)\b", text):
        return "loss"
    if re.search(r"\b(duplication|duplicated|trisomy|tetrasomy|gain)\b", text):
        return "gain"
    return "unknown"


def _collect_target_cytobands(features: dict[str, Any], medgen: dict[str, Any]) -> list[str]:
    """기존 feature와 MedGen gene location에서 ClinGen 검색용 cytoband를 수집한다."""
    bands: list[str] = []

    for region in features.get("core_regions", []):
        for key in ("feature_name", "cytoband"):
            value = str(region.get(key, "") or "").strip()
            if value:
                bands.append(value)

    for gene in features.get("core_genes", []):
        value = str(gene.get("cytoband", "") or "").strip()
        if value:
            bands.append(value)

    for gene in medgen.get("gene_locations", []) or []:
        value = str(gene.get("cytoband", "") or "").strip()
        if value:
            bands.append(value)

    for target in medgen.get("genomic_targets", []) or []:
        value = str(target.get("cytoband", "") or "").strip()
        if value:
            bands.append(value)

    seen: set[str] = set()
    output: list[str] = []
    for band in bands:
        key = band.upper().replace(" ", "")
        if key and key not in seen:
            seen.add(key)
            output.append(band)
    return output


def _chrom_from_record(record: dict[str, Any]) -> str:
    chrom = str(record.get("chrom", "") or record.get("chromosome", "") or "").strip()
    if chrom:
        return chrom if chrom.startswith("chr") else f"chr{chrom}"
    band = str(record.get("cytoband", "") or "").strip()
    match = re.match(r"^(?:chr)?(\d{1,2}|X|Y)", band, re.I)
    return f"chr{match.group(1).upper()}" if match else ""


def merge_clingen_discovery(
    features: dict[str, Any],
    clingen_genes: list[dict[str, Any]],
    clingen_regions: list[dict[str, Any]],
) -> None:
    """ClinGen 질환/cytoband 검색 결과를 신규 CoreGene/CoreRegion feature로 반영한다."""
    discovered = {
        "target_chromosomes": [],
        "core_regions": [],
        "core_genes": [],
    }

    for region in clingen_regions:
        chrom = _chrom_from_record(region)
        band = str(region.get("cytoband", "") or "").strip()
        name = band or str(region.get("region_name", "") or "").strip()
        if not chrom or not name:
            continue
        discovered["core_regions"].append({
            "feature_name": name,
            "chromosome": chrom,
            "start": int(region.get("start", 0) or 0),
            "end": int(region.get("end", 0) or 0),
            "size_mb": round((int(region.get("end", 0) or 0) - int(region.get("start", 0) or 0) + 1) / 1_000_000, 3)
                       if int(region.get("start", 0) or 0) > 0 and int(region.get("end", 0) or 0) > 0 else 0.0,
            "feature_type": "CoreRegion",
            "source": "ClinGen",
            "clingen_region": [dict(region)],
        })

    for gene in clingen_genes:
        symbol = str(gene.get("gene_symbol", "") or "").strip()
        if not symbol:
            continue
        chrom = _chrom_from_record(gene)
        discovered["core_genes"].append({
            "feature_name": symbol,
            "chromosome": chrom,
            "cytoband": str(gene.get("cytoband", "") or ""),
            "start": int(gene.get("start", 0) or 0),
            "end": int(gene.get("end", 0) or 0),
            "size_mb": 0.0,
            "feature_type": "CoreGene",
            "source": "ClinGen",
            "clingen_dosage": dict(gene),
        })

    chromosomes = {
        _chrom_from_record(item)
        for item in [*clingen_regions, *clingen_genes]
        if _chrom_from_record(item)
    }
    for chrom in sorted(chromosomes):
        discovered["target_chromosomes"].append({
            "feature_name": chrom,
            "chromosome": chrom,
            "start": 0,
            "end": 0,
            "size_mb": 0.0,
            "feature_type": "TargetChromosome",
            "source": "ClinGen disease/cytoband match",
        })

    merge_discovered_features(features, discovered)

# ════════════════════════════════════════════════════════════════
# 메인 annotation
# ════════════════════════════════════════════════════════════════

def annotate_syndrome(
    syndrome_name:      str | None,
    email:              str,
    nipt_id:            str | None = None,
    markerset_path:     str = "",
    clingen_gene_path:  str = "",
    clingen_region_path:str = "",
    gencc_path:         str = "",
    hpo_path:           str = "",
    uniprot_dat:        str = "",
    max_lit:            int = 10,
    verify_pmc:         bool = True,
) -> dict[str, Any]:

    # ── [1] 질환 및 선택적 마커셋 로드 ───────────────────────
    if not syndrome_name:
        raise ValueError("신규 질환 annotation에는 syndrome_name이 필요합니다.")

    actual_syndrome = syndrome_name.strip()
    rows: list[dict] = []
    if markerset_path:
        print(f"\n[1] 마커셋 로드(선택): {markerset_path}")
        all_rows = load_markerset(markerset_path)
        rows = find_syndrome_rows(all_rows, actual_syndrome, nipt_id)
        if not rows:
            print(f"  → 마커셋에 '{actual_syndrome}' 없음: 신규 질환 모드로 계속")
    else:
        print("\n[1] 마커셋 없음: 신규 질환 모드")

    features = parse_markerset_features(rows, fallback_syndrome=actual_syndrome)
    actual_syndrome = features["syndrome"]
    cbnipt_disease  = features["cbnipt_disease"]
    nipt_group      = features["nipt_group"]
    print(f"  → {actual_syndrome} ({cbnipt_disease}) [{nipt_group or '미지정'}]")

    # ── [2] MedGen + 신규 질환 feature discovery ─────────────
    print(f"\n[2] MedGen")
    medgen = fetch_medgen(actual_syndrome, email)

    discovered = medgen.get("discovered_features", {})
    discovered = filter_discovered_features(discovered, actual_syndrome)

    print(
        "  [MedGen feature] target="
        f"{[x.get('chromosome') for x in discovered.get('target_chromosomes', [])]} "
        "region="
        f"{[x.get('feature_name') for x in discovered.get('core_regions', [])]} "
        "gene="
        f"{[x.get('feature_name') for x in discovered.get('core_genes', [])]}"
    )

    merge_discovered_features(features, discovered)
    assign_feature_ids(features)

    print(f"     CoreGene={len(features['core_genes'])} "
          f"CoreRegion={len(features['core_regions'])} "
          f"TargetChromosome={len(features['target_chromosomes'])}")

    # ── [3] PubMed disease discovery + abstracts ─────────────
    print(f"\n[3] PubMed disease discovery + abstracts")

    # MedGen 화면에 노출된 PMID만 사용하지 않고, MedGen에서 얻은
    # canonical disease/synonym을 seed로 PubMed를 독립 검색한다.
    # discovery 단계에서는 PMC 역검증을 끄고 API 호출량을 줄인다.
    pubmed_discovery = discover_and_fetch_disease_literature(
        disease_name=actual_syndrome,
        synonyms=medgen.get("synonyms", []) or [],
        email=email,
        medgen_literature=medgen.get("literature", {}),
        bucket_limits={
            "review": min(max_lit, 5),
            "genotype_phenotype": max_lit,
            "genetics": max_lit,
            "cohort": min(max_lit, 5),
            "general": max_lit,
        },
        verify_pmc=False,
    )

    all_pmids = pubmed_discovery.get("pmids", [])
    pubmed_lookup = pubmed_discovery.get("articles", {})

    print(
        f"  [PubMed] merged PMID={len(all_pmids)} "
        f"articles={len(pubmed_lookup)}"
    )

    # MedGen literature에 포함된 paper에는 기존처럼 abstract를 병합한다.
    # 독립 discovery paper 전체는 result["pubmed_discovery"]에 별도 보존한다.
    for cat_data in medgen.get("literature", {}).values():
        if isinstance(cat_data, dict):
            for paper in cat_data.get("papers", []):
                pmid = paper.get("pmid", "")
                if pmid and pmid in pubmed_lookup:
                    ab = pubmed_lookup[pmid]
                    paper.update({
                        "abstract":           ab.get("abstract", ""),
                        "abstract_sections":  ab.get("abstract_sections", {}),
                        "mesh_terms":         ab.get("mesh_terms", []),
                        "pdf_links":          ab.get("pdf_links", []),
                        "full_text_available":ab.get("full_text_available", False),
                        "pmc_verify":         ab.get("pmc_verify"),
                        "doi":                ab.get("doi", ""),
                        "pmc_id":             ab.get("pmc_id", ""),
                    })

    # ── 외부 DB 로드 (선택) ───────────────────────────────────
    clingen_gene_db   = parse_clingen_gene(clingen_gene_path) if clingen_gene_path else {}
    clingen_region_db = parse_clingen_region(clingen_region_path) if clingen_region_path else []

    # ClinGen 질환명/xref/cytoband 검색 결과를 실제 feature에 반영
    clingen_matched_genes: list[dict[str, Any]] = []
    clingen_matched_regions: list[dict[str, Any]] = []
    if clingen_gene_db or clingen_region_db:
        xrefs = medgen.get("cross_references", {}) or {}
        synonyms = medgen.get("synonyms", []) or []
        definition = medgen.get("definition", "")
        event = _infer_cnv_event(actual_syndrome, " ".join(synonyms), definition)
        target_bands = _collect_target_cytobands(features, medgen)

        clingen_matched_genes = find_clingen_genes(
            genes=clingen_gene_db,
            syndrome_name=actual_syndrome,
            mondo_id=xrefs.get("MONDO", ""),
            omim_id=xrefs.get("OMIM", ""),
            orphanet_id=xrefs.get("Orphanet", ""),
            cytobands=target_bands,
            event=event,
        ) if clingen_gene_db else []

        clingen_matched_regions = find_clingen_regions(
            regions=clingen_region_db,
            syndrome_name=actual_syndrome,
            mondo_id=xrefs.get("MONDO", ""),
            omim_id=xrefs.get("OMIM", ""),
            orphanet_id=xrefs.get("Orphanet", ""),
            cytobands=target_bands,
            event=event,
        ) if clingen_region_db else []

        merge_clingen_discovery(
            features,
            clingen_matched_genes,
            clingen_matched_regions,
        )
        assign_feature_ids(features)

        print(
            f"  [ClinGen discovery] event={event} bands={target_bands} "
            f"genes={len(clingen_matched_genes)} regions={len(clingen_matched_regions)}"
        )

    gencc_db: dict = {}
    if gencc_path:
        from resource_parse_pipeline.parsers.gencc_parser import load_gencc
        gencc_db = load_gencc(gencc_path)

    hpo_db: dict = {}
    if hpo_path:
        from resource_parse_pipeline.parsers.hpo_parser import load_hpo
        hpo_db = load_hpo(hpo_path)

    # UniProt local DB
    uniprot_local_db = None
    if uniprot_dat:
        from resource_parse_pipeline.parsers.uniprot_parser import UniProtLocalDB
        uniprot_local_db = UniProtLocalDB(uniprot_dat)

    # ── [4] CoreGene annotation ──────────────────────────────
    if features["core_genes"]:
        print(f"\n[4] CoreGene annotation ({len(features['core_genes'])}개)")
        # ClinGen / GenCC / HPO 외부 DB에서 미리 hint 추출
        for gene in features["core_genes"]:
            sym = gene["feature_name"]
            cd  = clingen_gene_db.get(sym, {})
            if cd.get("cytoband") and not gene.get("cytoband"):
                gene["cytoband"] = cd["cytoband"]
            gene["clingen_dosage"] = cd
            gene["gencc"]          = gencc_db.get(sym, [])
            gene["hpo_terms"]      = hpo_db.get(sym, [])[:20]

        annotate_core_genes(features["core_genes"], email, uniprot_local_db)

        # feature_literature
        from resource_parse_pipeline.annotators.gene_annotator import fetch_ncbi_gene
        for gene in features["core_genes"]:
            sym = gene["feature_name"]
            print(f"  [Gene lit] '{sym}'")
            gene["feature_literature"] = _fetch_gene_literature(
                sym, actual_syndrome, email, max_lit)

    # ── [5] CoreRegion annotation ────────────────────────────
    if features["core_regions"]:
        print(f"\n[5] CoreRegion annotation ({len(features['core_regions'])}개)")
        annotate_regions(
            features["core_regions"],
            actual_syndrome, email,
            clingen_region_db=clingen_region_db,
            max_lit=max_lit,
        )

    # ── [6] TargetChromosome annotation ─────────────────────
    chromosomal_annotation = []
    if features["target_chromosomes"] or features["partial_chromosomes"]:
        print(f"\n[6] TargetChromosome annotation")
        chromosomal_annotation = annotate_all_chromosomal(
            features, cbnipt_disease, clingen_region_db)

    return {
        "cbnipt_disease": cbnipt_disease,
        "syndrome":   actual_syndrome,
        "nipt_group": nipt_group,
        "markerset":  features,
        "chromosomal_annotation": chromosomal_annotation,
        "medgen":     medgen,
        "pubmed_discovery": pubmed_discovery,
        "clingen_discovery": {
            "matched_genes": clingen_matched_genes,
            "matched_regions": clingen_matched_regions,
        },
        "annotation_sources": {
            "medgen":          "NCBI MedGen HTML (공개)",
            "pubmed":          "NCBI PubMed ESearch+EFetch disease discovery (공개)",
            "ncbi_gene":       "NCBI Gene ESearch+ESummary (공개)",
            "uniprot":         f"UniProt {'로컬파일' if uniprot_dat else 'ExPASy/REST'}",
            "clingen_dosage":  clingen_gene_path or "N/A",
            "gencc":           gencc_path or "N/A",
            "hpo":             hpo_path or "N/A",
        },
        "collected_at": __import__("datetime").datetime.utcnow().isoformat() + "Z",
    }


def _fetch_gene_literature(gene_symbol: str, syndrome: str,
                            email: str, max_results: int) -> dict:
    """CoreGene 관련 PubMed 논문 수집."""
    import xml.etree.ElementTree as ET
    from resource_parse_pipeline.parsers.pubmed_parser import _parse_single_article

    queries = {
        "gene_disease": (
            f'("{gene_symbol}"[Gene Name] OR "{gene_symbol}"[tiab]) '
            f'AND ("{syndrome}"[MeSH Terms] OR "{syndrome}"[tiab]) '
            f'AND ("haploinsufficiency"[tiab] OR "deletion"[tiab] OR '
            f'"pathogenic"[tiab] OR "critical"[tiab])'
        ),
        "gene_function": (
            f'"{gene_symbol}"[Gene Name] '
            f'AND ("function"[tiab] OR "mechanism"[tiab] OR "pathway"[tiab]) '
            f'AND "Homo sapiens"[Organism]'
        ),
    }
    result: dict = {}
    for qtype, query in queries.items():
        try:
            r = utils.http_get(f"{utils.EUTILS}/esearch.fcgi",
                               params={"db": "pubmed", "term": query,
                                       "retmax": max_results,
                                       "sort": "relevance",
                                       "retmode": "json", "email": email})
            pmids = r.json().get("esearchresult", {}).get("idlist", [])
            if not pmids:
                result[qtype] = []
                continue
            time.sleep(0.35)
            r2 = utils.http_get(f"{utils.EUTILS}/efetch.fcgi",
                                params={"db": "pubmed", "id": ",".join(pmids),
                                        "rettype": "abstract", "retmode": "xml",
                                        "email": email})
            root = ET.fromstring(r2.text.encode("utf-8"))
            articles = [_parse_single_article(pa)
                        for pa in root.findall(".//PubmedArticle")]
            result[qtype] = [a for a in articles if a]
            print(f"    [{qtype}] {len(result[qtype])}편")
        except Exception as e:
            print(f"    [{qtype}] 오류: {e}")
            result[qtype] = []
        time.sleep(0.35)
    return result


# ════════════════════════════════════════════════════════════════
# TSV 출력
# ════════════════════════════════════════════════════════════════

def export_genes_tsv(result: dict, output_path: str):
    rows = []
    syndrome = result["syndrome"]
    xrefs    = result.get("medgen", {}).get("cross_references", {})

    for gene in result["markerset"]["core_genes"]:
        sym = gene["feature_name"]
        ng  = gene.get("ncbi_gene", {})
        up  = gene.get("uniprot", {})
        cd  = gene.get("clingen_dosage", {})
        gc_list = gene.get("gencc", [])
        gc0 = gc_list[0] if gc_list else {}

        rows.append({
            "cbnipt_disease":      result["cbnipt_disease"],
            "syndrome":            syndrome,
            "nipt_group":          result["nipt_group"],
            "chromosome":          gene.get("chromosome", ""),
            "feature_id":          gene.get("feature_id", ""),
            "feature_name":        sym,
            "feature_type":        "CoreGene",
            "genomic_pos_start":   ng.get("chr_start", ""),
            "genomic_pos_end":     ng.get("chr_stop", ""),
            "cytoband":            gene.get("cytoband", ""),
            "mondo_id":            xrefs.get("MONDO", ""),
            "omim_id":             xrefs.get("OMIM", ""),
            "orphanet_id":         xrefs.get("Orphanet", ""),
            "ncbi_gene_id":        ng.get("ncbi_gene_id", ""),
            "gene_full_name":      ng.get("gene_full_name", ""),
            "hgnc_id":             ng.get("hgnc_id", ""),
            "gene_summary":        ng.get("gene_summary", "")[:400],
            "uniprot_accession":   up.get("uniprot_accession", ""),
            "protein_name":        up.get("protein_name", ""),
            "protein_function":    up.get("protein_function", "")[:400],
            "protein_domains":     "; ".join(
                f"{d['name']}({d.get('start','')}-{d.get('end','')})"
                for d in up.get("protein_domains", [])
            ),
            "hi_score":            cd.get("hi_score", ""),
            "hi_score_label":      cd.get("hi_score_label", ""),
            "hi_disease":          cd.get("hi_disease_name", ""),
            "ts_score":            cd.get("ts_score", ""),
            "gencc_classification":gc0.get("classification", ""),
            "gencc_moi":           gc0.get("moi", ""),
            "hpo_terms":           "; ".join(
                f"{h['hpo_id']}:{h['hpo_name']}"
                for h in gene.get("hpo_terms", [])[:10]
            ),
        })

    if not rows:
        return
    import csv as _csv
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = _csv.DictWriter(f, fieldnames=list(rows[0].keys()),
                                  delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"  → {output_path} ({len(rows)}행)")


def export_regions_tsv(result: dict, output_path: str):
    rows = []
    xrefs = result.get("medgen", {}).get("cross_references", {})
    all_regions = sorted(
        result["markerset"]["target_chromosomes"]
        + result["markerset"]["partial_chromosomes"]
        + result["markerset"]["core_regions"],
        key=lambda x: x.get("feature_id", "999"),
    )
    for reg in all_regions:
        ch = (reg.get("clingen_region") or [{}])[0]
        coords = reg.get("grch38_coords", {})  # TargetChromosome
        rows.append({
            "cbnipt_disease":      result["cbnipt_disease"],
            "syndrome":            result["syndrome"],
            "nipt_group":          result["nipt_group"],
            "chromosome":          reg.get("chromosome", ""),
            "feature_id":          reg.get("feature_id", ""),
            "feature_name":        reg.get("feature_name", ""),
            "feature_type":        reg.get("feature_type", ""),
            "genomic_pos_start":   reg.get("start") or coords.get("start", ""),
            "genomic_pos_end":     reg.get("end") or coords.get("end", ""),
            "size_mb":             reg.get("size_mb") or coords.get("length_mb", ""),
            "mechanism":           reg.get("mechanism", ""),
            "karyotype":           reg.get("karyotype", ""),
            "mondo_id":            xrefs.get("MONDO", ""),
            "omim_id":             xrefs.get("OMIM", ""),
            "clingen_region_name": ch.get("region_name", ""),
            "region_hi_score":     ch.get("hi_score", ""),
            "region_hi_label":     ch.get("hi_score_label", ""),
            "region_hi_disease":   ch.get("hi_disease_id", ""),
            "region_ts_score":     ch.get("ts_score", ""),
            "assembly":            coords.get("assembly", ""),
            "clinical_notes":      reg.get("clinical_notes", "")[:300],
        })

    if not rows:
        return
    import csv as _csv
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = _csv.DictWriter(f, fieldnames=list(rows[0].keys()),
                                  delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"  → {output_path} ({len(rows)}행)")


# ════════════════════════════════════════════════════════════════
# CLI
# ════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Syndrome annotation (질환명 기반, 마커셋 선택)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
예시:
  # 신규 질환: 마커셋 없이 실행
  python -m resource_parse_pipeline.syndrome_annotator \
      --syndrome "DiGeorge syndrome" \
      --email you@example.com \
      --output-dir ./output

  # 기존 마커셋을 함께 사용하는 경우
  python -m resource_parse_pipeline.syndrome_annotator \
      --syndrome "DiGeorge syndrome" \
      --markerset markerset.tsv \
      --email you@example.com
        """
    )
    parser.add_argument("--syndrome",        type=str, required=True)
    parser.add_argument("--nipt-id",         type=str, default=None,
                        help="선택 사항: 기존 마커셋 검색용 ID")
    parser.add_argument("--markerset",       type=str, default="",
                        help="선택 사항: 기존 마커셋 TSV")
    parser.add_argument("--email",           default="your@email.com")
    parser.add_argument("--clingen-gene",    default="")
    parser.add_argument("--clingen-region",  default="")
    parser.add_argument("--gencc",           default="")
    parser.add_argument("--hpo",             default="")
    parser.add_argument("--uniprot-dat",     default="",
                        help="UniProt Swiss-Prot human .dat[.gz] 파일")
    parser.add_argument("--max-lit",         type=int, default=10)
    parser.add_argument("--no-verify-pmc",   action="store_true")
    parser.add_argument("--output-dir",      default=".")
    args = parser.parse_args()

    utils.init_session(args.email)
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    result = annotate_syndrome(
        syndrome_name=args.syndrome,
        nipt_id=args.nipt_id,
        markerset_path=args.markerset,
        email=args.email,
        clingen_gene_path=args.clingen_gene,
        clingen_region_path=args.clingen_region,
        gencc_path=args.gencc,
        hpo_path=args.hpo,
        uniprot_dat=args.uniprot_dat,
        max_lit=args.max_lit,
        verify_pmc=not args.no_verify_pmc,
    )

    sid = result["cbnipt_disease"]
    out = Path(args.output_dir)

    json_path   = out / f"{sid}_annotated.json"
    genes_path  = out / f"{sid}_genes.tsv"
    regions_path= out / f"{sid}_regions.tsv"

    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)

    export_genes_tsv(result, str(genes_path))
    export_regions_tsv(result, str(regions_path))

    print(f"\n▶ 완료")
    print(f"  JSON   : {json_path}")
    print(f"  genes  : {genes_path}")
    print(f"  regions: {regions_path}")

    stats = {
        "core_genes":    len(result["markerset"]["core_genes"]),
        "core_regions":  len(result["markerset"]["core_regions"]),
        "target_chroms": len(result["chromosomal_annotation"]),
    }
    print(f"  stats: {stats}")


if __name__ == "__main__":
    main()