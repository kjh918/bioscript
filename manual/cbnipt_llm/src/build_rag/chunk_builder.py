"""
chunk_builder.py
────────────────
medgen_results.json + pubmed_abstracts.json → RAG chunk JSON

chunk 유형 (chunk_type):
  disease_overview        질환 정의 + synonyms + cross-refs
  clinical_guideline      professional guidelines 논문 1편
  etiology                etiology 카테고리 논문 1편
  diagnosis               diagnosis 카테고리 논문 1편
  therapy                 therapy 카테고리 논문 1편
  prognosis               prognosis 카테고리 논문 1편
  clinical_prediction     clinical prediction guides 논문 1편
  systematic_review       recent systematic reviews 논문 1편

각 chunk 공통 필드:
  chunk_id        "{disease_uid}_{type}_{pmid or seq}"
  chunk_type      위 유형 중 하나
  disease_name    질환명
  disease_uid     MedGen UID
  concept_id      UMLS Concept ID
  xrefs           MONDO/OMIM/Orphanet/SNOMED-CT
  text            임베딩할 본문 (핵심)
  metadata        검색/필터용 부가 정보

Usage:
  python chunk_builder.py \\
      --medgen medgen_results.json \\
      --pubmed pubmed_abstracts.json \\
      --output chunks.json
"""

import re
import json
import hashlib
import argparse
from typing import Any

# ── 카테고리 매핑 ──────────────────────────────────────────────
CLINICAL_CAT_MAP = {
    "etiology":                  "etiology",
    "diagnosis":                 "diagnosis",
    "therapy":                   "therapy",
    "prognosis":                 "prognosis",
    "clinical_prediction_guides":"clinical_prediction",
}


# ── PubMed lookup 빌드 ─────────────────────────────────────────

def build_pubmed_lookup(pubmed_records: list[dict]) -> dict[str, dict]:
    """pmid → pubmed record 딕셔너리"""
    return {str(r.get("pmid", "")): r for r in pubmed_records if r.get("pmid")}


# ── PDF 링크 추출 ─────────────────────────────────────────────

def get_best_pdf_links(pubmed_rec: dict | None) -> list[dict]:
    """
    type='pdf' 우선, 없으면 'html', 없으면 전체 반환.
    """
    if not pubmed_rec:
        return []
    links = pubmed_rec.get("pdf_links", [])
    pdf_only = [l for l in links if l.get("type") == "pdf"]
    return pdf_only if pdf_only else links


# ── 공통 메타데이터 빌더 ───────────────────────────────────────

def _base_meta(disease: dict) -> dict:
    return {
        "disease_name":  disease.get("disease_name", ""),
        "disease_uid":   disease.get("medgen_uid", ""),
        "concept_id":    disease.get("concept_id", ""),
        "disease_type":  disease.get("disease_type", ""),
        "xrefs":         disease.get("cross_references", {}),
        "medgen_url":    disease.get("url", ""),
    }


def _paper_meta(paper: dict, pubmed_rec: dict | None, category: str) -> dict:
    """논문 1편 메타데이터"""
    pmid = str(paper.get("pmid", ""))
    meta: dict[str, Any] = {
        "pmid":             pmid,
        "category":         category,
        "title":            paper.get("title", "") if pubmed_rec is None else pubmed_rec.get("title", paper.get("title", "")),
        "journal":          pubmed_rec.get("journal", "") if pubmed_rec else paper.get("journal", ""),
        "year":             pubmed_rec.get("year", "") if pubmed_rec else "",
        "authors":          pubmed_rec.get("authors", []) if pubmed_rec else [],
        "publication_types":pubmed_rec.get("publication_types", []) if pubmed_rec else [],
        "mesh_terms":       pubmed_rec.get("mesh_terms", []) if pubmed_rec else [],
        "keywords":         pubmed_rec.get("keywords", []) if pubmed_rec else [],
        "full_text_available": pubmed_rec.get("full_text_available", False) if pubmed_rec else False,
        "pdf_links":        get_best_pdf_links(pubmed_rec),
        "pubmed_url":       f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else "",
        "doi":              pubmed_rec.get("doi", "") if pubmed_rec else "",
        "pmc_id":           pubmed_rec.get("pmc_id", "") if pubmed_rec else "",
    }
    return meta


def _chunk_id(disease_uid: str, chunk_type: str, suffix: str) -> str:
    raw = f"{disease_uid}__{chunk_type}__{suffix}"
    # 길이 제한 + 특수문자 제거
    safe = re.sub(r"[^a-zA-Z0-9_\-]", "_", raw)[:120]
    return safe


# ── chunk 생성 함수들 ──────────────────────────────────────────

def make_overview_chunk(disease: dict) -> dict:
    """
    질환 개요 chunk.
    definition 전문 + synonyms + cross-refs 포함.
    """
    name       = disease.get("disease_name", "")
    uid        = disease.get("medgen_uid", "")
    concept_id = disease.get("concept_id", "")
    definition = disease.get("definition", "")
    synonyms   = disease.get("synonyms", [])
    xrefs      = disease.get("cross_references", {})

    xref_lines = "\n".join(f"  {k}: {v}" for k, v in xrefs.items())
    syn_str    = ", ".join(synonyms) if synonyms else "N/A"

    text = f"""Disease: {name}
MedGen UID: {uid} | Concept ID: {concept_id}
Synonyms: {syn_str}
Cross-references:
{xref_lines}

Definition:
{definition}"""

    return {
        "chunk_id":   _chunk_id(uid, "disease_overview", uid),
        "chunk_type": "disease_overview",
        "disease_name": name,
        "disease_uid":  uid,
        "concept_id":   concept_id,
        "text":         text.strip(),
        "metadata": {
            **_base_meta(disease),
            "synonyms":  synonyms,
            "has_definition": bool(definition),
        },
    }


def make_paper_chunk(
    disease: dict,
    paper: dict,
    chunk_type: str,
    pubmed_lookup: dict[str, dict],
) -> dict | None:
    """
    논문 1편 → chunk.
    abstract가 있으면 full abstract 포함, 없으면 title만.
    """
    pmid = str(paper.get("pmid", "")).strip()
    if not pmid:
        return None

    pubmed_rec  = pubmed_lookup.get(pmid)
    disease_name = disease.get("disease_name", "")
    uid          = disease.get("medgen_uid", "")

    # 제목: pubmed > medgen paper 순
    title = (
        (pubmed_rec.get("title", "") if pubmed_rec else "")
        or paper.get("title", "")
    )

    # Abstract 처리
    abstract = pubmed_rec.get("abstract", "") if pubmed_rec else ""
    abstract_sections = pubmed_rec.get("abstract_sections", {}) if pubmed_rec else {}

    # 저자 문자열
    authors = pubmed_rec.get("authors", []) if pubmed_rec else []
    author_str = ", ".join(authors[:6])
    if len(authors) > 6:
        author_str += f" et al."

    journal = pubmed_rec.get("journal", paper.get("journal", "")) if pubmed_rec else paper.get("journal", "")
    year    = pubmed_rec.get("year", "") if pubmed_rec else ""
    doi     = pubmed_rec.get("doi", "") if pubmed_rec else ""

    # 본문 텍스트 구성
    lines = [
        f"Disease: {disease_name}",
        f"Category: {chunk_type.replace('_', ' ').title()}",
        f"",
        f"Title: {title}",
    ]
    if author_str:
        lines.append(f"Authors: {author_str}")
    if journal or year:
        lines.append(f"Journal: {journal} ({year})")
    if doi:
        lines.append(f"DOI: {doi}")
    lines.append(f"PMID: {pmid}")

    if abstract_sections:
        # Structured abstract → 섹션별 출력
        lines.append("")
        lines.append("Abstract:")
        for label, content in abstract_sections.items():
            lines.append(f"  [{label}] {content}")
    elif abstract:
        lines.append("")
        lines.append(f"Abstract: {abstract}")
    else:
        lines.append("")
        lines.append("(Abstract not available)")

    mesh = pubmed_rec.get("mesh_terms", []) if pubmed_rec else []
    if mesh:
        lines.append("")
        lines.append(f"MeSH Terms: {'; '.join(mesh[:12])}")

    kw = pubmed_rec.get("keywords", []) if pubmed_rec else []
    if kw:
        lines.append(f"Keywords: {'; '.join(kw[:8])}")

    text = "\n".join(lines)

    return {
        "chunk_id":   _chunk_id(uid, chunk_type, pmid),
        "chunk_type": chunk_type,
        "disease_name": disease_name,
        "disease_uid":  uid,
        "concept_id":   disease.get("concept_id", ""),
        "text":         text.strip(),
        "metadata": {
            **_base_meta(disease),
            **_paper_meta(paper, pubmed_rec, chunk_type),
            "has_abstract": bool(abstract),
            "has_structured_abstract": bool(abstract_sections),
        },
    }


# ── 메인 빌더 ─────────────────────────────────────────────────

def build_chunks(
    medgen_records: list[dict],
    pubmed_lookup:  dict[str, dict],
) -> list[dict]:
    chunks: list[dict] = []
    seen_chunk_ids: set[str] = set()  # 중복 방지

    def _add(chunk: dict | None):
        if chunk is None:
            return
        cid = chunk["chunk_id"]
        if cid in seen_chunk_ids:
            return
        seen_chunk_ids.add(cid)
        chunks.append(chunk)

    for disease in medgen_records:
        uid = disease.get("medgen_uid", "")
        lit = disease.get("literature", {})

        # 1. 질환 개요
        _add(make_overview_chunk(disease))

        # 2. Professional guidelines
        gl = lit.get("professional_guidelines", {})
        for paper in gl.get("papers", []):
            _add(make_paper_chunk(disease, paper, "clinical_guideline", pubmed_lookup))

        # 3. Recent clinical studies (카테고리별)
        cs = lit.get("recent_clinical_studies", {})
        for cat_key, chunk_type in CLINICAL_CAT_MAP.items():
            cat_data = cs.get(cat_key, {})
            for paper in cat_data.get("papers", []):
                _add(make_paper_chunk(disease, paper, chunk_type, pubmed_lookup))

        # 4. Systematic reviews
        sr = lit.get("recent_systematic_reviews", {})
        for paper in sr.get("papers", []):
            _add(make_paper_chunk(disease, paper, "systematic_review", pubmed_lookup))

        # 5. ELink fallback
        for paper in lit.get("elink_related", []):
            _add(make_paper_chunk(disease, paper, "related_literature", pubmed_lookup))

    return chunks


# ── 통계 출력 ─────────────────────────────────────────────────

def print_stats(chunks: list[dict]) -> None:
    from collections import Counter
    type_counts = Counter(c["chunk_type"] for c in chunks)
    has_abstract = sum(1 for c in chunks if c.get("metadata", {}).get("has_abstract"))
    has_pdf = sum(
        1 for c in chunks
        if any(l.get("type") in ("pdf",) for l in c.get("metadata", {}).get("pdf_links", []))
    )
    diseases = set(c["disease_uid"] for c in chunks)

    print(f"\n{'='*55}")
    print(f"  총 chunk 수     : {len(chunks)}")
    print(f"  질환 수         : {len(diseases)}")
    print(f"  abstract 포함   : {has_abstract}")
    print(f"  PDF 링크 보유   : {has_pdf}")
    print(f"\n  chunk_type 분포:")
    for t, n in sorted(type_counts.items(), key=lambda x: -x[1]):
        print(f"    {t:<30} {n:>4}개")
    print(f"{'='*55}")


# ── CLI ───────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="MedGen + PubMed JSON → RAG chunk 빌더",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
예시:
  python chunk_builder.py \\
      --medgen medgen_results.json \\
      --pubmed pubmed_abstracts.json \\
      --output chunks.json
        """,
    )
    parser.add_argument("--medgen", required=True, help="medgen_results.json 경로")
    parser.add_argument("--pubmed", required=True, help="pubmed_abstracts.json 경로")
    parser.add_argument("--output", default="chunks.json", help="출력 파일명")
    args = parser.parse_args()

    print(f"▶ MedGen JSON 로드: {args.medgen}")
    with open(args.medgen, encoding="utf-8") as f:
        medgen_records: list[dict] = json.load(f)
    print(f"  {len(medgen_records)}개 질환 레코드")

    print(f"▶ PubMed JSON 로드: {args.pubmed}")
    with open(args.pubmed, encoding="utf-8") as f:
        pubmed_records: list[dict] = json.load(f)
    print(f"  {len(pubmed_records)}개 abstract 레코드")

    pubmed_lookup = build_pubmed_lookup(pubmed_records)
    print(f"  PMID lookup 완료: {len(pubmed_lookup)}개")

    print("\n▶ Chunk 생성 중...")
    chunks = build_chunks(medgen_records, pubmed_lookup)

    print_stats(chunks)

    with open(args.output, "w", encoding="utf-8") as f:
        json.dump(chunks, f, ensure_ascii=False, indent=2)
    print(f"\n▶ 저장 완료: {args.output}")

    # 샘플 미리보기
    print("\n▶ 샘플 chunk 미리보기:")
    for chunk in chunks[:3]:
        print(f"\n  chunk_id  : {chunk['chunk_id']}")
        print(f"  type      : {chunk['chunk_type']}")
        print(f"  disease   : {chunk['disease_name']}")
        print(f"  text[:200]:\n{chunk['text'][:200]}")
        print(f"  ...")


if __name__ == "__main__":
    main()