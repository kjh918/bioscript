#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path
from typing import Any

SPACE_RE = re.compile(r"\s+")
SENTENCE_RE = re.compile(r"(?<=[.!?])\s+(?=[A-Z0-9\[(])")


def clean(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, list):
        return "; ".join(x for x in (clean(v) for v in value) if x)
    if isinstance(value, dict):
        return "; ".join(f"{k}: {clean(v)}" for k, v in value.items() if clean(v))
    return SPACE_RE.sub(" ", str(value)).strip()


def read_json(path: Path) -> Any:
    raw = path.read_text(encoding="utf-8-sig").strip()
    if not raw:
        raise ValueError(f"Empty input file: {path}")
    try:
        return json.loads(raw)
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid JSON: {path}\n{exc}") from exc


def top_records(data: Any) -> list[dict[str, Any]]:
    if isinstance(data, list):
        return [x for x in data if isinstance(x, dict)]
    if isinstance(data, dict):
        # Common wrappers, then fallback to single record.
        for key in ("records", "results", "items", "articles", "data"):
            value = data.get(key)
            if isinstance(value, list):
                return [x for x in value if isinstance(x, dict)]
        return [data]
    return []


def detect_source(records: list[dict[str, Any]]) -> str:
    if not records:
        raise ValueError("No object records found")
    sample = records[0]
    keys = set(sample)
    pubmed_score = sum(k in keys for k in ("pmid", "abstract", "abstract_sections", "rag_text", "mesh_terms"))
    medgen_score = sum(k in keys for k in ("medgen_uid", "concept_id", "disease_name", "definition", "literature"))
    if pubmed_score > medgen_score and pubmed_score >= 2:
        return "pubmed"
    if medgen_score > pubmed_score and medgen_score >= 2:
        return "medgen"
    raise ValueError(f"Could not detect source type from keys: {sorted(keys)}")


def split_text(text: str, max_chars: int, overlap_chars: int) -> list[str]:
    text = clean(text)
    if not text:
        return []
    if len(text) <= max_chars:
        return [text]

    sentences = [s.strip() for s in SENTENCE_RE.split(text) if s.strip()]
    chunks: list[str] = []
    current = ""

    for sentence in sentences:
        if len(sentence) > max_chars:
            if current:
                chunks.append(current)
                current = ""
            step = max(1, max_chars - overlap_chars)
            for start in range(0, len(sentence), step):
                piece = sentence[start:start + max_chars].strip()
                if piece:
                    chunks.append(piece)
            continue

        candidate = f"{current} {sentence}".strip() if current else sentence
        if len(candidate) <= max_chars:
            current = candidate
            continue

        chunks.append(current)
        overlap = current[-overlap_chars:].strip() if overlap_chars > 0 else ""
        current = f"{overlap} {sentence}".strip()

    if current:
        chunks.append(current)

    return chunks


def make_id(*parts: Any) -> str:
    raw = "|".join(clean(p) for p in parts)
    return hashlib.sha1(raw.encode("utf-8")).hexdigest()[:16]


def year_value(value: Any) -> int | None:
    text = clean(value)
    match = re.search(r"\b(19|20)\d{2}\b", text)
    return int(match.group(0)) if match else None


def pubmed_chunks(record: dict[str, Any], max_chars: int, overlap_chars: int) -> list[dict[str, Any]]:
    pmid = clean(record.get("pmid"))
    if not pmid:
        return []
    title = clean(record.get("title")) or f"PubMed {pmid}"
    sections = record.get("abstract_sections")

    content_sections: list[tuple[str, str]] = []
    if isinstance(sections, dict) and any(clean(v) for v in sections.values()):
        content_sections = [(clean(k), clean(v)) for k, v in sections.items() if clean(v)]
    else:
        abstract = clean(record.get("abstract"))
        rag_text = clean(record.get("rag_text"))
        if abstract:
            content_sections = [("ABSTRACT", abstract)]
        elif rag_text:
            content_sections = [("RAG_TEXT", rag_text)]
        else:
            # Metadata-only papers still get one searchable chunk.
            metadata_text = "\n".join(x for x in [
                f"Title: {title}",
                f"Journal: {clean(record.get('journal'))}" if clean(record.get('journal')) else "",
                f"MeSH terms: {clean(record.get('mesh_terms'))}" if clean(record.get('mesh_terms')) else "",
                f"Keywords: {clean(record.get('keywords'))}" if clean(record.get('keywords')) else "",
            ] if x)
            if metadata_text:
                content_sections = [("METADATA", metadata_text)]

    out: list[dict[str, Any]] = []
    for section, body in content_sections:
        pieces = split_text(body, max_chars, overlap_chars)
        for idx, piece in enumerate(pieces, 1):
            text = f"Title: {title}\nPMID: {pmid}\nSection: {section}\n\n{piece}"
            chunk_id = f"pubmed:{pmid}:{section.lower().replace(' ', '_')}:{idx:03d}:{make_id(text)}"
            out.append({
                "chunk_id": chunk_id,
                "document_id": f"pubmed:{pmid}",
                "source_type": "pubmed",
                "title": title,
                "text": text,
                "metadata": {
                    "pmid": pmid,
                    "doi": clean(record.get("doi")),
                    "pmc_id": clean(record.get("pmc_id")),
                    "journal": clean(record.get("journal")),
                    "year": year_value(record.get("year")),
                    "pub_date": clean(record.get("pub_date")),
                    "publication_types": record.get("publication_types") or [],
                    "authors": record.get("authors") or [],
                    "mesh_terms": record.get("mesh_terms") or [],
                    "keywords": record.get("keywords") or [],
                    "section": section,
                    "chunk_index": idx,
                    "chunk_count": len(pieces),
                    "full_text_available": bool(record.get("full_text_available")),
                },
                "source": {
                    "source_id": pmid,
                    "url": clean(record.get("pubmed_url")) or f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                },
            })
    return out


def collect_linked_pmids(node: Any) -> list[str]:
    found: set[str] = set()
    def walk(x: Any) -> None:
        if isinstance(x, dict):
            pmid = clean(x.get("pmid"))
            if pmid:
                found.add(pmid)
            for value in x.values():
                walk(value)
        elif isinstance(x, list):
            for value in x:
                walk(value)
    walk(node)
    return sorted(found)


def collect_curated_resources(node: Any) -> list[dict[str, str]]:
    result: list[dict[str, str]] = []
    def walk(x: Any) -> None:
        if isinstance(x, dict):
            title, url = clean(x.get("title")), clean(x.get("url"))
            if title and url and "pmid" not in x:
                result.append({"title": title, "url": url})
            for value in x.values():
                walk(value)
        elif isinstance(x, list):
            for value in x:
                walk(value)
    walk(node)
    unique = {(r["title"], r["url"]): r for r in result}
    return list(unique.values())


def medgen_chunks(record: dict[str, Any], max_chars: int, overlap_chars: int) -> list[dict[str, Any]]:
    uid = clean(record.get("medgen_uid"))
    concept_id = clean(record.get("concept_id"))
    if not uid and not concept_id:
        return []
    primary_id = concept_id or uid
    title = clean(record.get("disease_name")) or f"MedGen {primary_id}"
    definition = clean(record.get("definition"))
    synonyms = record.get("synonyms") or []
    xrefs = record.get("cross_references") or {}
    linked_pmids = collect_linked_pmids(record.get("literature") or {})
    curated = collect_curated_resources(record.get("literature") or {})

    out: list[dict[str, Any]] = []
    if definition:
        pieces = split_text(definition, max_chars, overlap_chars)
        for idx, piece in enumerate(pieces, 1):
            text = (
                f"Disease: {title}\nMedGen UID: {uid}\nConcept ID: {concept_id}\n"
                f"Disease type: {clean(record.get('disease_type'))}\n"
                f"Synonyms: {clean(synonyms)}\n\nDefinition: {piece}"
            )
            out.append({
                "chunk_id": f"medgen:{primary_id}:definition:{idx:03d}:{make_id(text)}",
                "document_id": f"medgen:{primary_id}",
                "source_type": "medgen",
                "title": title,
                "text": text,
                "metadata": {
                    "medgen_uid": uid,
                    "concept_id": concept_id,
                    "disease_name": title,
                    "disease_type": clean(record.get("disease_type")),
                    "synonyms": synonyms,
                    "cross_references": xrefs,
                    "linked_pmids": linked_pmids,
                    "section": "definition",
                    "chunk_index": idx,
                    "chunk_count": len(pieces),
                },
                "source": {"source_id": uid or concept_id, "url": clean(record.get("url"))},
            })

    if linked_pmids or curated:
        curated_text = "; ".join(f"{x['title']} ({x['url']})" for x in curated)
        text = "\n".join(x for x in [
            f"Disease: {title}",
            f"Concept ID: {concept_id}",
            f"Related PubMed IDs: {'; '.join(linked_pmids)}" if linked_pmids else "",
            f"Curated resources: {curated_text}" if curated_text else "",
        ] if x)
        out.append({
            "chunk_id": f"medgen:{primary_id}:literature:001:{make_id(text)}",
            "document_id": f"medgen:{primary_id}",
            "source_type": "medgen",
            "title": title,
            "text": text,
            "metadata": {
                "medgen_uid": uid,
                "concept_id": concept_id,
                "disease_name": title,
                "disease_type": clean(record.get("disease_type")),
                "synonyms": synonyms,
                "cross_references": xrefs,
                "linked_pmids": linked_pmids,
                "curated_resources": curated,
                "section": "literature_links",
                "chunk_index": 1,
                "chunk_count": 1,
            },
            "source": {"source_id": uid or concept_id, "url": clean(record.get("url"))},
        })
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description="Auto-detect PubMed and MedGen JSON files and create RAG JSONL chunks")
    parser.add_argument("input1", type=Path)
    parser.add_argument("input2", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-chars", type=int, default=2400)
    parser.add_argument("--overlap-chars", type=int, default=250)
    args = parser.parse_args()

    loaded: dict[str, list[dict[str, Any]]] = {}
    for path in (args.input1, args.input2):
        records = top_records(read_json(path))
        source = detect_source(records)
        if source in loaded:
            raise ValueError(f"Two {source} files were provided; expected one PubMed and one MedGen file")
        loaded[source] = records
        print(f"Detected {source}: {path} ({len(records)} records)")

    if set(loaded) != {"pubmed", "medgen"}:
        raise ValueError(f"Expected PubMed and MedGen inputs, detected: {sorted(loaded)}")

    chunks: list[dict[str, Any]] = []
    pub_count = 0
    med_count = 0
    for record in loaded["pubmed"]:
        built = pubmed_chunks(record, args.max_chars, args.overlap_chars)
        chunks.extend(built)
        pub_count += len(built)
    for record in loaded["medgen"]:
        built = medgen_chunks(record, args.max_chars, args.overlap_chars)
        chunks.extend(built)
        med_count += len(built)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for chunk in chunks:
            handle.write(json.dumps(chunk, ensure_ascii=False) + "\n")

    print(f"PubMed chunks written: {pub_count}")
    print(f"MedGen chunks written: {med_count}")
    print(f"Total chunks written: {len(chunks)}")
    print(f"Output: {args.output.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())