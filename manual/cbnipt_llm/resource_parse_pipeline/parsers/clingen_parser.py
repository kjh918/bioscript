"""
clingen_parser.py
─────────────────
ClinGen Dosage Sensitivity TSV parser.

Key behavior
------------
1. Parses columns by normalized header name, not by fixed column index.
2. Supports both Gene and Region curation TSV files.
3. Preserves '#ISCA ID' header rows instead of treating them as comments.
4. Falls back to the historical 25-column layout only when a usable header
   cannot be found.
5. Splits PMID fields safely and prevents dates/MONDO IDs from being stored
   as PMIDs.

Expected input examples
-----------------------
  ClinGen_gene_curation_list_GRCh38.tsv
  ClinGen_region_curation_list_GRCh38.tsv
"""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any, Iterable


# -----------------------------------------------------------------------------
# Score labels
# -----------------------------------------------------------------------------

DOSAGE_LABEL = {
    "0": "No evidence available",
    "1": "Little evidence for dosage pathogenicity",
    "2": "Emerging evidence for dosage pathogenicity",
    "3": "Sufficient evidence for dosage pathogenicity",
    "30": "Gene associated with autosomal recessive phenotype",
    "40": "Dosage sensitivity unlikely",
}


# -----------------------------------------------------------------------------
# Header aliases
# -----------------------------------------------------------------------------

# All keys and aliases are normalized by _normalize_header().
HEADER_ALIASES: dict[str, tuple[str, ...]] = {
    # Identity
    "gene_symbol": (
        "gene symbol",
        "gene",
        "symbol",
    ),
    "gene_id": (
        "gene id",
        "entrez gene id",
        "ncbi gene id",
        "geneid",
    ),
    "isca_id": (
        "isca id",
        "isca region id",
        "region id",
    ),
    "region_name": (
        "isca region name",
        "region name",
        "name",
    ),
    "cytoband": (
        "cytoband",
        "cyto band",
        "cytogenetic location",
    ),
    "genomic_location": (
        "genomic location grch38",
        "grch38 genomic location",
        "genomic location",
        "location grch38",
        "grch38 location",
    ),

    # Haploinsufficiency
    "hi_score": (
        "haploinsufficiency score",
        "haploinsufficiency",
        "hi score",
    ),
    "hi_score_desc": (
        "haploinsufficiency score description",
        "hi score description",
        "haploinsufficiency description",
    ),
    "hi_disease_id": (
        "haploinsufficiency disease id",
        "hi disease id",
        "haploinsufficiency phenotype id",
        "hi phenotype id",
    ),
    "hi_disease_name": (
        "haploinsufficiency disease name",
        "hi disease name",
        "haploinsufficiency phenotype name",
        "hi phenotype name",
    ),
    "hi_description": (
        "haploinsufficiency evidence description",
        "hi evidence description",
        "haploinsufficiency evidence",
        "hi evidence",
    ),
    "hi_pmid1": (
        "haploinsufficiency pmid1",
        "hi pmid1",
        "haploinsufficiency pmid 1",
        "hi pmid 1",
    ),
    "hi_pmid2": (
        "haploinsufficiency pmid2",
        "hi pmid2",
        "haploinsufficiency pmid 2",
        "hi pmid 2",
    ),
    "hi_pmid3": (
        "haploinsufficiency pmid3",
        "hi pmid3",
        "haploinsufficiency pmid 3",
        "hi pmid 3",
    ),
    "hi_pmids": (
        "haploinsufficiency pmids",
        "hi pmids",
        "haploinsufficiency PMID",
        "hi PMID",
    ),
    "hi_notes": (
        "haploinsufficiency notes",
        "hi notes",
    ),

    # Triplosensitivity
    "ts_score": (
        "triplosensitivity score",
        "triplosensitivity",
        "ts score",
    ),
    "ts_score_desc": (
        "triplosensitivity score description",
        "ts score description",
        "triplosensitivity description",
    ),
    "ts_disease_id": (
        "triplosensitivity disease id",
        "ts disease id",
        "triplosensitivity phenotype id",
        "ts phenotype id",
    ),
    "ts_disease_name": (
        "triplosensitivity disease name",
        "ts disease name",
        "triplosensitivity phenotype name",
        "ts phenotype name",
    ),
    "ts_description": (
        "triplosensitivity evidence description",
        "ts evidence description",
        "triplosensitivity evidence",
        "ts evidence",
    ),
    "ts_pmid1": (
        "triplosensitivity pmid1",
        "ts pmid1",
        "triplosensitivity pmid 1",
        "ts pmid 1",
    ),
    "ts_pmid2": (
        "triplosensitivity pmid2",
        "ts pmid2",
        "triplosensitivity pmid 2",
        "ts pmid 2",
    ),
    "ts_pmid3": (
        "triplosensitivity pmid3",
        "ts pmid3",
        "triplosensitivity pmid 3",
        "ts pmid 3",
    ),
    "ts_pmids": (
        "triplosensitivity pmids",
        "ts pmids",
        "triplosensitivity PMID",
        "ts PMID",
    ),
    "ts_notes": (
        "triplosensitivity notes",
        "ts notes",
    ),

    # Metadata
    "date_evaluated": (
        "date last evaluated",
        "last evaluated",
        "date evaluated",
        "evaluation date",
    ),
    "loss_omim_id": (
        "loss phenotype omim id",
        "haploinsufficiency phenotype omim id",
        "hi phenotype omim id",
        "loss omim id",
    ),
    "ts_omim_id": (
        "ts phenotype omim id",
        "triplosensitivity phenotype omim id",
        "gain phenotype omim id",
        "gain omim id",
    ),
}


# Historical 25-column layout. Used only when no header can be resolved.
LEGACY_INDEX = {
    "primary_0": 0,
    "primary_1": 1,
    "cytoband": 2,
    "genomic_location": 3,
    "hi_score": 4,
    "hi_score_desc": 5,
    "hi_disease_id": 6,
    "hi_disease_name": 7,
    "hi_description": 8,
    "hi_pmid1": 9,
    "hi_pmid2": 10,
    "hi_pmid3": 11,
    "hi_notes": 12,
    "ts_score": 13,
    "ts_score_desc": 14,
    "ts_disease_id": 15,
    "ts_disease_name": 16,
    "ts_description": 17,
    "ts_pmid1": 18,
    "ts_pmid2": 19,
    "ts_pmid3": 20,
    "ts_notes": 21,
    "date_evaluated": 22,
    "loss_omim_id": 23,
    "ts_omim_id": 24,
}


# -----------------------------------------------------------------------------
# General helpers
# -----------------------------------------------------------------------------

def _clean(value: Any) -> str:
    if value is None:
        return ""
    return str(value).strip().strip("\ufeff")


def _normalize_header(value: str) -> str:
    """Normalize a header for robust alias matching."""
    value = _clean(value).lstrip("#").lower()
    value = value.replace("_", " ").replace("-", " ")
    value = re.sub(r"[()\[\]{}:/]", " ", value)
    value = re.sub(r"\s+", " ", value).strip()
    return value


def _canonical_alias_map() -> dict[str, str]:
    mapping: dict[str, str] = {}
    for canonical, aliases in HEADER_ALIASES.items():
        mapping[_normalize_header(canonical)] = canonical
        for alias in aliases:
            mapping[_normalize_header(alias)] = canonical
    return mapping


ALIAS_TO_CANONICAL = _canonical_alias_map()


def _looks_like_header(row: list[str]) -> bool:
    normalized = {_normalize_header(v) for v in row if _clean(v)}
    recognized = sum(1 for v in normalized if v in ALIAS_TO_CANONICAL)

    first = _normalize_header(row[0]) if row else ""
    identity_header = first in {
        "gene symbol",
        "gene",
        "isca id",
        "isca region id",
        "region id",
    }
    return identity_header or recognized >= 3


def _header_index(header: list[str]) -> dict[str, int]:
    """Create canonical field -> column index mapping from the real header."""
    index: dict[str, int] = {}
    for idx, raw_name in enumerate(header):
        normalized = _normalize_header(raw_name)
        canonical = ALIAS_TO_CANONICAL.get(normalized)
        if canonical and canonical not in index:
            index[canonical] = idx
    return index


def _safe_index(row: list[str], idx: int | None) -> str:
    if idx is None or idx < 0 or idx >= len(row):
        return ""
    return _clean(row[idx])


def _field(
    row: list[str],
    index: dict[str, int],
    name: str,
    *,
    legacy_index: int | None = None,
) -> str:
    if name in index:
        return _safe_index(row, index[name])
    if legacy_index is not None:
        return _safe_index(row, legacy_index)
    return ""


def _parse_loc(loc_str: str) -> tuple[str, int, int]:
    """
    Parse forms such as:
      chr22:19756703-19783593
      22:19756703-19783593
      chrX:1,000-2,000
    """
    value = _clean(loc_str).replace(",", "")
    match = re.search(
        r"(?:chr)?(?P<chrom>\d{1,2}|X|Y|M|MT)\s*:\s*"
        r"(?P<start>\d+)\s*[-–—]\s*(?P<end>\d+)",
        value,
        flags=re.IGNORECASE,
    )
    if not match:
        return "", 0, 0

    chrom_token = match.group("chrom").upper()
    if chrom_token == "MT":
        chrom_token = "M"

    start = int(match.group("start"))
    end = int(match.group("end"))
    if start > end:
        start, end = end, start

    return f"chr{chrom_token}", start, end


def _is_valid_pmid(value: str) -> bool:
    """A PMID is numeric. Dates and MONDO/OMIM IDs must not enter PMID lists."""
    return bool(re.fullmatch(r"\d{4,10}", _clean(value)))


def _split_pmids(*values: str) -> list[str]:
    output: list[str] = []
    seen: set[str] = set()

    for raw in values:
        for token in re.split(r"[;,|\s]+", _clean(raw)):
            token = token.strip()
            if not _is_valid_pmid(token):
                continue
            if token not in seen:
                seen.add(token)
                output.append(token)

    return output


def _score_label(score: str, original_description: str) -> str:
    score = _clean(score)
    original_description = _clean(original_description)
    return DOSAGE_LABEL.get(score, original_description)


def _read_tsv(path: str) -> tuple[list[str], list[list[str]], bool]:
    """
    Read ClinGen TSV while retaining a '#ISCA ID' header.

    Returns
    -------
    header, rows, has_resolved_header
    """
    header: list[str] = []
    rows: list[list[str]] = []
    has_resolved_header = False

    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")

        for raw_row in reader:
            row = [_clean(v) for v in raw_row]
            if not row or not any(row):
                continue

            first = row[0]

            # Keep a possible '#ISCA ID' header, but skip metadata comments.
            if first.startswith("#"):
                if _looks_like_header(row):
                    header = row
                    has_resolved_header = True
                elif re.match(r"^#ISCA-", first, flags=re.IGNORECASE):
                    rows.append(row)
                continue

            if not has_resolved_header and _looks_like_header(row):
                header = row
                has_resolved_header = True
                continue

            rows.append(row)

    return header, rows, has_resolved_header


def _warn_schema(path: str, header: list[str], index: dict[str, int], kind: str) -> None:
    identity_required = "gene_symbol" if kind == "gene" else "isca_id"
    required = {identity_required, "genomic_location", "hi_score", "ts_score"}
    missing = sorted(required - set(index))

    if missing:
        print(
            f"  [ClinGen {kind.title()}] 경고: 헤더에서 필수 컬럼 일부를 "
            f"찾지 못했습니다: {missing}"
        )
        print(f"  [ClinGen {kind.title()}] 실제 헤더: {header}")
        print(
            "  [ClinGen] 헤더 기반 매핑이 불완전하여 해당 필드만 "
            "legacy 25-column index로 보완합니다."
        )


def _common_record(
    row: list[str],
    index: dict[str, int],
    *,
    allow_legacy: bool,
) -> dict[str, Any]:
    legacy = LEGACY_INDEX if allow_legacy else {}

    cytoband = _field(
        row,
        index,
        "cytoband",
        legacy_index=legacy.get("cytoband"),
    )
    genomic_location = _field(
        row,
        index,
        "genomic_location",
        legacy_index=legacy.get("genomic_location"),
    )
    chrom, start, end = _parse_loc(genomic_location)

    hi_score = _field(
        row,
        index,
        "hi_score",
        legacy_index=legacy.get("hi_score"),
    )
    hi_score_desc = _field(
        row,
        index,
        "hi_score_desc",
        legacy_index=legacy.get("hi_score_desc"),
    )
    ts_score = _field(
        row,
        index,
        "ts_score",
        legacy_index=legacy.get("ts_score"),
    )
    ts_score_desc = _field(
        row,
        index,
        "ts_score_desc",
        legacy_index=legacy.get("ts_score_desc"),
    )

    hi_pmids = _split_pmids(
        _field(row, index, "hi_pmids"),
        _field(row, index, "hi_pmid1", legacy_index=legacy.get("hi_pmid1")),
        _field(row, index, "hi_pmid2", legacy_index=legacy.get("hi_pmid2")),
        _field(row, index, "hi_pmid3", legacy_index=legacy.get("hi_pmid3")),
    )
    ts_pmids = _split_pmids(
        _field(row, index, "ts_pmids"),
        _field(row, index, "ts_pmid1", legacy_index=legacy.get("ts_pmid1")),
        _field(row, index, "ts_pmid2", legacy_index=legacy.get("ts_pmid2")),
        _field(row, index, "ts_pmid3", legacy_index=legacy.get("ts_pmid3")),
    )

    return {
        "cytoband": cytoband,
        "genomic_location": genomic_location,
        "chrom": chrom,
        "start": start,
        "end": end,

        "hi_score": hi_score,
        "hi_score_label": _score_label(hi_score, hi_score_desc),
        "hi_score_desc": hi_score_desc,
        "hi_disease_id": _field(
            row,
            index,
            "hi_disease_id",
            legacy_index=legacy.get("hi_disease_id"),
        ),
        "hi_disease_name": _field(
            row,
            index,
            "hi_disease_name",
            legacy_index=legacy.get("hi_disease_name"),
        ),
        "hi_description": _field(
            row,
            index,
            "hi_description",
            legacy_index=legacy.get("hi_description"),
        ),
        "hi_pmids": hi_pmids,
        "hi_notes": _field(
            row,
            index,
            "hi_notes",
            legacy_index=legacy.get("hi_notes"),
        ),

        "ts_score": ts_score,
        "ts_score_label": _score_label(ts_score, ts_score_desc),
        "ts_score_desc": ts_score_desc,
        "ts_disease_id": _field(
            row,
            index,
            "ts_disease_id",
            legacy_index=legacy.get("ts_disease_id"),
        ),
        "ts_disease_name": _field(
            row,
            index,
            "ts_disease_name",
            legacy_index=legacy.get("ts_disease_name"),
        ),
        "ts_description": _field(
            row,
            index,
            "ts_description",
            legacy_index=legacy.get("ts_description"),
        ),
        "ts_pmids": ts_pmids,
        "ts_notes": _field(
            row,
            index,
            "ts_notes",
            legacy_index=legacy.get("ts_notes"),
        ),

        "date_evaluated": _field(
            row,
            index,
            "date_evaluated",
            legacy_index=legacy.get("date_evaluated"),
        ),
        "loss_omim_id": _field(
            row,
            index,
            "loss_omim_id",
            legacy_index=legacy.get("loss_omim_id"),
        ),
        "ts_omim_id": _field(
            row,
            index,
            "ts_omim_id",
            legacy_index=legacy.get("ts_omim_id"),
        ),
    }


# -----------------------------------------------------------------------------
# Gene parser
# -----------------------------------------------------------------------------

def parse_clingen_gene(path: str) -> dict[str, dict[str, Any]]:
    """Parse a ClinGen gene dosage curation TSV into a gene-symbol lookup."""
    if not path or not Path(path).is_file():
        print(f"  [ClinGen Gene] 파일 없음: {path}")
        return {}

    header, rows, has_header = _read_tsv(path)
    index = _header_index(header) if has_header else {}
    if has_header:
        _warn_schema(path, header, index, "gene")
    else:
        print(
            "  [ClinGen Gene] 경고: 헤더를 찾지 못해 legacy 25-column "
            "index를 사용합니다."
        )

    lookup: dict[str, dict[str, Any]] = {}

    for line_number, row in enumerate(rows, start=2 if has_header else 1):
        symbol = _field(
            row,
            index,
            "gene_symbol",
            legacy_index=(LEGACY_INDEX["primary_0"] if not has_header else None),
        )
        if not symbol or symbol.startswith("#"):
            continue

        common = _common_record(row, index, allow_legacy=not has_header)
        if not common["chrom"]:
            print(
                f"  [ClinGen Gene] 좌표 파싱 실패: line={line_number}, "
                f"gene={symbol}, location={common['genomic_location']!r}"
            )

        lookup[symbol] = {
            "gene_symbol": symbol,
            "ncbi_gene_id": _field(
                row,
                index,
                "gene_id",
                legacy_index=(LEGACY_INDEX["primary_1"] if not has_header else None),
            ),
            **common,
        }

    print(f"  [ClinGen Gene] {len(lookup)}개 유전자 로드")
    return lookup


# -----------------------------------------------------------------------------
# Region parser
# -----------------------------------------------------------------------------

def parse_clingen_region(path: str) -> list[dict[str, Any]]:
    """Parse a ClinGen region dosage curation TSV."""
    if not path or not Path(path).is_file():
        print(f"  [ClinGen Region] 파일 없음: {path}")
        return []

    header, rows, has_header = _read_tsv(path)
    index = _header_index(header) if has_header else {}
    if has_header:
        _warn_schema(path, header, index, "region")
    else:
        print(
            "  [ClinGen Region] 경고: 헤더를 찾지 못해 legacy 25-column "
            "index를 사용합니다."
        )

    regions: list[dict[str, Any]] = []

    for line_number, row in enumerate(rows, start=2 if has_header else 1):
        isca_id = _field(
            row,
            index,
            "isca_id",
            legacy_index=(LEGACY_INDEX["primary_0"] if not has_header else None),
        ).lstrip("#").strip()
        if not isca_id:
            continue

        common = _common_record(row, index, allow_legacy=not has_header)
        if not common["chrom"]:
            print(
                f"  [ClinGen Region] 좌표 파싱 실패로 제외: "
                f"line={line_number}, isca_id={isca_id}, "
                f"location={common['genomic_location']!r}"
            )
            continue

        regions.append({
            "isca_id": isca_id,
            "region_name": _field(
                row,
                index,
                "region_name",
                legacy_index=(LEGACY_INDEX["primary_1"] if not has_header else None),
            ),
            **common,
        })

    print(f"  [ClinGen Region] {len(regions)}개 영역 로드")
    return regions


# -----------------------------------------------------------------------------
# Overlap
# -----------------------------------------------------------------------------

def get_region_overlaps(
    chrom: str,
    start: int,
    end: int,
    regions: Iterable[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Return ClinGen regions overlapping a genomic interval."""
    if not chrom or start <= 0 or end <= 0:
        return []
    if start > end:
        start, end = end, start

    normalized_chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
    return [
        region
        for region in regions
        if region.get("chrom") == normalized_chrom
        and int(region.get("start", 0)) <= end
        and int(region.get("end", 0)) >= start
    ]


# -----------------------------------------------------------------------------
# Disease / cytoband matching helpers
# -----------------------------------------------------------------------------

def _normalize_identifier(value: str) -> str:
    """Normalize MONDO/OMIM/ORPHA identifiers for exact token matching."""
    value = _clean(value).upper().replace(" ", "")
    value = value.replace("ORPHANET:", "ORPHA:")
    if value.isdigit():
        return value
    return value


def _identifier_tokens(value: str) -> set[str]:
    """Extract identifiers without confusing PMIDs with disease identifiers."""
    text = _clean(value).upper()
    tokens: set[str] = set()
    for prefix, number in re.findall(
        r"\b(MONDO|OMIM|ORPHA|ORPHANET)\s*[:#]?\s*(\d+)\b", text
    ):
        normalized_prefix = "ORPHA" if prefix == "ORPHANET" else prefix
        tokens.add(f"{normalized_prefix}:{number}")
        tokens.add(number)
    return tokens


def _normalize_disease_text(value: str) -> str:
    value = _clean(value).lower()
    value = re.sub(r"\bsyndrome\b", " ", value)
    value = re.sub(r"[^a-z0-9]+", " ", value)
    return re.sub(r"\s+", " ", value).strip()


def _disease_text_matches(query: str, candidate: str) -> bool:
    q = _normalize_disease_text(query)
    c = _normalize_disease_text(candidate)
    if not q or not c:
        return False
    if q == c or q in c or c in q:
        return True
    q_tokens = set(q.split())
    c_tokens = set(c.split())
    if not q_tokens or not c_tokens:
        return False
    # Require meaningful overlap, avoiding generic one-word matches.
    overlap = q_tokens & c_tokens
    return len(overlap) >= 2 and len(overlap) / len(q_tokens) >= 0.6


def _record_matches_disease(
    record: dict[str, Any],
    syndrome_name: str = "",
    mondo_id: str = "",
    omim_id: str = "",
    orphanet_id: str = "",
) -> bool:
    candidate_ids = set()
    for key in (
        "hi_disease_id", "ts_disease_id", "loss_omim_id", "ts_omim_id"
    ):
        candidate_ids |= _identifier_tokens(str(record.get(key, "")))

    requested_ids: set[str] = set()
    for prefix, raw in (
        ("MONDO", mondo_id), ("OMIM", omim_id), ("ORPHA", orphanet_id)
    ):
        raw = _clean(raw)
        if not raw:
            continue
        m = re.search(r"(\d+)", raw)
        if m:
            requested_ids.add(f"{prefix}:{m.group(1)}")
            requested_ids.add(m.group(1))

    if requested_ids and candidate_ids & requested_ids:
        return True

    if syndrome_name:
        return any(
            _disease_text_matches(syndrome_name, str(record.get(key, "")))
            for key in ("hi_disease_name", "ts_disease_name", "region_name")
        )
    return False


def normalize_cytoband(value: str) -> str:
    """Normalize chr17p13.3 / 17p13.3 / 17P13.3 to 17P13.3."""
    value = _clean(value).upper().replace("CHR", "")
    value = value.replace("–", "-").replace("—", "-")
    value = re.sub(r"\s+", "", value)
    return value


def cytoband_matches(record_band: str, target_band: str) -> bool:
    """Match whole chromosome, arm, band, or band range conservatively."""
    record = normalize_cytoband(record_band)
    target = normalize_cytoband(target_band)
    if not record or not target:
        return False

    if re.fullmatch(r"(?:\d{1,2}|X|Y)", target):
        return record.startswith(target + "P") or record.startswith(target + "Q")

    if "-" in target:
        left, right = target.split("-", 1)
        chrom_match = re.match(r"^(\d{1,2}|X|Y)", left)
        if chrom_match and re.match(r"^[PQ]", right):
            right = chrom_match.group(1) + right
        return (
            record.startswith(left)
            or record.startswith(right)
            or left.startswith(record)
            or right.startswith(record)
        )

    return record.startswith(target) or target.startswith(record)


def find_clingen_genes(
    genes: dict[str, dict[str, Any]],
    syndrome_name: str = "",
    mondo_id: str = "",
    omim_id: str = "",
    orphanet_id: str = "",
    cytobands: Iterable[str] = (),
    event: str = "unknown",
) -> list[dict[str, Any]]:
    """Find curated genes by disease identity and/or target cytoband."""
    results: dict[str, dict[str, Any]] = {}
    bands = [b for b in cytobands if _clean(b)]

    for symbol, record in genes.items():
        disease_match = _record_matches_disease(
            record, syndrome_name, mondo_id, omim_id, orphanet_id
        )
        band_match = any(
            cytoband_matches(str(record.get("cytoband", "")), band)
            for band in bands
        )

        score_ok = True
        if event == "loss":
            score_ok = str(record.get("hi_score", "")) in {"1", "2", "3"}
        elif event == "gain":
            score_ok = str(record.get("ts_score", "")) in {"1", "2", "3"}

        if disease_match or (band_match and score_ok):
            item = dict(record)
            item["disease_match"] = disease_match
            item["cytoband_match"] = band_match
            results[symbol] = item

    return list(results.values())


def find_clingen_regions(
    regions: Iterable[dict[str, Any]],
    syndrome_name: str = "",
    mondo_id: str = "",
    omim_id: str = "",
    orphanet_id: str = "",
    cytobands: Iterable[str] = (),
    event: str = "unknown",
) -> list[dict[str, Any]]:
    """Find curated regions by disease identity and/or cytoband."""
    output: list[dict[str, Any]] = []
    seen: set[str] = set()
    bands = [b for b in cytobands if _clean(b)]

    for record in regions:
        disease_match = _record_matches_disease(
            record, syndrome_name, mondo_id, omim_id, orphanet_id
        )
        band_match = any(
            cytoband_matches(str(record.get("cytoband", "")), band)
            for band in bands
        )

        score_ok = True
        if event == "loss":
            score_ok = str(record.get("hi_score", "")) in {"1", "2", "3"}
        elif event == "gain":
            score_ok = str(record.get("ts_score", "")) in {"1", "2", "3"}

        if disease_match or (band_match and score_ok):
            key = str(record.get("isca_id", "")) or (
                f"{record.get('chrom')}:{record.get('start')}-{record.get('end')}"
            )
            if key in seen:
                continue
            seen.add(key)
            item = dict(record)
            item["disease_match"] = disease_match
            item["cytoband_match"] = band_match
            output.append(item)

    return output


# -----------------------------------------------------------------------------
# Optional validation utility
# -----------------------------------------------------------------------------

def validate_clingen_file(path: str, kind: str) -> None:
    """
    Print parsed header mapping for troubleshooting.

    kind must be 'gene' or 'region'.
    """
    if kind not in {"gene", "region"}:
        raise ValueError("kind must be 'gene' or 'region'")
    if not Path(path).is_file():
        raise FileNotFoundError(path)

    header, rows, has_header = _read_tsv(path)
    index = _header_index(header) if has_header else {}

    print(f"[ClinGen validate] path={path}")
    print(f"[ClinGen validate] kind={kind}")
    print(f"[ClinGen validate] header_detected={has_header}")
    print(f"[ClinGen validate] n_columns={len(header) if header else 'unknown'}")
    print(f"[ClinGen validate] n_data_rows={len(rows)}")
    print("[ClinGen validate] canonical mapping:")
    for field_name, column_index in sorted(index.items(), key=lambda item: item[1]):
        print(
            f"  {column_index:>3}: {field_name:<22} "
            f"<- {header[column_index]!r}"
        )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Validate or parse ClinGen dosage sensitivity TSV files."
    )
    parser.add_argument("path", help="ClinGen TSV path")
    parser.add_argument("--kind", choices=["gene", "region"], required=True)
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Only print detected header mapping",
    )
    args = parser.parse_args()

    validate_clingen_file(args.path, args.kind)
    if not args.validate_only:
        if args.kind == "gene":
            parse_clingen_gene(args.path)
        else:
            parse_clingen_region(args.path)