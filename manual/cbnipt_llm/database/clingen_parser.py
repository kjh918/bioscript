"""
clingen_parser.py
─────────────────
ClinGen Dosage Sensitivity TSV 파서 (인덱스 기반, 헤더 유무 무관)

파일:
  ClinGen_gene_curation_list_GRCh38.tsv
  ClinGen_region_curation_list_GRCh38.tsv

컬럼 구조 (25컬럼, 공통):
  Gene TSV   [0]=Gene Symbol  [1]=Gene ID
  Region TSV [0]=#ISCA ID     [1]=ISCA Region Name

  [2]  cytoBand
  [3]  Genomic Location GRCh38   "chr22:19756703-19783593"
  [4]  Haploinsufficiency Score  숫자 (0,1,2,3,30,40)
  [5]  HI Score Description      "Sufficient evidence for dosage pathogenicity"
  [6]  HI Disease ID             "OMIM:188400;ORPHA:567"
  [7]  HI Disease Name
  [8]  HI Description
  [9]  HI PMID1
  [10] HI PMID2
  [11] HI PMID3
  [12] HI Notes
  [13] Triplosensitivity Score
  [14] TS Score Description
  [15] TS Disease ID
  [16] TS Disease Name
  [17] TS Description
  [18] TS PMID1
  [19] TS PMID2
  [20] TS PMID3
  [21] TS Notes
  [22] Date Last Evaluated
  [23] Loss phenotype OMIM ID
  [24] TS phenotype OMIM ID
"""

import re
import csv
from io import StringIO
from pathlib import Path

# HI/TS Score 숫자 → 레이블
HI_LABEL = {
    "0":  "No evidence available",
    "1":  "Little evidence for dosage pathogenicity",
    "2":  "Emerging evidence for dosage pathogenicity",
    "3":  "Sufficient evidence for dosage pathogenicity",
    "30": "Gene associated with autosomal recessive phenotype",
    "40": "Dosage sensitivity unlikely",
}


def _col(row: list[str], idx: int) -> str:
    """인덱스 안전 접근 + 공백 제거"""
    return row[idx].strip() if idx < len(row) else ""


def _parse_loc(loc_str: str) -> tuple[str, int, int]:
    """'chr22:19756703-19783593' → ('chr22', 19756703, 19783593)"""
    m = re.match(r"(chr[\dXY]+):(\d+)-(\d+)", loc_str.strip())
    if m:
        return m.group(1), int(m.group(2)), int(m.group(3))
    return "", 0, 0


def _is_header_row(row: list[str]) -> bool:
    """첫 번째 컬럼이 헤더인지 판단"""
    first = _col(row, 0).lstrip("#").strip().lower()
    return first in ("isca id", "gene symbol", "#isca id")


def _read_rows(path: str) -> list[list[str]]:
    """TSV 읽기: # 주석 제거, 헤더 행 제거, 빈 행 제거"""
    rows = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            stripped = line.rstrip("\n")
            # 주석 라인 스킵 (#으로 시작하되 데이터 행이 아닌 것)
            if stripped.startswith("#") and not re.match(r"#ISCA-", stripped):
                continue
            cols = stripped.split("\t")
            if not any(c.strip() for c in cols):
                continue
            # 헤더 행 스킵
            if _is_header_row(cols):
                continue
            rows.append(cols)
    return rows


# ════════════════════════════════════════════════════════════════
# Gene TSV 파서
# ════════════════════════════════════════════════════════════════

def parse_clingen_gene(path: str) -> dict[str, dict]:
    """
    ClinGen_gene_curation_list_GRCh38.tsv → {gene_symbol: {...}}

    반환 필드:
      gene_symbol, ncbi_gene_id, cytoband, chrom, start, end
      hi_score, hi_score_label, hi_disease_id, hi_disease_name
      hi_description, hi_pmids
      ts_score, ts_score_label, ts_disease_id, ts_disease_name
      ts_description, ts_pmids
      date_evaluated, loss_omim_id, ts_omim_id
    """
    if not path or not Path(path).exists():
        print(f"  [ClinGen Gene] 파일 없음: {path}")
        return {}

    lookup: dict[str, dict] = {}
    rows = _read_rows(path)

    for row in rows:
        sym = _col(row, 0)
        if not sym or sym.startswith("#"):
            continue

        chrom, start, end = _parse_loc(_col(row, 3))
        hi_s  = _col(row, 4)
        ts_s  = _col(row, 13)

        # PMIDs 수집 (col 9,10,11)
        hi_pmids = [p for p in [_col(row, 9), _col(row, 10), _col(row, 11)] if p]
        ts_pmids = [p for p in [_col(row, 18), _col(row, 19), _col(row, 20)] if p]

        lookup[sym] = {
            "gene_symbol":      sym,
            "ncbi_gene_id":     _col(row, 1),
            "cytoband":         _col(row, 2),
            "chrom":            chrom,
            "start":            start,
            "end":              end,
            # Haploinsufficiency
            "hi_score":         hi_s,
            "hi_score_label":   HI_LABEL.get(hi_s, _col(row, 5)),
            "hi_score_desc":    _col(row, 5),   # 원문 그대로
            "hi_disease_id":    _col(row, 6),   # "OMIM:188400;ORPHA:567"
            "hi_disease_name":  _col(row, 7),
            "hi_description":   _col(row, 8),
            "hi_pmids":         hi_pmids,
            "hi_notes":         _col(row, 12),
            # Triplosensitivity
            "ts_score":         ts_s,
            "ts_score_label":   HI_LABEL.get(ts_s, _col(row, 14)),
            "ts_score_desc":    _col(row, 14),
            "ts_disease_id":    _col(row, 15),
            "ts_disease_name":  _col(row, 16),
            "ts_description":   _col(row, 17),
            "ts_pmids":         ts_pmids,
            "ts_notes":         _col(row, 21),
            # 메타
            "date_evaluated":   _col(row, 22),
            "loss_omim_id":     _col(row, 23),
            "ts_omim_id":       _col(row, 24),
        }

    print(f"  [ClinGen Gene] {len(lookup)}개 유전자 로드")
    return lookup


# ════════════════════════════════════════════════════════════════
# Region TSV 파서
# ════════════════════════════════════════════════════════════════

def parse_clingen_region(path: str) -> list[dict]:
    """
    ClinGen_region_curation_list_GRCh38.tsv → 목록

    반환 필드 (gene과 동일 구조, isca_id / region_name 추가):
      isca_id, region_name, cytoband, chrom, start, end
      hi_score, hi_score_label, hi_disease_id, hi_disease_name
      hi_description, hi_pmids
      ts_score, ts_score_label, ts_disease_id, ts_disease_name
      ts_description, ts_pmids
      date_evaluated
    """
    from pathlib import Path
    
    if not path or not Path(path).exists():
        print(f"  [ClinGen Region] 파일 없음: {path}")
        return []

    regions = []
    rows = _read_rows(path)

    for row in rows:
        isca_id = _col(row, 0).lstrip("#").strip()
        if not isca_id:
            continue

        chrom, start, end = _parse_loc(_col(row, 3))
        if not chrom:
            continue

        # [MODIFIED] 변경 이유: 주석으로 제공된 실제 TSV 컬럼 순서(0~22)에 맞게 파싱 인덱스를 전면 재배치
        hi_s = _col(row, 4)   # Haploinsufficiency Score
        ts_s = _col(row, 12)  # Triplosensitivity Score

        # PMID 컬럼들 (HI: 6~11, TS: 14~19)
        hi_pmids = [
            p for p in [
                _col(row, 6), _col(row, 7), _col(row, 8), 
                _col(row, 9), _col(row, 10), _col(row, 11)
            ] if p
        ]
        
        ts_pmids = [
            p for p in [
                _col(row, 14), _col(row, 15), _col(row, 16), 
                _col(row, 17), _col(row, 18), _col(row, 19)
            ] if p
        ]
        
        # ISCA ID(0), ISCA Region Name(1), cytoBand(2), Genomic Location(3)
        # Haploinsufficiency Score(4), Haploinsufficiency Description(5)
        # Triplosensitivity Score(12), Triplosensitivity Description(13)
        # Date Last Evaluated(20), Haploinsufficiency Disease ID(21), Triplosensitivity Disease ID(22)

        regions.append({
            "isca_id":          isca_id,
            "region_name":      _col(row, 1),
            "cytoband":         _col(row, 2),
            "chrom":            chrom,
            "start":            start,
            "end":              end,
            
            # Haploinsufficiency
            "hi_score":         hi_s,
            "hi_score_label":   HI_LABEL.get(hi_s, hi_s), # Label이 없으면 원본 Score 표시
            "hi_disease_id":    _col(row, 21),
            "hi_disease_name":  "", # 헤더에 Disease Name이 별도로 없으므로 빈 문자열 처리
            "hi_description":   _col(row, 5),
            "hi_pmids":         hi_pmids,
            "hi_notes":         "",
            
            # Triplosensitivity
            "ts_score":         ts_s,
            "ts_score_label":   HI_LABEL.get(ts_s, ts_s),
            "ts_disease_id":    _col(row, 22),
            "ts_disease_name":  "",
            "ts_description":   _col(row, 13),
            "ts_pmids":         ts_pmids,
            "ts_notes":         "",
            
            "date_evaluated":   _col(row, 20),
        })

    print(f"  [ClinGen Region] {len(regions)}개 영역 로드")
    return regions


def get_region_overlaps(chrom: str, start: int, end: int,
                         regions: list[dict]) -> list[dict]:
    """좌표 overlap으로 ClinGen region 매칭"""
    return [
        r for r in regions
        if r["chrom"] == chrom and r["start"] <= end and r["end"] >= start
    ]