"""
utils.py
─────────
공통 유틸리티: HTTP 세션, rate-limit GET, 공유 캐시
"""

import time
import functools
from typing import Any

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# ── API 엔드포인트 ─────────────────────────────────────────────
EUTILS   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
PUBMED   = "https://pubmed.ncbi.nlm.nih.gov"
PMC_BASE = "https://www.ncbi.nlm.nih.gov/pmc/articles"
UNIPROT  = "https://rest.uniprot.org/uniprotkb"

# ── GRCh38 염색체 전체 좌표 (하드코딩) ───────────────────────
# 출처: NCBI GRCh38.p14 (GCA_000001405.29)
GRCH38_CHROM_COORDS: dict[str, dict] = {
    "chr1":  {"start": 1, "end": 248956422, "length_mb": 248.96},
    "chr2":  {"start": 1, "end": 242193529, "length_mb": 242.19},
    "chr3":  {"start": 1, "end": 198295559, "length_mb": 198.30},
    "chr4":  {"start": 1, "end": 190214555, "length_mb": 190.21},
    "chr5":  {"start": 1, "end": 181538259, "length_mb": 181.54},
    "chr6":  {"start": 1, "end": 170805979, "length_mb": 170.81},
    "chr7":  {"start": 1, "end": 159345973, "length_mb": 159.35},
    "chr8":  {"start": 1, "end": 145138636, "length_mb": 145.14},
    "chr9":  {"start": 1, "end": 138394717, "length_mb": 138.39},
    "chr10": {"start": 1, "end": 133797422, "length_mb": 133.80},
    "chr11": {"start": 1, "end": 135086622, "length_mb": 135.09},
    "chr12": {"start": 1, "end": 133275309, "length_mb": 133.28},
    "chr13": {"start": 1, "end": 114364328, "length_mb": 114.36},
    "chr14": {"start": 1, "end": 107043718, "length_mb": 107.04},
    "chr15": {"start": 1, "end": 101991189, "length_mb": 101.99},
    "chr16": {"start": 1, "end": 90338345,  "length_mb": 90.34},
    "chr17": {"start": 1, "end": 83257441,  "length_mb": 83.26},
    "chr18": {"start": 1, "end": 80373285,  "length_mb": 80.37},
    "chr19": {"start": 1, "end": 58617616,  "length_mb": 58.62},
    "chr20": {"start": 1, "end": 64444167,  "length_mb": 64.44},
    "chr21": {"start": 1, "end": 46709983,  "length_mb": 46.71},
    "chr22": {"start": 1, "end": 50818468,  "length_mb": 50.82},
    "chrX":  {"start": 1, "end": 156040895, "length_mb": 156.04},
    "chrY":  {"start": 1, "end": 57227415,  "length_mb": 57.23},
}

# p arm / q arm centromere 대략 위치 (UCSC hg38 기준, Mb)
# centromere 앞이 p arm, 뒤가 q arm
GRCH38_CENTROMERE: dict[str, int] = {
    "chr1": 123400000,  "chr2": 92500000,  "chr3": 90900000,
    "chr4": 49700000,   "chr5": 46800000,  "chr6": 59800000,
    "chr7": 60600000,   "chr8": 45200000,  "chr9": 43000000,
    "chr10": 39800000,  "chr11": 53400000, "chr12": 35500000,
    "chr13": 17700000,  "chr14": 17200000, "chr15": 19000000,
    "chr16": 36800000,  "chr17": 25100000, "chr18": 18500000,
    "chr19": 26200000,  "chr20": 28100000, "chr21": 12000000,
    "chr22": 15000000,  "chrX": 61000000,  "chrY": 10400000,
}


def get_chrom_marker(chrom: str) -> dict:
    """
    염색체 이름 → GRCh38 전체 좌표 마커 반환.
    TargetChromosome feature_type에 사용.
    """
    key = chrom if chrom.startswith("chr") else f"chr{chrom}"
    coords = GRCH38_CHROM_COORDS.get(key, {})
    cen    = GRCH38_CENTROMERE.get(key, 0)
    if not coords:
        return {}
    return {
        "chromosome":   key,
        "start":        coords["start"],
        "end":          coords["end"],
        "length_mb":    coords["length_mb"],
        "centromere":   cen,
        "p_arm_end":    cen,
        "q_arm_start":  cen + 1,
        "assembly":     "GRCh38",
    }


# ── HTTP 세션 ─────────────────────────────────────────────────

_session: requests.Session | None = None
_email: str = "your@email.com"


def init_session(email: str) -> None:
    global _session, _email
    _email = email
    s = requests.Session()
    retry = Retry(
        total=5, backoff_factor=1.5,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
    )
    s.mount("https://", HTTPAdapter(max_retries=retry))
    s.mount("http://",  HTTPAdapter(max_retries=retry))
    s.headers.update({"User-Agent": f"nipt-pipeline/1.0 (mailto:{email})"})
    _session = s


def get_session() -> requests.Session:
    if _session is None:
        init_session(_email)
    return _session


def http_get(url: str, params: dict | None = None,
             delay: float = 0.35) -> requests.Response:
    time.sleep(delay)
    r = get_session().get(url, params=params, timeout=20)
    r.raise_for_status()
    return r


# ── 간단 메모이제이션 캐시 ────────────────────────────────────

_cache: dict[str, Any] = {}


def cached(key: str, fn):
    """key가 캐시에 없으면 fn() 호출 후 저장."""
    if key not in _cache:
        _cache[key] = fn()
    return _cache[key]


def clear_cache():
    _cache.clear()
