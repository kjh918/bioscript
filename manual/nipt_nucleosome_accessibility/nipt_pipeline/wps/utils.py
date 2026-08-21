"""
wps/utils.py
============
WPS 모듈 공통 유틸리티
- hg38 chrom sizes
- manifest.json 읽기/쓰기
- NPY track 로더 (manifest 경로 우선, fallback prefix 기반)
"""
from __future__ import annotations

import json
import logging
import os
from typing import Optional

import numpy as np

log = logging.getLogger(__name__)

# ── hg38 standard chrom sizes ─────────────────────────────────────────
HG38_CHROM_SIZES: dict[str, int] = {
    "chr1":  248956422, "chr2":  242193529, "chr3":  198295559,
    "chr4":  190214555, "chr5":  181538259, "chr6":  170805979,
    "chr7":  159345973, "chr8":  145138636, "chr9":  138394717,
    "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
    "chr16":  90338345, "chr17":  83257441, "chr18":  80373285,
    "chr19":  58617616, "chr20":  64444167, "chr21":  46709983,
    "chr22":  50818468, "chrX":  156040895, "chrY":  57227415,
}

STANDARD_CHROMS: list[str] = list(HG38_CHROM_SIZES.keys())

# WPS 파라미터 (mode별)
WPS_PARAMS: dict[str, dict] = {
    "L": {"frag_min": 120, "frag_max": 180, "window": 120},
    "S": {"frag_min":  35, "frag_max":  80, "window":  16},
}

MAX_FRAG = 1000

METRIC_SUFFIX: dict[str, str] = {
    "raw":      "raw.npy",
    "cov":      "cov.npy",
    "frag_cov": "frag_cov.npy",
    "norm":     "norm.npy",
}


# ── manifest ──────────────────────────────────────────────────────────
def write_manifest(
    track_dir: str,
    mode: str,
    chrom_stats: dict[str, dict],
    extra: dict | None = None,
) -> str:
    prm = WPS_PARAMS[mode]
    data: dict = {
        "mode":        mode,
        "frag_min":    prm["frag_min"],
        "frag_max":    prm["frag_max"],
        "window":      prm["window"],
        "chromosomes": chrom_stats,
    }
    if extra:
        data.update(extra)
    path = os.path.join(track_dir, "manifest.json")
    with open(path, "w") as fh:
        json.dump(data, fh, indent=2)
    return path


def read_manifest(wps_dir: str, mode: str) -> dict:
    """
    {wps_dir}/{mode}/manifest.json 로드.
    없으면 빈 dict 반환.
    """
    path = os.path.join(wps_dir, mode, "manifest.json")
    if not os.path.isfile(path):
        log.warning("manifest 없음: %s", path)
        return {}
    with open(path) as fh:
        return json.load(fh)


# ── NPY track 경로 조립 ────────────────────────────────────────────────
def track_path(
    wps_dir: str,
    mode: str,
    chrom: str,
    metric: str,
    prefix: str = "",
) -> str:
    """
    {wps_dir}/{mode}/{prefix}.{chrom}.{metric}.npy 경로 반환.
    wps_dir 은 항상 루트 (L/ S/ 상위)를 가리켜야 합니다.
    """
    fname = f"{prefix}.{chrom}" if prefix else chrom
    return os.path.join(wps_dir, mode, f"{fname}.{METRIC_SUFFIX[metric]}")


# ── NPY track 로더 ────────────────────────────────────────────────────
def load_track(
    wps_dir: str,
    mode: str,
    chrom: str,
    metric: str,
    prefix: str = "",
    manifest: dict | None = None,
) -> Optional[np.ndarray]:
    """
    NPY 1-D 배열 로드.

    우선순위
    --------
    1. manifest["chromosomes"][chrom][metric] 절대경로 (파일 존재 확인)
    2. track_path(wps_dir, mode, chrom, metric, prefix) fallback

    Parameters
    ----------
    wps_dir  : compute.run() out_dir 루트 (L/ S/ 의 상위)
    mode     : "L" | "S"
    """
    def _load(p: str) -> Optional[np.ndarray]:
        if not os.path.isfile(p):
            return None
        try:
            arr = np.load(p, mmap_mode="r", allow_pickle=False)
            return arr if arr.ndim == 1 else None
        except Exception as e:
            log.warning("NPY 로드 실패 %s: %s", p, e)
            return None

    # 1. manifest 절대경로 우선
    if manifest:
        entry = manifest.get("chromosomes", {}).get(chrom, {})
        if isinstance(entry, dict):
            p = entry.get(metric, "")
            if p:
                arr = _load(p)
                if arr is not None:
                    return arr
                log.debug("manifest 경로 miss (파일 없음): %s", p)

    # 2. 경로 직접 조립 fallback
    p   = track_path(wps_dir, mode, chrom, metric, prefix)
    arr = _load(p)
    if arr is None:
        log.debug("fallback 경로 miss: %s", p)
    return arr