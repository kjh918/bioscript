"""
CNV data loader.

두 가지 입력 형식 지원:
  1. 실제 파이프라인: nipt_loader.load_cnv_tsv() 사용 (전체 genome 단일 TSV)
  2. 구버전: 디렉토리 내 cnv_chr{N}.tsv 파일들
"""

from __future__ import annotations
from pathlib import Path
from typing import Optional
import re
import pandas as pd
from .reference import ALL_CHROMS

# 실제 파이프라인 컬럼 → 표준 컬럼
_COL_ALIASES = {
    "copy_number_signal": "cn",
    "bin_baf":            "baf",
    "bin_BAF":            "baf",
    "observed_cn":        "cn",
    "start":              "pos",
}


def _normalize_df(df: pd.DataFrame) -> pd.DataFrame:
    """컬럼명 정규화."""
    rename = {}
    for col in df.columns:
        key = col.lower()
        if key in _COL_ALIASES:
            rename[col] = _COL_ALIASES[key]
        elif col in _COL_ALIASES:
            rename[col] = _COL_ALIASES[col]
    df = df.rename(columns=rename)
    if "pos" not in df.columns and "start" in df.columns:
        df["pos"] = df["start"]
    if "cn" not in df.columns and "copy_number" in df.columns:
        df["cn"] = df["copy_number"]
    return df


def load_cnv_dir(directory: "str | Path") -> dict[str, pd.DataFrame]:
    """디렉토리에서 염색체별 TSV 로드 → {chrom: DataFrame}."""
    d = Path(directory)
    if not d.is_dir():
        raise ValueError(f"Not a directory: {d}")

    result: dict[str, pd.DataFrame] = {}
    chrom_re = re.compile(r"chr([0-9]{1,2}|X|Y)", re.IGNORECASE)

    for f in sorted(d.glob("*.tsv")):
        m = chrom_re.search(f.stem)
        if not m:
            continue
        chrom = m.group(1).upper()
        if chrom not in ALL_CHROMS:
            continue
        try:
            df = pd.read_csv(f, sep="\t")
            df = _normalize_df(df)
            result[chrom] = df
            print(f"[cnv_loader] {f.name} → chr{chrom}  ({len(df):,} rows)")
        except Exception as e:
            print(f"[cnv_loader] SKIP {f.name}: {e}")

    return result


def get_chrom_df(
    cnv_data: dict[str, pd.DataFrame],
    chrom: str,
    fallback_fn=None,
) -> pd.DataFrame:
    chrom = str(chrom).replace("chr", "").upper()
    if chrom in cnv_data:
        return cnv_data[chrom]
    if fallback_fn:
        return fallback_fn()
    return pd.DataFrame(columns=["pos", "cn", "baf"])
