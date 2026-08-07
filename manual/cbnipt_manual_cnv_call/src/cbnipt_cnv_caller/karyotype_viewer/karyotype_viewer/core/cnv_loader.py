"""
CNV data loader: --cnv 디렉토리에서 chr별 TSV를 패턴 매칭으로 로드.
{dir}/cnv_chr{chrom}.tsv  또는  {dir}/*chr{chrom}*.tsv  패턴 지원.
"""

from __future__ import annotations
from pathlib import Path
from typing import Optional
import re
import pandas as pd
from .parser import load_tsv
from .reference import ALL_CHROMS


def load_cnv_dir(directory: "str | Path") -> dict[str, pd.DataFrame]:
    """
    디렉토리에서 염색체별 CN/BAF TSV를 로드해 {chrom: DataFrame} dict 반환.
    파일명에 chr{N} 패턴이 있으면 자동 매핑.
    """
    d = Path(directory)
    if not d.is_dir():
        raise ValueError(f"Not a directory: {d}")

    result: dict[str, pd.DataFrame] = {}
    chrom_re = re.compile(r"chr([0-9]{1,2}|X|Y)", re.IGNORECASE)


    for f in sorted(d.glob("*.tsv")) :
        m = chrom_re.search(f.stem)
        if not m:
            continue
        chrom = m.group(1).upper()
        if chrom not in ALL_CHROMS:
            continue
        try:
            df = load_tsv(f)
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
    """
    cnv_data dict에서 chrom DataFrame 반환.
    없으면 fallback_fn() 호출 (demo_dataframe 등).
    """
    chrom = str(chrom).replace("chr", "").upper()
    if chrom in cnv_data:
        return cnv_data[chrom]
    if fallback_fn:
        return fallback_fn()
    return pd.DataFrame(columns=["pos", "cn", "baf"])
