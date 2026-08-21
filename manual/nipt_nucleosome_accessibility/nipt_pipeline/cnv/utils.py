"""
cnv/utils.py
공통 유틸리티: hg38 chrom sizes, bin 생성, TSV I/O
"""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterator

# hg38 standard chrom sizes (UCSC)
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


def iter_bins(
    bin_size: int = 100_000,
    chroms: list[str] | None = None,
) -> Iterator[tuple[str, int, int]]:
    """
    (chrom, start, end) 튜플을 순서대로 yield.
    end는 chrom 경계에서 clamp.
    """
    target = chroms or STANDARD_CHROMS
    for chrom in target:
        chrom_len = HG38_CHROM_SIZES[chrom]
        start = 0
        while start < chrom_len:
            end = min(start + bin_size, chrom_len)
            yield chrom, start, end
            start = end


def read_tsv(path: str | Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(rows: list[dict], path: str | Path, fieldnames: list[str]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)