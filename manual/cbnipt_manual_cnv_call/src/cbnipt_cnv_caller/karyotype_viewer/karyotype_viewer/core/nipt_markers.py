"""
NIPT marker DB loader.
marker TSV(GCX_NIPT_NO, NIPT_ID, SYNDROME, NIPT_GROUP, CHROMOSOME,
           FEATURE_NAME, FEATURE_TYPE, GENOMIC_POS_START, GENOMIC_POS_END)
를 읽어 syndrome별 구조체로 제공.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
import pandas as pd


CALL_COLORS = {
    "LOW_RISK":     "#38A169",
    "SUSPECTED": "#DD6B20",
    "HIGH_RISK":   "#E53E3E",
    "UNKNOWN":    "#718096",
}

GROUP_COLORS = {
    "Autosome Abnormality":       "#FC8181",
    "Sex Chromosome Abnormality": "#90CDF4",
    "Micro Deletion":             "#F6AD55",
}


@dataclass
class MarkerFeature:
    feature_name: str
    feature_type: str   # TargetChromosome / PrimaryTargetRegion / CoreRegion / CoreGene
    chromosome:   str
    start:        int
    end:          int

    @property
    def size_mb(self) -> float:
        return (self.end - self.start) / 1e6


@dataclass
class NiptSyndrome:
    nipt_id:   str
    syndrome:  str
    group:     str          # Autosome Abnormality / Sex Chromosome Abnormality / Micro Deletion
    features:  list[MarkerFeature] = field(default_factory=list)
    call:      str = "UNKNOWN"   # LOW_RISK / SUSPECTED / HIGH_RISK / UNKNOWN
    cn_value:  Optional[float] = None

    @property
    def primary_chrom(self) -> str:
        """첫 번째 TargetChromosome 또는 PrimaryTargetRegion 의 염색체."""
        for f in self.features:
            if f.feature_type in ("TargetChromosome", "PrimaryTargetRegion"):
                return f.chromosome.replace("chr", "")
        return self.features[0].chromosome.replace("chr", "") if self.features else "?"

    @property
    def call_color(self) -> str:
        return CALL_COLORS.get(self.call, CALL_COLORS["UNKNOWN"])

    @property
    def group_color(self) -> str:
        return GROUP_COLORS.get(self.group, "#CBD5E0")


def load_marker_tsv(path: "str | Path") -> dict[str, NiptSyndrome]:
    """
    marker TSV → {nipt_id: NiptSyndrome} dict.
    """
    df = pd.read_csv(Path(path), sep="\t")
    df.columns = [c.strip() for c in df.columns]

    syndromes: dict[str, NiptSyndrome] = {}
    for _, row in df.iterrows():
        nid = str(row["NIPT_ID"]).strip()
        if nid not in syndromes:
            syndromes[nid] = NiptSyndrome(
                nipt_id  = nid,
                syndrome = str(row["SYNDROME"]).strip(),
                group    = str(row["NIPT_GROUP"]).strip(),
            )
        syndromes[nid].features.append(MarkerFeature(
            feature_name = str(row["FEATURE_NAME"]).strip(),
            feature_type = str(row["FEATURE_TYPE"]).strip(),
            chromosome   = str(row["CHROMOSOME"]).strip(),
            start        = int(row["GENOMIC_POS_START"]),
            end          = int(row["GENOMIC_POS_END"]),
        ))
    return syndromes
