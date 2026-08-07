"""
Typed data models for karyotype report input.
No UI / Dash dependency.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Literal

from .reference import CHROM_SIZES, FEMALE_CHROMS, MALE_CHROMS


Sex = Literal["female", "male"]
CnvType = Literal["trisomy", "monosomy", "partial_gain", "partial_loss"]


@dataclass
class CnvEvent:
    chr: str
    type: CnvType
    cn: int
    iscn: str
    start: int
    stop: int
    color: str = "#FC8181"

    def __post_init__(self) -> None:
        self.chr = str(self.chr).replace("chr", "")
        if self.chr not in CHROM_SIZES:
            raise ValueError(f"Unknown chromosome: {self.chr!r}")
        if self.start < 1:
            self.start = 1
        if self.stop > CHROM_SIZES[self.chr]:
            self.stop = CHROM_SIZES[self.chr]

    @property
    def length_bp(self) -> int:
        return self.stop - self.start


@dataclass
class GeneAnnotation:
    chr: str
    name: str
    start: int
    stop: int
    color: str = "#68D391"

    def __post_init__(self) -> None:
        self.chr = str(self.chr).replace("chr", "")
        if self.chr not in CHROM_SIZES:
            raise ValueError(f"Unknown chromosome: {self.chr!r}")


@dataclass
class SampleData:
    id: str
    sex: Sex
    events: list[CnvEvent] = field(default_factory=list)
    genes: list[GeneAnnotation] = field(default_factory=list)

    @property
    def display_chroms(self) -> list[str]:
        return FEMALE_CHROMS if self.sex == "female" else MALE_CHROMS

    def events_for_chrom(self, chrom: str) -> list[CnvEvent]:
        chrom = str(chrom).replace("chr", "")
        return [e for e in self.events if e.chr == chrom]

    def genes_for_chrom(self, chrom: str) -> list[GeneAnnotation]:
        chrom = str(chrom).replace("chr", "")
        return [g for g in self.genes if g.chr == chrom]
