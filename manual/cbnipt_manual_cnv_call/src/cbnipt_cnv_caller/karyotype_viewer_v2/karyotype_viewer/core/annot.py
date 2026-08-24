"""
Builders that convert SampleData → Ideogram rawAnnots / dropdown options / ISCN string.
No Dash / Plotly imports.
"""

from __future__ import annotations
import json

from .models import SampleData
from .reference import ALL_CHROMS


def build_iscn(sample: SampleData) -> str:
    """Build an ISCN-like summary string from sample events."""
    total = 46
    parts: list[str] = []
    for ev in sample.events:
        if ev.type == "trisomy":
            total += 1
        elif ev.type == "monosomy":
            total -= 1
        parts.append(ev.iscn)
    sex = "XX" if sample.sex == "female" else "XY"
    suffix = "," + ",".join(parts) if parts else ""
    return f"{total},{sex}{suffix}"


def build_raw_annots(sample: SampleData) -> dict:
    """
    Return a rawAnnots dict compatible with Dash-Bio Ideogram's
    ``annotations`` prop.

    Structure
    ---------
    {
        "keys": ["name", "start", "length", "trackIndex", "color"],
        "annots": [
            {"chr": "1", "annots": [[name, start, length, trackIndex, color], ...]},
            ...
        ]
    }
    """
    by_chr: dict[str, list] = {c: [] for c in ALL_CHROMS}

    for ev in sample.events:
        by_chr[ev.chr].append([
            ev.iscn,
            ev.start,
            max(1, ev.stop - ev.start),
            0,          # trackIndex 0 = CNV track
            ev.color,
        ])

    for gene in sample.genes:
        by_chr[gene.chr].append([
            gene.name,
            gene.start,
            max(50_000, gene.stop - gene.start),
            1,          # trackIndex 1 = Gene track
            gene.color,
        ])

    return {
        "keys": ["name", "start", "length", "trackIndex", "color"],
        "annots": [
            {"chr": c, "annots": by_chr[c]}
            for c in ALL_CHROMS
            if by_chr[c]
        ],
    }


def annotation_options_for_chrom(sample: SampleData, chrom: str) -> list[dict]:
    """
    Build dcc.Dropdown options list for the annotation selector.
    Each option value is a JSON string ``{"start": int, "end": int}``.
    """
    chrom = str(chrom).replace("chr", "")
    options: list[dict] = []

    for ev in sample.events_for_chrom(chrom):
        options.append({
            "label": (
                f"CNV | {ev.iscn} | "
                f"{ev.start/1e6:.3f}–{ev.stop/1e6:.3f} Mb"
            ),
            "value": json.dumps({"start": ev.start, "end": ev.stop}),
        })

    for gene in sample.genes_for_chrom(chrom):
        options.append({
            "label": (
                f"Gene | {gene.name} | "
                f"{gene.start/1e6:.3f}–{gene.stop/1e6:.3f} Mb"
            ),
            "value": json.dumps({"start": gene.start, "end": gene.stop}),
        })

    return options
