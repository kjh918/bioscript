"""
TSV / CSV input parsers.

Expected columns (case-insensitive, flexible aliases):
  chrom / chromosome / chr      – chromosome name
  pos / position / start        – genomic position (bp)
  cn / copy_number / ratio      – copy number value
  baf / vaf / b_allele_freq     – B-allele frequency (optional)

Sample metadata TSV columns (for --sample flag):
  sample_id, sex, chr, type, cn, iscn, start, stop, color (CNV rows)
  sample_id, chr, gene, start, stop, color                 (gene rows; type == "gene")
"""

from __future__ import annotations

import io
import base64
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .models import SampleData, CnvEvent, GeneAnnotation
from .reference import CHROM_SIZES


# ---------------------------------------------------------------------------
# Column alias normalisation
# ---------------------------------------------------------------------------
_CHROM_ALIASES = {"chromosome", "chr", "chrom"}
_POS_ALIASES   = {"pos", "position", "start", "chromstart", "chrom_start"}
_CN_ALIASES    = {"cn", "copy_number", "copynumber", "ratio", "log2"}
_BAF_ALIASES   = {"baf", "vaf", "b_allele_freq", "allele_freq"}


def _normalise_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Lower-case and alias-map column names in-place."""
    rename: dict[str, str] = {}
    cols = {c.strip().lower(): c for c in df.columns}

    for alias in _CHROM_ALIASES:
        if alias in cols and "chrom" not in rename.values():
            rename[cols[alias]] = "chrom"
            break
    for alias in _POS_ALIASES:
        if alias in cols and "pos" not in rename.values():
            rename[cols[alias]] = "pos"
            break
    for alias in _CN_ALIASES:
        if alias in cols and "cn" not in rename.values():
            rename[cols[alias]] = "cn"
            break
    for alias in _BAF_ALIASES:
        if alias in cols and "baf" not in rename.values():
            rename[cols[alias]] = "baf"
            break

    return df.rename(columns=rename)


def validate_dataframe(df: pd.DataFrame) -> tuple[bool, list[str]]:
    """Return (ok, list_of_error_messages)."""
    errors: list[str] = []
    for col in ("pos", "cn"):
        if col not in df.columns:
            errors.append(f"필수 컬럼 없음: '{col}'")
    return (len(errors) == 0, errors)


def load_tsv(path_or_str: "str | Path", *, sep: Optional[str] = None) -> pd.DataFrame:
    """
    Load a CN/BAF TSV or CSV from a file path.
    sep=None → auto-detect from extension (.tsv → tab, else comma).
    """
    path = Path(path_or_str)
    if sep is None:
        sep = "\t" if path.suffix.lower() in (".tsv", ".txt") else ","
    df = pd.read_csv(path, sep=sep)
    df = _normalise_columns(df)
    ok, errs = validate_dataframe(df)
    if not ok:
        raise ValueError("; ".join(errs))
    keep = [c for c in ("chrom", "pos", "cn", "baf") if c in df.columns]
    df = df[keep].dropna(subset=["pos", "cn"])
    df["pos"] = df["pos"].astype(int)
    return df


def load_base64_upload(contents: str, filename: str) -> pd.DataFrame:
    """Parse a Dash dcc.Upload base64 payload."""
    _, payload = contents.split(",", 1)
    decoded = base64.b64decode(payload).decode("utf-8")
    sep = "\t" if str(filename).lower().endswith(".tsv") else ","
    df = pd.read_csv(io.StringIO(decoded), sep=sep)
    df = _normalise_columns(df)
    ok, errs = validate_dataframe(df)
    if not ok:
        raise ValueError("; ".join(errs))
    keep = [c for c in ("chrom", "pos", "cn", "baf") if c in df.columns]
    df = df[keep].dropna(subset=["pos", "cn"])
    df["pos"] = df["pos"].astype(int)
    return df


# ---------------------------------------------------------------------------
# Sample metadata TSV loader
# ---------------------------------------------------------------------------
def load_sample_tsv(path: "str | Path") -> SampleData:
    """
    Load sample metadata (events + genes) from a structured TSV.

    Required header row:
        sample_id  sex  chr  type  cn  iscn  start  stop  color  gene

    Rows with type in (trisomy|monosomy|partial_gain|partial_loss) → CnvEvent
    Rows with type == "gene"                                        → GeneAnnotation
    """
    path = Path(path)
    df = pd.read_csv(path, sep="\t")
    df.columns = [c.strip().lower() for c in df.columns]

    sample_id = str(df["sample_id"].iloc[0])
    sex_raw = str(df["sex"].iloc[0]).strip().lower()
    sex = "female" if sex_raw.startswith("f") else "male"

    events: list[CnvEvent] = []
    genes:  list[GeneAnnotation] = []

    for _, row in df.iterrows():
        t = str(row.get("type", "")).strip().lower()
        if t == "gene":
            genes.append(GeneAnnotation(
                chr=str(row["chr"]),
                name=str(row.get("gene", row.get("iscn", "?"))),
                start=int(row["start"]),
                stop=int(row["stop"]),
                color=str(row.get("color", "#68D391")),
            ))
        elif t in ("trisomy", "monosomy", "partial_gain", "partial_loss"):
            events.append(CnvEvent(
                chr=str(row["chr"]),
                type=t,
                cn=int(row.get("cn", 2)),
                iscn=str(row.get("iscn", "")),
                start=int(row["start"]),
                stop=int(row["stop"]),
                color=str(row.get("color", "#FC8181")),
            ))

    return SampleData(id=sample_id, sex=sex, events=events, genes=genes)
