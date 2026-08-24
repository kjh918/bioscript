"""
Demo CN/BAF DataFrame generator.
Used when no real data file is supplied (preview / testing).
"""

from __future__ import annotations
from typing import Optional

import numpy as np
import pandas as pd

from .models import SampleData
from .reference import CHROM_SIZES


def demo_dataframe(
    sample: SampleData,
    chrom: str,
    start_bp: Optional[int] = None,
    end_bp: Optional[int] = None,
) -> pd.DataFrame:
    """
    Generate synthetic CN + BAF points for *chrom* using event definitions
    from *sample*.

    Parameters
    ----------
    sample   : SampleData instance (provides events for CN/BAF shaping)
    chrom    : chromosome name without "chr" prefix
    start_bp : start position (bp); None → chromosome start
    end_bp   : end position (bp);   None → chromosome end

    Returns
    -------
    pd.DataFrame with columns: chrom, pos, cn, baf
    """
    chrom = str(chrom).replace("chr", "")
    size  = CHROM_SIZES[chrom]

    start_bp = max(1, int(start_bp)) if start_bp is not None else 1
    end_bp   = min(size, int(end_bp)) if end_bp is not None else size
    if end_bp <= start_bp:
        end_bp = min(size, start_bp + 1_000_000)

    region_mode = (start_bp > 1) or (end_bp < size)
    n   = 420 if region_mode else 1_600
    rng = np.random.default_rng(sum(ord(x) for x in chrom) + 111)
    pos = np.linspace(start_bp, end_bp, n).astype(int)

    cn  = np.full(n, 2.0)
    baf = np.full(n, 0.5)

    for ev in sample.events_for_chrom(chrom):
        mask = (pos >= ev.start) & (pos <= ev.stop)
        if ev.type in ("trisomy", "partial_gain"):
            cn[mask]  = 3.0
            baf[mask] = rng.choice([0.33, 0.67], mask.sum())
        elif ev.type in ("monosomy", "partial_loss"):
            cn[mask]  = 1.0
            baf[mask] = rng.choice([0.03, 0.97], mask.sum())

    return pd.DataFrame({
        "chrom": chrom,
        "pos":   pos,
        "cn":    np.clip(cn  + rng.normal(0, 0.18, n), 0, 6),
        "baf":   np.clip(baf + rng.normal(0, 0.035, n), 0, 1),
    })
