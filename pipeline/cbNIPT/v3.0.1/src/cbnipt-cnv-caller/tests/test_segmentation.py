import numpy as np
import pandas as pd

from cbnipt_cnv_caller.calling.segmentation import assign_cn_state, rolling_micro_cnv_segmentation


def test_assign_cn_state_autosome_and_sex_ploidy():
    seg = pd.DataFrame([
        {"chrom": "chr1", "seg_mean": 0.0},    # ratio 1 -> cn 2 -> NEUT
        {"chrom": "chr1", "seg_mean": 1.0},    # ratio 2 -> cn 4 -> AMP
        {"chrom": "chr1", "seg_mean": -2.0},   # ratio .25 -> cn 0 -> DEL
        {"chrom": "chrX", "seg_mean": 0.0},    # male X: expected ploidy 1 -> cn 1 -> NEUT
    ])
    out = assign_cn_state(seg, baseline_ploidy=2, sex_tag="XY")
    assert list(out["cnv_call"]) == ["NEUT", "AMP", "DEL", "NEUT"]
    assert out.loc[3, "expected_ploidy"] == 1


def test_rolling_micro_cnv_detects_constant_deletion_block():
    n = 8
    df = pd.DataFrame({
        "chrom": ["chr7"] * n,
        "start": [i * 100_000 for i in range(n)],
        "end": [(i + 1) * 100_000 for i in range(n)],
        "normalized_count": [1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
    })
    segs = rolling_micro_cnv_segmentation(df, signal_col="normalized_count")
    assert len(segs) == 1
    assert segs.iloc[0]["cnv_call"] == "DEL"
    assert segs.iloc[0]["length"] >= 400_000
