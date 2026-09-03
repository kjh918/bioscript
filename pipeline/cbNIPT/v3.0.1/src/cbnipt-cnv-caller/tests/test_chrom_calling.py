import numpy as np
import pandas as pd

from cbnipt_cnv_caller.calling.chrom_calling import compute_chrom_summary


def test_compute_chrom_summary_aggregates_present_columns_and_sorts():
    df = pd.DataFrame({
        "chrom": ["chr2"] * 3 + ["chr1"] * 2 + ["chrX"] * 4,
        "log2_chrom_norm": [0.1, 0.2, 0.3, 0.0, 0.0, -1.0, -1.0, -1.0, -1.0],
        "bin_BAF": [0.5] * 9,
    })
    summary = compute_chrom_summary(df)

    # sorted by chromosome order
    assert list(summary["chrom"]) == ["chr1", "chr2", "chrX"]
    # per-chrom bin counts
    assert list(summary["n_bins"]) == [2, 3, 4]
    # median/mean columns generated only for present candidate cols
    assert "log2_chrom_norm_median" in summary.columns
    assert "bin_BAF_mean" in summary.columns
    assert np.isclose(summary.loc[summary.chrom == "chr2", "log2_chrom_norm_median"].iloc[0], 0.2)
