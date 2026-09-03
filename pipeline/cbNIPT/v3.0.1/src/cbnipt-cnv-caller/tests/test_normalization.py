import numpy as np
import pandas as pd

from cbnipt_cnv_caller.calling.normalization import apply_qc_bin, normalize_by_chrom_with_sex


def test_apply_qc_bin_drops_filtered_and_low_coverage_but_protects_sex():
    df = pd.DataFrame([
        {"chrom": "chr1", "is_filtered": False, "qc_pass_fragments": 10, "breadth_ratio": 0.9},
        {"chrom": "chr1", "is_filtered": True,  "qc_pass_fragments": 10, "breadth_ratio": 0.9},
        {"chrom": "chr1", "is_filtered": False, "qc_pass_fragments": 0,  "breadth_ratio": 0.0},
        {"chrom": "chrX", "is_filtered": False, "qc_pass_fragments": 0,  "breadth_ratio": 0.0},
    ])
    out = apply_qc_bin(df, min_depth=1, min_coverage=0.1)
    assert list(out["chrom"]) == ["chr1", "chrX"]


def test_normalize_by_chrom_centres_balanced_diploid_on_two():
    # 5 autosomes, 20 flat bins each, log2_corrected == 0  => perfectly diploid
    rows = []
    for c in ["chr1", "chr2", "chr3", "chr4", "chr5"]:
        for i in range(20):
            rows.append({"chrom": c, "start": i * 100_000, "end": (i + 1) * 100_000,
                         "log2_corrected": 0.0})
    df = pd.DataFrame(rows)

    df_norm, sex_tag, x_ratio, y_ratio = normalize_by_chrom_with_sex(df, value_col="log2_corrected")

    assert np.allclose(df_norm["copy_number_signal"], 2.0, atol=1e-6)
    assert np.allclose(df_norm["log2_chrom_norm"], 0.0, atol=1e-6)
