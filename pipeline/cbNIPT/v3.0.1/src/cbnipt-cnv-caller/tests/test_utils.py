import numpy as np

from cbnipt_cnv_caller.common.utils import chrom_key, sort_chroms, safe_log2fc


def test_chrom_key_orders_autosomes_then_sex():
    assert chrom_key("chr1") == 1
    assert chrom_key("2") == 2
    assert chrom_key("chrX") == 23
    assert chrom_key("chrY") == 24
    assert chrom_key("chrUn_junk") == 99


def test_sort_chroms():
    got = sort_chroms(["chrX", "chr2", "chr10", "chr1", "chrY"])
    assert got == ["chr1", "chr2", "chr10", "chrX", "chrY"]


def test_safe_log2fc_no_zero_division():
    assert safe_log2fc(4.0, 1.0) == np.log2((4.0 + 1e-9) / (1.0 + 1e-9))
    # ref == 0 must not raise
    assert np.isfinite(safe_log2fc(1.0, 0.0))
