import pandas as pd

from cbnipt_cnv_caller.calling.calculate_tstv import calculate_titv_ratio


def test_titv_counts_transitions_and_transversions():
    df = pd.DataFrame([
        {"ref": "A", "alt": "G"},   # Ti
        {"ref": "C", "alt": "T"},   # Ti
        {"ref": "T", "alt": "C"},   # Ti
        {"ref": "G", "alt": "A"},   # Ti
        {"ref": "A", "alt": "C"},   # Tv
        {"ref": "G", "alt": "T"},   # Tv
        {"ref": "A", "alt": "AT"},  # indel -> ignored
        {"ref": "A", "alt": "A"},   # not a variant -> ignored
    ])
    out = calculate_titv_ratio(df)
    assert out["ti_count"] == 4
    assert out["tv_count"] == 2
    assert out["titv_ratio"] == 2.0


def test_titv_empty_input():
    assert calculate_titv_ratio(pd.DataFrame()) == {"ti_count": 0, "tv_count": 0, "titv_ratio": 0.0}
