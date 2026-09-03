from cbnipt_cnv_caller.calling.sex_calling import (
    determine_fetal_sex,
    classify_sex_chromosome,
    _normalize_sex_tag,
)


def test_determine_fetal_sex():
    assert determine_fetal_sex(x_cn=2.0, y_cn=0.0) == "XX"
    assert determine_fetal_sex(x_cn=1.0, y_cn=1.0) == "XY"
    # ambiguous -> Unknown (both above/below thresholds)
    assert determine_fetal_sex(x_cn=1.0, y_cn=0.0) == "Unknown"


def test_normalize_sex_tag():
    assert _normalize_sex_tag("Male (XY)") == "XY"
    assert _normalize_sex_tag("Female (XX)") == "XX"
    assert _normalize_sex_tag("weird") == "Unknown"


def test_classify_sex_chromosome_hardcalls():
    # normal female X (x~2, y~0) matches no aneuploidy rule
    assert classify_sex_chromosome("chrX", x_cn=2.0, y_cn=0.0) is None

    # Turner (X0): x in [0.5, 1.4], y <= 0.4
    turner = classify_sex_chromosome("chrX", x_cn=1.0, y_cn=0.1)
    assert turner is not None and turner["pattern_name"].startswith("Turner")
    assert turner["call"] == "ABNORMAL"

    # Klinefelter (XXY): x in [1.6, 2.4], y in [0.4, 1.5]
    kline = classify_sex_chromosome("chrX", x_cn=2.0, y_cn=1.0)
    assert kline is not None and kline["pattern_name"].startswith("Klinefelter")

    # non sex chromosome -> always None
    assert classify_sex_chromosome("chr21", x_cn=1.0, y_cn=1.0) is None
