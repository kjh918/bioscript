"""Smoke test: the karyotype v2 report module builds a self-contained HTML report."""

from karyotype_viewer.demo_sample import make_demo_sample
from karyotype_viewer.reporter import save_report_dir


def test_save_report_dir_writes_index_html(tmp_path):
    sample = make_demo_sample()
    index = save_report_dir(tmp_path, sample, report_id="TEST_001", affected_only=False)

    assert index.exists()
    assert index.name == "index.html"

    html = index.read_text(encoding="utf-8")
    assert "<html" in html.lower()

    # manifest + per-chrom CNV payloads were emitted for the affected chromosomes
    assert (tmp_path / "data" / "manifest.js").exists()
    assert (tmp_path / "data" / "cnv" / "chr21.js").exists()


def test_demo_sample_shape():
    s = make_demo_sample()
    assert s.id == "DEMO_001"
    assert {e.chr for e in s.events} == {"21", "X", "5", "17"}
