"""
karyotype_viewer CLI entry point.

Usage
-----
# 내장 데모
python -m karyotype_viewer

# 실제 샘플 + CNV 디렉토리 + 마커 DB
python -m karyotype_viewer \\
    --sample data/DEMO_NIPT/sample.tsv \\
    --cnv    data/DEMO_NIPT/ \\
    --markers data/nipt_markers.tsv \\
    --chrom 21
"""

from __future__ import annotations
import argparse
import sys
import uvicorn


def main(argv=None):
    p = argparse.ArgumentParser(prog="karyotype_viewer")
    p.add_argument("--sample",  metavar="TSV",  default=None,
                   help="Sample metadata TSV. 생략 시 내장 데모 사용.")
    p.add_argument("--cnv",     metavar="DIR",  default=None,
                   help="cnv_chr{N}.tsv 파일이 있는 디렉토리.")
    p.add_argument("--markers", metavar="TSV",  default=None,
                   help="NIPT marker DB TSV (nipt_markers.tsv).")
    p.add_argument("--chrom",   default="21",
                   help="초기 표시 염색체 (default: 21).")
    p.add_argument("--host",    default="0.0.0.0")
    p.add_argument("--port",    default=8050, type=int)
    p.add_argument("--reload",  action="store_true")
    args = p.parse_args(argv)

    # ── Sample ──────────────────────────────────────────────────────────
    if args.sample:
        from .core.parser import load_sample_tsv
        try:
            sample = load_sample_tsv(args.sample)
            print(f"[kv] sample: {sample.id}  "
                  f"({len(sample.events)} events, {len(sample.genes)} genes)")
        except Exception as e:
            print(f"[kv] ERROR sample TSV: {e}", file=sys.stderr); sys.exit(1)
    else:
        from .demo_sample import make_demo_sample
        sample = make_demo_sample()
        print("[kv] 내장 데모 샘플 사용 (DEMO_001)")

    # ── CNV data ─────────────────────────────────────────────────────────
    cnv_data = {}
    if args.cnv:
        from .core.cnv_loader import load_cnv_dir
        try:
            cnv_data = load_cnv_dir(args.cnv)
            print(f"[kv] CNV loaded: {len(cnv_data)} chromosomes")
        except Exception as e:
            print(f"[kv] WARN cnv: {e}")

    # ── Marker DB + call 계산 ────────────────────────────────────────────
    syndromes = {}
    if args.markers:
        from .core.nipt_markers import load_marker_tsv
        try:
            syndromes = load_marker_tsv(args.markers)
            # CN 값으로 call 판정 (cnv_data 있을 때만)
            if cnv_data:
                _assign_calls(syndromes, cnv_data)
            print(f"[kv] markers: {len(syndromes)} syndromes")
        except Exception as e:
            print(f"[kv] WARN markers: {e}")

    # ── Initial chrom ─────────────────────────────────────────────────────
    chrom = args.chrom.replace("chr", "").upper()
    if chrom not in sample.display_chroms:
        chrom = sample.display_chroms[0]

    # ── Launch ────────────────────────────────────────────────────────────
    from .server import create_app
    app = create_app(
        sample, initial_chrom=chrom,
        title=f"Karyotype – {sample.id}",
        cnv_data=cnv_data,
        syndromes=syndromes,
    )
    print(f"[kv] → http://{args.host}:{args.port}/")
    uvicorn.run(app, host=args.host, port=args.port, reload=args.reload)


def _assign_calls(syndromes, cnv_data):
    """
    CNV DataFrame에서 syndrome primary region의 median CN을 계산해
    간단한 룰로 ABNORMAL / SUSPICIOUS / NORMAL 판정.
    """
    import numpy as np
    from .core.nipt_markers import NiptSyndrome

    RULES = {
        "Autosome Abnormality":       {"pos_min": 2.5, "sus_min": 2.3},
        "Sex Chromosome Abnormality": {"pos_min": 2.5, "sus_min": 2.3,
                                       "pos_max": 1.4, "sus_max": 1.6},
        "Micro Deletion":             {"pos_max": 1.5, "sus_max": 1.65},
    }

    for syn in syndromes.values():
        chrom = syn.primary_chrom
        if chrom not in cnv_data:
            continue
        df = cnv_data[chrom]

        # primary feature region
        primary = next(
            (f for f in syn.features
             if f.feature_type in ("TargetChromosome", "PrimaryTargetRegion")),
            syn.features[0] if syn.features else None,
        )
        if primary is None:
            continue

        mask = (df["pos"] >= primary.start) & (df["pos"] <= primary.end)
        if mask.sum() < 3:
            continue

        cn_med = float(df.loc[mask, "cn"].median())
        syn.cn_value = cn_med

        rule = RULES.get(syn.group, {})
        if syn.group == "Micro Deletion":
            if cn_med <= rule.get("pos_max", 1.5):
                syn.call = "ABNORMAL"
            elif cn_med <= rule.get("sus_max", 1.65):
                syn.call = "SUSPICIOUS"
            else:
                syn.call = "NORMAL"
        elif syn.group in ("Autosome Abnormality", "Sex Chromosome Abnormality"):
            if cn_med >= rule.get("pos_min", 2.5):
                syn.call = "ABNORMAL"
            elif cn_med >= rule.get("sus_min", 2.3):
                syn.call = "SUSPICIOUS"
            elif "pos_max" in rule and cn_med <= rule["pos_max"]:
                syn.call = "ABNORMAL"
            elif "sus_max" in rule and cn_med <= rule["sus_max"]:
                syn.call = "SUSPICIOUS"
            else:
                syn.call = "NORMAL"


if __name__ == "__main__":
    main()
