"""
karyotype_viewer CLI entry point.

실제 파이프라인 output 기준:
  cnv.tsv         — 전체 genome bin (chrom, start, end, bin_BAF, copy_number_signal)
  nipt_results.tsv — syndrome 판정 결과 (SYNDROME, DIAGNOSIS, OBSERVED_CN, ...)

Usage
-----
# Web viewer
python -m karyotype_viewer \\
    --cnv     results/cnv.tsv \\
    --nipt    results/nipt_results.tsv \\
    --sample  results/sample_info.tsv   # optional

# Report 생성
python -m karyotype_viewer \\
    --cnv     results/cnv.tsv \\
    --nipt    results/nipt_results.tsv \\
    --report  output/SAMPLE_001

# 구버전 호환 (cnv 디렉토리 + marker DB)
python -m karyotype_viewer \\
    --cnv-dir  data/DEMO_NIPT/ \\
    --markers  data/nipt_markers.tsv \\
    --sample   data/DEMO_NIPT/sample.tsv
"""

from __future__ import annotations
import argparse
import sys
from pathlib import Path
import uvicorn


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="karyotype_viewer",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # ── 실제 파이프라인 output ────────────────────────────────────────────
    p.add_argument("--cnv",    metavar="TSV",  default=None,
                   help="전체 genome CNV TSV (chrom/start/end/bin_BAF/copy_number_signal)")
    p.add_argument("--nipt",   metavar="TSV",  default=None,
                   help="NIPT 판정 결과 TSV (SYNDROME/DIAGNOSIS/OBSERVED_CN/...)")
    p.add_argument("--markers", metavar="TSV", default=None,
               help="nipt_markers.tsv — feature 좌표 annotation용")
    # ── 구버전 호환 ───────────────────────────────────────────────────────
    p.add_argument("--cnv-dir",  metavar="DIR",  default=None, dest="cnv_dir",
                   help="[구버전] cnv_chr{N}.tsv 파일 디렉토리")
    p.add_argument("--sample",   metavar="TSV",  default=None,
                   help="Sample 메타데이터 TSV (선택)")
    # ── 공통 ─────────────────────────────────────────────────────────────
    p.add_argument("--report",    metavar="DIR", default=None,
                   help="지정 시 report 디렉토리 생성 후 종료")
    p.add_argument("--report-id", metavar="ID",  default=None, dest="report_id",
                   help="Report ID (기본: sample ID 또는 날짜)")
    p.add_argument("--chrom",     default="21",
                   help="초기 표시 염색체 (default: 21)")
    p.add_argument("--host",      default="0.0.0.0")
    p.add_argument("--port",      default=8050, type=int)

    args = p.parse_args(argv)

    # ── 입력 모드 결정 ────────────────────────────────────────────────────
    use_pipeline = bool(args.cnv or args.nipt)
    use_legacy   = bool(args.cnv_dir or args.markers)

    if not use_pipeline and not use_legacy:
        # 데모 모드
        print("[kv] 입력 파일 없음 — 내장 데모 샘플 사용")

    # ── CNV 로드 ──────────────────────────────────────────────────────────
    cnv_data   = {}
    syndromes  = {}
    sample     = None

    if use_pipeline and args.cnv:
        from .core.nipt_loader import load_cnv_tsv
        try:
            cnv_data = load_cnv_tsv(args.cnv)
        except Exception as e:
            print(f"[kv] ERROR cnv TSV: {e}", file=sys.stderr); sys.exit(1)

    elif use_legacy and args.cnv_dir:
        from .core.cnv_loader import load_cnv_dir
        try:
            cnv_data = load_cnv_dir(args.cnv_dir)
            print(f"[kv] CNV dir: {len(cnv_data)} chromosomes")
        except Exception as e:
            print(f"[kv] WARN cnv_dir: {e}")

    # ── Syndrome 판정 로드 ────────────────────────────────────────────────
    if use_pipeline and args.nipt:
        from .core.nipt_loader import load_nipt_results, build_events_from_results
        try:
            syndromes = load_nipt_results(args.nipt, marker_path=args.markers)
        except Exception as e:
            print(f"[kv] ERROR nipt TSV: {e}", file=sys.stderr); sys.exit(1)
    
        print(syndromes)

    elif use_legacy and args.markers:
        from .core.nipt_markers import load_marker_tsv
        try:
            syndromes = load_marker_tsv(args.markers)
            if cnv_data:
                _assign_calls_legacy(syndromes, cnv_data)
            print(f"[kv] markers: {len(syndromes)} syndromes")
        except Exception as e:
            print(f"[kv] WARN markers: {e}")

    # ── Sample 메타데이터 ─────────────────────────────────────────────────
    if args.sample:
        from .core.parser import load_sample_tsv
        try:
            sample = load_sample_tsv(args.sample)
            print(f"[kv] sample: {sample.id}")
        except Exception as e:
            print(f"[kv] WARN sample: {e}")

    # sample 없으면 nipt_results / cnv에서 자동 생성
    if sample is None:
        sample = _make_sample(syndromes, cnv_data, args)

    # nipt_results에서 events 자동 생성
    if use_pipeline and args.nipt and not sample.events:
        from .core.nipt_loader import build_events_from_results
        sample.events = build_events_from_results(syndromes)

    print(f"[kv] Sample: {sample.id}  events={len(sample.events)}  syndromes={len(syndromes)}")

    # ── Report mode ───────────────────────────────────────────────────────
    if args.report:
        from .reporter import save_report_dir
        import datetime
        rid = args.report_id or f"{sample.id}_{datetime.date.today().isoformat()}"
        out = save_report_dir(
            args.report, sample,
            syndromes=syndromes,
            cnv_data=cnv_data,
            report_id=rid,
        )
        print(f"[kv] Report → {out.parent}")
        return

    # ── Web viewer ────────────────────────────────────────────────────────
    chrom = args.chrom.replace("chr", "").upper()
    if chrom not in sample.display_chroms:
        chrom = sample.display_chroms[0]

    from .server import create_app
    app = create_app(
        sample, initial_chrom=chrom,
        title=f"Karyotype – {sample.id}",
        cnv_data=cnv_data,
        syndromes=syndromes,
    )
    print(f"[kv] → http://{args.host}:{args.port}/")
    uvicorn.run(app, host=args.host, port=args.port)


# ---------------------------------------------------------------------------
# Sample 자동 생성 (nipt_results에서)
# ---------------------------------------------------------------------------
def _make_sample(syndromes, cnv_data, args):
    from .core.models import SampleData
    import datetime

    # sex 추정: chrY CN 있으면 male
    sex = "female"
    if "Y" in cnv_data:
        df_y = cnv_data["Y"]
        cn_col = next((c for c in ["cn","copy_number_signal"] if c in df_y.columns), None)
        if cn_col:
            y_cn = float(df_y[cn_col].median())
            sex = "male" if y_cn > 0.5 else "female"

    sample_id = (
        args.report_id
        or (str(Path(args.cnv).stem).split('.')[0] if args.cnv else None)
        or (str(Path(args.nipt).stem).split('.')[0] if args.nipt else None)
        or f"SAMPLE_{datetime.date.today().isoformat()}"
    )

    return SampleData(
        id  = sample_id,
        sex = sex,
    )


# ---------------------------------------------------------------------------
# 구버전 호환 call 판정
# ---------------------------------------------------------------------------
def _assign_calls_legacy(syndromes, cnv_data):
    """marker DB + cnv_dir 방식 (구버전 호환)."""
    from pathlib import Path

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

        # cn 컬럼 탐색
        cn_col = next(
            (c for c in ["cn", "copy_number_signal", "ratio"] if c in df.columns),
            None
        )
        if cn_col is None:
            continue

        primary = next(
            (f for f in syn.features
             if f.feature_type in ("TargetChromosome", "PrimaryTargetRegion")),
            syn.features[0] if syn.features else None,
        )
        if primary is None:
            continue

        pos_col = "pos" if "pos" in df.columns else "start"
        mask = (df[pos_col] >= primary.start) & (df[pos_col] <= primary.end)
        if mask.sum() < 3:
            continue

        cn_med = float(df.loc[mask, cn_col].median())
        syn.cn_value = cn_med

        rule = RULES.get(syn.group, {})
        if syn.group == "Micro Deletion":
            if cn_med <= rule.get("pos_max", 1.5):
                syn.call = "HIGH_RISK"
            elif cn_med <= rule.get("sus_max", 1.65):
                syn.call = "SUSPECTED"
            else:
                syn.call = "LOW_RISK"
        else:
            if cn_med >= rule.get("pos_min", 2.5):
                syn.call = "HIGH_RISK"
            elif cn_med >= rule.get("sus_min", 2.3):
                syn.call = "SUSPECTED"
            elif "pos_max" in rule and cn_med <= rule["pos_max"]:
                syn.call = "HIGH_RISK"
            elif "sus_max" in rule and cn_med <= rule["sus_max"]:
                syn.call = "SUSPECTED"
            else:
                syn.call = "LOW_RISK"


if __name__ == "__main__":
    main()
