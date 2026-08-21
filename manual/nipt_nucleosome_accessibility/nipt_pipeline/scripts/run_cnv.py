#!/usr/bin/env python3
"""
scripts/run_cnv.py

통합 NIPT CNV 파이프라인
  step 1. binning    : BAM → bin count TSV
  step 2. normalize  : bin count → GC correction + Z-score TSV
  step 3. cnv call   : Z-score → 염색체 CNV call TSV
  step 4. plot       : short/long Z-score scatter plot (PNG)

Usage:
    python run_cnv.py \
        --bam    sample.bam \
        --bw     hg38.mappability.100mer.bw \
        --fasta  hg38.fa \
        --outdir output/sample_id \
        --sample sample_id \
        [--config config/default.yaml] \
        [--skip-binning]   # bincount TSV가 이미 있을 때
        [--skip-normalize] # normalized TSV가 이미 있을 때
        [--no-plot]
"""
import argparse
import logging
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cnv.binning       import run_binning
from cnv.normalize     import run_normalize
from cnv.calculate_cnv import run_calculate_cnv
from cnv.plotting      import plot_cnv_scatter
import fetal_fraction

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="NIPT CNV 통합 파이프라인")
    p.add_argument("--bam",    required=True,  help="입력 BAM (sorted, indexed)")
    p.add_argument("--bw",     required=True,  help="mappability bigWig")
    p.add_argument("--fasta",  required=True,  help="hg38 reference FASTA")
    p.add_argument("--outdir", required=True,  help="출력 디렉토리")
    p.add_argument("--sample", required=True,  help="샘플 ID (출력 파일명 prefix)")
    p.add_argument("--config", default=None,   help="YAML config 파일")
    p.add_argument("--skip-binning",   action="store_true", help="binning 스킵 (기존 TSV 재사용)")
    p.add_argument("--skip-normalize", action="store_true", help="normalization 스킵 (기존 TSV 재사용)")
    p.add_argument("--no-plot",        action="store_true", help="scatter plot 생성 안 함")
    p.add_argument("--no-ff",          action="store_true", help="fetal fraction 추정 스킵")
    return p.parse_args()


def load_config(path: str | None) -> dict:
    default_cfg = Path(__file__).resolve().parents[1] / "config" / "default.yaml"
    cfg_path = Path(path) if path else default_cfg
    with open(cfg_path) as fh:
        return yaml.safe_load(fh)


def main() -> None:
    args   = parse_args()
    cfg    = load_config(args.config)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    sid = args.sample
    b   = cfg["binning"]
    n   = cfg["normalization"]
    c   = cfg["cnv"]

    # 출력 파일 경로 정의
    bincount_tsv   = outdir / f"{sid}.bincount.tsv"
    normalized_tsv = outdir / f"{sid}.normalized.tsv"
    cnv_tsv        = outdir / f"{sid}.cnv.tsv"
    plot_short     = outdir / f"{sid}.cnv_short.png"
    plot_long      = outdir / f"{sid}.cnv_long.png"
    ff_json        = outdir / f"{sid}.fetal_fraction.json"

    # ── step 1. binning ──────────────────────────────────────────────────────
    if args.skip_binning:
        logger.info("[SKIP] binning → %s 재사용", bincount_tsv)
    else:
        logger.info("=== STEP 1 / 3 : binning ===")
        run_binning(
            bam_path              = args.bam,
            bw_path               = args.bw,
            fasta_path            = args.fasta,
            output_path           = str(bincount_tsv),
            bin_size              = b["bin_size"],
            short_long_cutoff     = b["short_long_cutoff"],
            mapq_min              = b["mapq_min"],
            mappability_threshold = b["mappability_threshold"],
        )

    # ── step 2. normalize ────────────────────────────────────────────────────
    if args.skip_normalize:
        logger.info("[SKIP] normalize → %s 재사용", normalized_tsv)
    else:
        logger.info("=== STEP 2 / 3 : normalize ===")
        run_normalize(
            bin_tsv_path    = str(bincount_tsv),
            output_path     = str(normalized_tsv),
            excluded_chroms = set(n["excluded_chroms"]),
            gc_correction   = n["gc_correction"],
            gc_loess_span   = n["gc_loess_span"],
        )

    # ── step 3. cnv call ─────────────────────────────────────────────────────
    logger.info("=== STEP 3 / 3 : CNV call ===")
    run_calculate_cnv(
        norm_tsv_path         = str(normalized_tsv),
        output_path           = str(cnv_tsv),
        zscore_gain_threshold = c["zscore_gain_threshold"],
        zscore_loss_threshold = c["zscore_loss_threshold"],
        min_bins_per_chrom    = c["min_bins_per_chrom"],
    )

    # ── step 4. plot ─────────────────────────────────────────────────────────
    if not args.no_plot:
        logger.info("=== PLOT : scatter plot 생성 ===")
        plot_cnv_scatter(
            norm_tsv_path         = str(normalized_tsv),
            cnv_tsv_path          = str(cnv_tsv),
            output_short          = str(plot_short),
            output_long           = str(plot_long),
            sample_id             = sid,
            zscore_gain_threshold = c["zscore_gain_threshold"],
            zscore_loss_threshold = c["zscore_loss_threshold"],
        )

    # ── step 5. fetal fraction ───────────────────────────────────────────
    if not args.no_ff:
        logger.info("=== STEP 4 / 4 : Fetal Fraction 추정 ===")
        fetal_fraction.run(
            bincount_tsv   = str(bincount_tsv),
            normalized_tsv = str(normalized_tsv),
            output_json    = str(ff_json),
            sample_id      = sid,
        )

    logger.info("=== 완료 ===")
    logger.info("  bincount   : %s", bincount_tsv)
    logger.info("  normalized : %s", normalized_tsv)
    logger.info("  cnv call   : %s", cnv_tsv)
    if not args.no_plot:
        logger.info("  plot short : %s", plot_short)
        logger.info("  plot long  : %s", plot_long)
    if not args.no_ff:
        logger.info("  fetal frac : %s", ff_json)


if __name__ == "__main__":
    main()