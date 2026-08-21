#!/usr/bin/env python3
"""
scripts/run_wps.py

통합 WPS 파이프라인
  step 1. compute   : BAM → chromosome별 NPY (multiprocess)
  step 2. aggregate : marker BED → aggregate profile parquet/npy  [선택]
  step 3. plot      : genome-wide line plot + marker aggregate plot

Usage:
    python run_wps.py \
        --bam    sample.bam \
        --outdir output/sample_id/wps \
        --sample sample_id \
        --mode   L \
        [--mode  S]              # 복수 지정 가능
        [--marker-bed markers.bed]
        [--jobs  8]
        [--resume]
        [--skip-compute]
        [--no-plot]

출력 구조:
    output/sample_id/wps/
    ├── L/
    │   ├── sample_id.chr1.raw.npy
    │   ├── sample_id.chr1.cov.npy
    │   ├── sample_id.chr1.frag_cov.npy
    │   ├── sample_id.chr1.norm.npy
    │   └── manifest.json
    ├── S/  (--mode S 지정 시)
    ├── plots/
    │   ├── sample_id.L.genome_wide.png
    │   ├── sample_id.S.genome_wide.png
    │   ├── sample_id.L.marker_aggregate.png
    │   └── sample_id.S.marker_aggregate.png
    └── aggregate/
        ├── sample_id.L.marker_stats.parquet
        ├── sample_id.L.marker_stats_profiles.npy
        ├── sample_id.S.marker_stats.parquet
        └── sample_id.S.marker_stats_profiles.npy
"""
import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from wps.compute   import run as wps_compute
from wps.aggregate import run as wps_aggregate
from wps.plotting  import plot_genome_wide, plot_marker_aggregate

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="NIPT WPS 통합 파이프라인")
    p.add_argument("--bam",        required=True,  help="입력 BAM (sorted, indexed)")
    p.add_argument("--outdir",     required=True,  help="출력 루트 디렉토리")
    p.add_argument("--sample",     required=True,  help="샘플 ID (파일명 prefix)")
    p.add_argument("--mode",       nargs="+", choices=["L", "S"], default=["L", "S"],
                   help="WPS mode (L / S 복수 지정 가능, default: L, S)")
    p.add_argument("--marker-bed", default=None,   help="marker BED (aggregate 선택)")
    p.add_argument("--jobs",       type=int, default=8, help="병렬 프로세스 수")
    p.add_argument("--min-mapq",   type=int, default=20, dest="min_mapq")
    p.add_argument("--win-size",   type=int, default=1000, dest="win_size",
                   help="blockMedian window 크기 (bp)")
    p.add_argument("--flank",      type=int, default=1000,
                   help="marker aggregate ±flank bp")
    p.add_argument("--resume",     action="store_true", help="기존 NPY 재사용")
    p.add_argument("--skip-compute", action="store_true", help="compute step 스킵")
    p.add_argument("--no-plot",    action="store_true", help="plot 생성 안 함")
    return p.parse_args()


def main() -> None:
    args    = parse_args()
    outdir  = Path(args.outdir)
    sid     = args.sample
    plot_dir = outdir  / "plots"
    agg_dir  = outdir / "aggregate"

    for mode in args.mode:
        logger.info("=" * 60)
        logger.info("MODE: %s", mode)
        logger.info("=" * 60)

        # ── step 1. compute ───────────────────────────────────────────
        if args.skip_compute:
            logger.info("[SKIP] compute — 기존 NPY 재사용")
        else:
            logger.info("=== STEP 1 : compute (multiprocess, jobs=%d) ===", args.jobs)
            wps_compute(
                bam_path = args.bam,
                out_dir  = str(outdir),
                mode     = mode,
                min_mapq = args.min_mapq,
                win_size = args.win_size,
                n_jobs   = args.jobs,
                resume   = args.resume,
                prefix   = sid,
            )

        # ── step 2. aggregate (marker BED 있을 때만) ──────────────────
        agg_parquet  = agg_dir / f"{sid}.{mode}.marker_stats.parquet"
        agg_profiles = agg_dir / f"{sid}.{mode}.marker_stats_profiles.npy"

        if args.marker_bed:
            logger.info("=== STEP 2 : aggregate [%s] ===", args.marker_bed)
            agg_dir.mkdir(parents=True, exist_ok=True)
            wps_aggregate(
                marker_bed    = args.marker_bed,
                out_path      = str(agg_parquet),
                wps_dir       = str(outdir),
                profile_flank = args.flank,
                prefix        = sid,
                save_profiles = True,
            )
        else:
            logger.info("[SKIP] aggregate — marker BED 미지정")

        # ── step 3. plot ──────────────────────────────────────────────
        if args.no_plot:
            logger.info("[SKIP] plot")
            continue

        logger.info("=== STEP 3 : plot ===")
        plot_dir.mkdir(parents=True, exist_ok=True)

        # genome-wide line plot
        gw_png = plot_dir / f"{sid}.{mode}.genome_wide.png"
        plot_genome_wide(
            wps_dir     = str(outdir),
            output_path = str(gw_png),
            mode        = mode,
            sample_id   = sid,
            prefix      = sid,
        )
        logger.info("genome-wide plot → %s", gw_png)

        # marker aggregate plot
        if args.marker_bed and agg_profiles.is_file():
            agg_png = plot_dir / f"{sid}.{mode}.marker_aggregate.png"
            plot_marker_aggregate(
                profiles_npy = str(agg_profiles),
                output_path  = str(agg_png),
                sample_id    = sid,
                title_suffix = f"mode={mode}",
            )
            logger.info("aggregate plot → %s", agg_png)

    logger.info("=== WPS 완료 ===")


if __name__ == "__main__":
    main()