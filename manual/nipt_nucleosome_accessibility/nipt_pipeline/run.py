#!/usr/bin/env python3
"""
run.py  —  NIPT 통합 파이프라인 진입점

subcommand
----------
  cnv   CNV + Fetal Fraction
  wps   WPS compute + aggregate + plot
  all   cnv → wps 순서 전체

skip 정책 (--force 없으면 결과물 존재 시 자동 skip)
---------------------------------------------------
  bincount.tsv    존재 → binning  skip
  normalized.tsv  존재 → normalize skip
  cnv.tsv         존재 → cnv call skip
  ff json         존재 → fetal fraction skip
  wps manifest    존재 → compute  skip  (mode별)
  agg parquet     존재 → aggregate skip (mode별)

--force 지정 시 모든 단계 재실행
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from cnv.binning       import run_binning
from cnv.normalize     import run_normalize
from cnv.calculate_cnv import run_calculate_cnv
from cnv.region_cnv    import run_region_cnv
from cnv.plotting      import plot_cnv_scatter, plot_cnv_scatter_by_chrom
import fetal_fraction
from wps.compute       import run as wps_compute
from wps.aggregate     import run as wps_aggregate
from wps.plotting      import plot_genome_wide, plot_genome_wide_by_chrom, plot_marker_aggregate, plot_group_overlay
from wps.peak_calling  import run as wps_call_peaks

import yaml

logging.basicConfig(
    level    = logging.INFO,
    format   = "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    handlers = [logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)


# ── 유틸 ──────────────────────────────────────────────────────────────
def _load_config(path: str | None) -> dict:
    default  = Path(__file__).resolve().parent / "config" / "default.yaml"
    cfg_path = Path(path) if path else default
    with open(cfg_path) as fh:
        return yaml.safe_load(fh)


def _skip(path: Path, label: str, force: bool) -> bool:
    """결과물이 이미 있고 force 아니면 True(skip) 반환"""
    if not force and path.exists():
        logger.info("[SKIP] %s — 결과물 존재: %s", label, path)
        return True
    return False


# ── argparse 헬퍼 ─────────────────────────────────────────────────────
def _add_common(p: argparse.ArgumentParser) -> None:
    p.add_argument("--bam",    required=True, help="입력 BAM (sorted, indexed)")
    p.add_argument("--outdir", required=True, help="출력 루트 디렉토리")
    p.add_argument("--sample", required=True, help="샘플 ID (파일명 prefix)")
    p.add_argument("--config", default=None,  help="YAML config 파일")
    p.add_argument("--force",  action="store_true",
                   help="결과물 존재해도 모든 단계 재실행")


def _add_cnv_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--bw",       required=True, help="mappability bigWig")
    p.add_argument("--fasta",    required=True, help="hg38 reference FASTA")
    p.add_argument("--cnv-jobs", type=int, default=4, dest="cnv_jobs",
                   help="binning 병렬 프로세스 수 (default: 4)")
    p.add_argument("--no-plot",  action="store_true", help="CNV scatter plot 생성 안 함")
    p.add_argument("--no-ff",    action="store_true", help="Fetal Fraction 추정 스킵")


def _add_wps_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--mode",       nargs="+", choices=["L", "S"], default=["L"])
    p.add_argument("--marker-bed", default=None, dest="marker_bed",
                   help="marker BED (aggregate 선택)")
    p.add_argument("--jobs",     type=int, default=4)
    p.add_argument("--min-mapq", type=int, default=20,   dest="min_mapq")
    p.add_argument("--win-size", type=int, default=1000, dest="win_size")
    p.add_argument("--flank",    type=int, default=1000,
                   help="marker aggregate ±flank bp")
    p.add_argument("--no-wps-plot", action="store_true", dest="no_wps_plot",
                   help="WPS plot 생성 안 함")


# ── CNV 파이프라인 ─────────────────────────────────────────────────────
def _run_cnv(args: argparse.Namespace, cfg: dict) -> None:
    outdir   = Path(args.outdir)
    plot_dir = outdir / "plots"
    sid      = args.sample
    force    = args.force
    b, n, c  = cfg["binning"], cfg["normalization"], cfg["cnv"]

    outdir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    bincount_tsv   = outdir / f"{sid}.bincount.tsv"
    normalized_tsv = outdir / f"{sid}.normalized.tsv"
    cnv_tsv        = outdir / f"{sid}.cnv.tsv"
    ff_json        = outdir / f"{sid}.fetal_fraction.json"
    plot_short     = plot_dir / f"{sid}.cnv_short.png"
    plot_long      = plot_dir / f"{sid}.cnv_long.png"
    plot_combined  = plot_dir / f"{sid}.cnv_combined.png"
    plot_ff_corr   = plot_dir / f"{sid}.ff_correlation.png"

    # step 1. binning
    if not _skip(bincount_tsv, "binning", force):
        logger.info("=== CNV 1/3 : binning ===")
        run_binning(
            bam_path              = args.bam,
            bw_path               = args.bw,
            fasta_path            = args.fasta,
            output_path           = str(bincount_tsv),
            bin_size              = b["bin_size"],
            short_long_cutoff     = b["short_long_cutoff"],
            mapq_min              = b["mapq_min"],
            mappability_threshold = b["mappability_threshold"],
            n_jobs                = getattr(args, "cnv_jobs", 4),
        )

    # step 2. normalize
    if not _skip(normalized_tsv, "normalize", force):
        logger.info("=== CNV 2/3 : normalize ===")
        run_normalize(
            bin_tsv_path    = str(bincount_tsv),
            output_path     = str(normalized_tsv),
            excluded_chroms = set(n["excluded_chroms"]),
            gc_correction   = n["gc_correction"],
            gc_loess_span   = n["gc_loess_span"],
        )

    # step 3. cnv call (염색체 단위)
    if not _skip(cnv_tsv, "cnv call", force):
        logger.info("=== CNV 3/3 : CNV call ===")
        run_calculate_cnv(
            norm_tsv_path         = str(normalized_tsv),
            output_path           = str(cnv_tsv),
            zscore_gain_threshold = c["zscore_gain_threshold"],
            zscore_loss_threshold = c["zscore_loss_threshold"],
            min_bins_per_chrom    = c["min_bins_per_chrom"],
        )

    # step 3-1. region CNV call (PELT segmentation, microdeletion 검출)
    region_cnv_tsv = outdir / f"{sid}.region_cnv.tsv"
    if not _skip(region_cnv_tsv, "region cnv call", force):
        logger.info("=== CNV : Region CNV call (PELT) ===")
        # ff_percent: fetal_fraction.json이 있으면 로드
        ff_percent = None
        if ff_json.exists():
            import json
            with open(ff_json) as fh:
                ff_data = json.load(fh)
            ff_percent = ff_data.get("ff_seqff") or ff_data.get("ff_chry")
        run_region_cnv(
            norm_tsv_path         = str(normalized_tsv),
            output_path           = str(region_cnv_tsv),
            zscore_gain_threshold = c["zscore_gain_threshold"],
            zscore_loss_threshold = c["zscore_loss_threshold"],
            min_size_mb           = c.get("min_region_size_mb", 1.0),
            penalty_multiplier    = c.get("pelt_penalty_multiplier", 2.0),
            ff_percent            = ff_percent,
        )

    # step 4. plot
    if not getattr(args, "no_plot", False):
        # genome-wide scatter
        if not _skip(plot_short, "cnv scatter plot", force):
            logger.info("=== CNV : genome-wide scatter plot ===")
            plot_cnv_scatter(
                norm_tsv_path         = str(normalized_tsv),
                cnv_tsv_path          = str(cnv_tsv),
                output_short          = str(plot_short),
                output_long           = str(plot_long),
                output_combined       = str(plot_combined),
                sample_id             = sid,
                zscore_gain_threshold = c["zscore_gain_threshold"],
                zscore_loss_threshold = c["zscore_loss_threshold"],
            )

        # chromosome별 scatter
        chrom_plot_dir = plot_dir / "by_chrom"
        if not force and chrom_plot_dir.exists() and any(chrom_plot_dir.iterdir()):
            logger.info("[SKIP] CNV chrom plot — 결과물 존재: %s", chrom_plot_dir)
        else:
            logger.info("=== CNV : chromosome별 scatter plot ===")
            plot_cnv_scatter_by_chrom(
                norm_tsv_path         = str(normalized_tsv),
                cnv_tsv_path          = str(cnv_tsv),
                output_dir            = str(chrom_plot_dir),
                sample_id             = sid,
                zscore_gain_threshold = c["zscore_gain_threshold"],
                zscore_loss_threshold = c["zscore_loss_threshold"],
                region_cnv_path       = str(region_cnv_tsv) if region_cnv_tsv.exists() else None,
            )

    # step 5. fetal fraction
    if not getattr(args, "no_ff", False):
        if not _skip(ff_json, "fetal fraction", force):
            logger.info("=== CNV : Fetal Fraction ===")
            fetal_fraction.run(
                bincount_tsv   = str(bincount_tsv),
                normalized_tsv = str(normalized_tsv),
                output_json    = str(ff_json),
                sample_id      = sid,
                plot_path      = str(plot_ff_corr),
            )

    _log_outputs("CNV", {
        "bincount":       bincount_tsv,
        "normalized":     normalized_tsv,
        "cnv call":       cnv_tsv,
        "region cnv":     region_cnv_tsv,
        "ff json":        ff_json,
        "plot short":     plot_short,
        "plot long":      plot_long,
        "plot combined":  plot_combined,
        "ff correlation": plot_ff_corr,
    })


# ── WPS 파이프라인 ─────────────────────────────────────────────────────
def _run_wps(args: argparse.Namespace) -> None:
    wps_root = Path(args.outdir) / "wps"
    sid      = args.sample
    force    = args.force
    modes    = getattr(args, "mode", ["L"])

    for mode in modes:
        logger.info("=" * 55)
        logger.info("WPS MODE: %s", mode)
        logger.info("=" * 55)

        plot_dir     = wps_root / "plots"
        agg_dir      = wps_root / "aggregate"
        manifest_path = wps_root / mode / "manifest.json"
        agg_parquet   = agg_dir / f"{sid}.{mode}.marker_stats.parquet"
        agg_profiles  = agg_dir / f"{sid}.{mode}.marker_stats_profiles.npy"

        # step 1. compute
        if not _skip(manifest_path, f"WPS compute [{mode}]", force):
            logger.info("=== WPS compute [%s] (jobs=%d) ===", mode, args.jobs)
            wps_compute(
                bam_path = args.bam,
                out_dir  = str(wps_root),
                mode     = mode,
                min_mapq = args.min_mapq,
                win_size = args.win_size,
                n_jobs   = args.jobs,
                resume   = False,          # skip은 manifest 기준으로 이미 처리
                prefix   = sid,
            )

        # step 2. aggregate
        if getattr(args, "marker_bed", None):
            if not _skip(agg_parquet, f"WPS aggregate [{mode}]", force):
                logger.info("=== WPS aggregate [%s] ===", mode)
                agg_dir.mkdir(parents=True, exist_ok=True)
                wps_aggregate(
                    marker_bed    = args.marker_bed,
                    out_path      = str(agg_parquet),
                    wps_dir       = str(wps_root),
                    profile_flank = args.flank,
                    prefix        = sid,
                    save_profiles = True,
                )
        else:
            logger.info("[SKIP] WPS aggregate — marker BED 미지정")

        # step 3. plot
        if not getattr(args, "no_wps_plot", False):
            plot_dir.mkdir(parents=True, exist_ok=True)
            agg_png = plot_dir / f"{sid}.{mode}.marker_aggregate.png"

            if getattr(args, "marker_bed", None) and agg_profiles.is_file():
                if not _skip(agg_png, f"WPS aggregate plot [{mode}]", force):
                    logger.info("=== WPS aggregate plot [%s] ===", mode)
                    plot_marker_aggregate(
                        profiles_npy = str(agg_profiles),
                        output_path  = str(agg_png),
                        sample_id    = sid,
                        title_suffix = f"mode={mode}",
                    )

            # group overlay plot (그룹별 개별 PNG)
            group_dir     = agg_dir / "by_group"
            group_plot_dir = plot_dir / f"{mode}_groups"
            if getattr(args, "marker_bed", None) and group_dir.is_dir():
                # 디렉토리 존재 + force 없으면 skip
                if not force and group_plot_dir.exists() and any(group_plot_dir.iterdir()):
                    logger.info("[SKIP] WPS group plot [%s] — 결과물 존재: %s",
                                mode, group_plot_dir)
                else:
                    logger.info("=== WPS group plot [%s] ===", mode)
                    plot_group_overlay(
                        group_dir     = str(group_dir),
                        output_dir    = str(group_plot_dir),
                        sample_id     = sid,
                        profile_flank = args.flank,
                    )

        # peak calling (L mode만)
        peaks_bed = wps_root / f"{sid}.{mode}.peaks.bed"
        if mode == "L":
            if not _skip(peaks_bed, f"WPS peak calling [{mode}]", force):
                logger.info("=== WPS peak calling [%s] ===", mode)
                from wps.peak_calling import run as wps_call_peaks
                wps_call_peaks(
                    wps_dir    = str(wps_root),
                    output_bed = str(peaks_bed),
                    prefix     = sid,
                    mode       = mode,
                )

        _log_outputs(f"WPS [{mode}]", {
            "wps dir":   wps_root / mode,
            "manifest":  manifest_path,
            "peaks bed": peaks_bed,
            "aggregate": agg_parquet,
        })

    # ── genome-wide plot: L+S 통합 (mode 루프 밖에서 1회) ────────────
    if not getattr(args, "no_wps_plot", False):
        plot_dir = wps_root / "plots"
        plot_dir.mkdir(parents=True, exist_ok=True)

        # genome-wide 통합 PNG
        gw_png = plot_dir / f"{sid}.genome_wide.png"
        if not _skip(gw_png, "WPS genome-wide plot (L+S)", force):
            logger.info("=== WPS genome-wide plot (L+S 통합) ===")
            plot_genome_wide(
                wps_dir     = str(wps_root),
                output_path = str(gw_png),
                sample_id   = sid,
                prefix      = sid,
            )

        # chromosome별 개별 PNG
        chrom_plot_dir = plot_dir / "by_chrom"
        if not force and chrom_plot_dir.exists() and any(chrom_plot_dir.iterdir()):
            logger.info("[SKIP] WPS chrom plot — 결과물 존재: %s", chrom_plot_dir)
        else:
            logger.info("=== WPS chromosome별 plot ===")
            plot_genome_wide_by_chrom(
                wps_dir    = str(wps_root),
                output_dir = str(chrom_plot_dir),
                sample_id  = sid,
                prefix     = sid,
            )

        _log_outputs("WPS plots", {
            "genome-wide": gw_png,
            "by_chrom":    chrom_plot_dir,
        })


# ── 출력 로그 ─────────────────────────────────────────────────────────
def _log_outputs(label: str, paths: dict) -> None:
    logger.info("=== %s 완료 ===", label)
    for k, v in paths.items():
        exists = "✓" if Path(v).exists() else "✗"
        logger.info("  %s %-18s: %s", exists, k, v)


# ── argparse ──────────────────────────────────────────────────────────
def build_parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(
        prog            = "run.py",
        description     = "NIPT 통합 파이프라인",
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
    )
    sub = root.add_subparsers(dest="command", required=True)

    p_cnv = sub.add_parser("cnv", help="CNV + Fetal Fraction")
    _add_common(p_cnv);  _add_cnv_args(p_cnv)

    p_wps = sub.add_parser("wps", help="WPS + aggregate + plot")
    _add_common(p_wps);  _add_wps_args(p_wps)

    p_all = sub.add_parser("all", help="CNV + FF + WPS 전체")
    _add_common(p_all);  _add_cnv_args(p_all);  _add_wps_args(p_all)

    return root


def main() -> None:
    args = build_parser().parse_args()
    cfg  = _load_config(args.config)

    if args.command == "cnv":
        _run_cnv(args, cfg)

    elif args.command == "wps":
        _run_wps(args)

    elif args.command == "all":
        logger.info("▶ CNV + FF")
        _run_cnv(args, cfg)
        logger.info("▶ WPS")
        _run_wps(args)

    logger.info("✓ 파이프라인 완료")


if __name__ == "__main__":
    main()