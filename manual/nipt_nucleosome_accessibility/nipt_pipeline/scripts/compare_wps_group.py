#!/usr/bin/env python3
"""
wps_group_compare.py
====================
파이프라인 NPY 중간파일(L/, S/ 서브폴더)을 직접 읽어
marker 영역별 L-WPS / S-WPS aggregate profile을 샘플×그룹별로 정리.

출력:
  {outdir}/
    {marker_type}.L.wps_compare.png   ← L mode 그룹 비교 plot
    {marker_type}.S.wps_compare.png   ← S mode 그룹 비교 plot
    {marker_type}.per_sample.tsv      ← 샘플별 L/S 통계값
    {marker_type}.group_summary.tsv   ← 그룹별 L/S 통계 요약

samples.tsv 형식 (헤더 필요, mode 컬럼 불필요):
  sample_id       group             wps_dir
  SRR19224105     cnv_pregnant      /path/to/SRR19224105/wps
  SRR5497625      normal_pregnant   /path/to/SRR5497625/wps

marker BED 형식:
  chrom  start  end  marker_id  marker_type  [flank]
  flank 컬럼은 --marker-flank-col 사용 시 활성화

flank 우선순위:
  1. --flank-map JSON (marker_type별)
  2. BED 6번째 컬럼 (--marker-flank-col 사용 시)
  3. --flank 전역값 (기본 1000)
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from scipy.signal import savgol_filter
import matplotlib.patches as mpatches

sys.path.insert(0, str(Path(__file__).resolve().parent))

logging.basicConfig(
    level  = logging.INFO,
    format = "%(asctime)s [%(levelname)s] %(message)s",
)
log = logging.getLogger(__name__)

# ── 상수 ─────────────────────────────────────────────────────────────
MODES   = ["L", "S"]
METRICS = {"L": "norm", "S": "norm"}   # mode별 기본 metric

METRIC_SUFFIX = {
    "raw":      "raw.npy",
    "cov":      "cov.npy",
    "frag_cov": "frag_cov.npy",
    "norm":     "norm.npy",
}

MODE_COLORS = {
    "L": {"normal_pregnant": "#1E88E5", "normal_nonpregnant": "#43A047",
          "cnv_pregnant": "#E53935", "disease_pregnant": "#FB8C00", "other": "#8E24AA"},
    "S": {"normal_pregnant": "#90CAF9", "normal_nonpregnant": "#A5D6A7",
          "cnv_pregnant": "#EF9A9A", "disease_pregnant": "#FFCC80", "other": "#CE93D8"},
}

def _get_color(mode: str, group: str) -> str:
    return MODE_COLORS.get(mode, MODE_COLORS["L"]).get(group, "#999999")


# ── NPY 탐색 ─────────────────────────────────────────────────────────
def _find_npy(
    wps_dir: str,
    prefix:  str,
    mode:    str,
    chrom:   str,
    metric:  str,
) -> Optional[Path]:
    suffix = METRIC_SUFFIX.get(metric, f"{metric}.npy")
    base   = Path(wps_dir)

    candidates = [
        base / mode / f"{prefix}.{chrom}.{suffix}",
        base / mode / f"{chrom}.{suffix}",
        base / f"{prefix}.{chrom}.{suffix}",
        base / f"{chrom}.{suffix}",
    ]

    # manifest.json 참조
    manifest_path = base / mode / "manifest.json"
    if manifest_path.exists():
        try:
            m = json.loads(manifest_path.read_text())
            p = m.get("chromosomes", {}).get(chrom, {}).get(metric, "")
            if p and Path(p).exists():
                return Path(p)
        except Exception:
            pass

    for p in candidates:
        if p.exists():
            return p
    return None


def _load_npy(path: Path) -> Optional[np.ndarray]:
    try:
        arr = np.load(path, mmap_mode="r", allow_pickle=False)
        return arr if arr.ndim == 1 else None
    except Exception as e:
        log.warning("NPY 로드 실패 %s: %s", path, e)
        return None


def load_chrom_track(
    wps_dir: str, prefix: str, mode: str, chrom: str, metric: str,
) -> Optional[np.ndarray]:
    p = _find_npy(wps_dir, prefix, mode, chrom, metric)
    if p is None:
        return None
    return _load_npy(p)


# ── BED 로더 ─────────────────────────────────────────────────────────
def load_marker_bed(bed_path: str, use_flank_col: bool = False) -> pd.DataFrame:
    df = pd.read_csv(bed_path, sep="\t", header=None, comment="#")

    base_cols = ["chrom", "start", "end", "marker_id", "marker_type"]
    extra = []
    if use_flank_col and df.shape[1] >= 6:
        extra.append("flank")
    if df.shape[1] >= len(base_cols) + len(extra) + 1:
        extra.append("strand")
    n_named = len(base_cols) + len(extra)
    df.columns = base_cols + extra + [f"_c{i}" for i in range(df.shape[1] - n_named)]

    df["start"] = pd.to_numeric(df["start"], errors="coerce").fillna(0).astype("int64")
    df["end"]   = pd.to_numeric(df["end"],   errors="coerce").fillna(0).astype("int64")
    df["marker_type"] = df["marker_type"].fillna("marker").astype(str)
    df["marker_id"]   = df["marker_id"].fillna(
        df["chrom"].astype(str) + ":" + df["start"].astype(str)
    ).astype(str)
    if "flank" in df.columns:
        df["flank"] = pd.to_numeric(df["flank"], errors="coerce").fillna(0).astype("int64")
    if "strand" not in df.columns:
        df["strand"] = "+"

    df = df[df["end"] > df["start"]].copy()
    log.info("BED 로드: %d markers, types=%s",
             len(df), df["marker_type"].unique().tolist())
    return df


# ── 샘플 aggregate ────────────────────────────────────────────────────
def aggregate_sample_mode(
    wps_dir:      str,
    prefix:       str,
    mode:         str,
    metric:       str,
    markers:      pd.DataFrame,
    global_flank: int,
    min_valid:    float = 0.7,
) -> Optional[np.ndarray]:
    """
    단일 샘플 × 단일 mode의 aggregate mean profile.
    chrom별 NPY를 캐싱해서 반복 로드 방지.
    """
    has_per_flank = "flank" in markers.columns
    target_len    = global_flank * 2

    # chrom별 캐시
    chrom_cache: dict[str, Optional[np.ndarray]] = {}
    for chrom in markers["chrom"].unique():
        chrom_cache[chrom] = load_chrom_track(wps_dir, prefix, mode, chrom, metric)

    profiles = []
    for _, mk in markers.iterrows():
        chrom  = mk["chrom"]
        arr    = chrom_cache.get(chrom)
        if arr is None:
            continue

        center  = (int(mk["start"]) + int(mk["end"])) // 2
        strand  = mk.get("strand", "+")
        m_flank = int(mk["flank"]) if has_per_flank and mk["flank"] > 0 else global_flank

        start = center - m_flank
        end   = center + m_flank
        lo, hi = max(0, start), min(len(arr), end)
        if lo >= hi:
            continue

        profile = np.full(end - start, np.nan, dtype=np.float32)
        profile[lo - start: lo - start + (hi - lo)] = arr[lo:hi]

        if strand == "-":
            profile = profile[::-1]

        if np.isfinite(profile).mean() < min_valid:
            continue

        # target_len으로 통일
        cur = len(profile)
        if cur == target_len:
            profiles.append(profile)
        elif cur > target_len:
            d = (cur - target_len) // 2
            profiles.append(profile[d: d + target_len])
        else:
            pl = (target_len - cur) // 2
            pr = target_len - cur - pl
            profiles.append(np.pad(profile, (pl, pr), constant_values=np.nan))

    if not profiles:
        return None

    mat = np.stack(profiles, axis=0)
    return np.nanmean(mat, axis=0).astype(np.float32)


# ── 통계값 계산 ───────────────────────────────────────────────────────
def compute_stats(profile: np.ndarray) -> dict:
    n = len(profile)
    c = n // 2
    v = np.isfinite(profile)

    if v.sum() < n * 0.5:
        return {k: np.nan for k in
                ("peak_amp","trough_depth","oscillation","period_score")}

    p = np.where(v, profile, 0.0)

    return {
        "peak_amp":     float(np.max(p[max(0,c-150): min(n,c+150)])),
        "trough_depth": float(np.min(p[max(0,c-50):  min(n,c+50)])),
        "oscillation":  float(np.std(p[v])),
        "period_score": float(
            np.abs(np.fft.rfft(p))[
                int(np.argmin(np.abs(np.fft.rfftfreq(n, 1.0) - 0.1)))
            ] / (np.abs(np.fft.rfft(p)).mean() + 1e-12)
        ),
    }


# ── 비교 plot (단일 mode) ─────────────────────────────────────────────
def plot_mode_compare(
    group_profiles: dict[str, list[np.ndarray]],
    group_labels:   dict[str, str],
    mode:           str,
    marker_type:    str,
    flank:          int,
    metric:         str,
    output_path:    str,
    smooth_window:  int  = 21,
    smooth_poly:    int  = 3,
) -> None:
    n_g = len(group_profiles)
    n_panels = 2 if n_g >= 2 else 1
    fig, axes = plt.subplots(n_panels, 1,
                             figsize=(14, 4.5 * n_panels), sharex=True)
    if n_panels == 1:
        axes = [axes]

    pos = np.arange(-flank, flank)
    group_means: dict[str, np.ndarray] = {}

    ax = axes[0]
    ax.set_facecolor("white")
    ax.axvline(0, color="#F44336", lw=0.9, ls="--", alpha=0.7, zorder=5)
    ax.axhline(0, color="#9E9E9E", lw=0.6, alpha=0.4)

    for group, profiles in group_profiles.items():
        if not profiles:
            continue
        color = _get_color(mode, group)
        label = group_labels.get(group, group)
        mat   = np.stack(profiles, axis=0)
        mean  = np.nanmean(mat, axis=0)
        ci    = 1.96 * stats.sem(mat, axis=0, nan_policy="omit")

        sw = min(smooth_window, len(mean) // 4 * 2 + 1)
        ms = savgol_filter(mean, sw, smooth_poly) if sw >= 3 else mean
        group_means[group] = ms

        ax.fill_between(pos, ms - ci, ms + ci, color=color, alpha=0.15)
        ax.plot(pos, ms, color=color, lw=1.8,
                label=f"{label} (n={len(profiles)})")

    ax.set_ylabel(f"{mode}-WPS  [{metric}]", fontsize=10)
    ax.set_title(
        f"{marker_type} — {mode}-WPS  (±{flank}bp, metric={metric})",
        fontsize=11, fontweight="bold",
    )
    ax.legend(fontsize=9, framealpha=0.85, loc="upper right")
    ax.spines[["top","right"]].set_visible(False)

    if n_panels == 2 and len(group_means) >= 2:
        ax2 = axes[1]
        ax2.set_facecolor("white")
        ax2.axvline(0, color="#F44336", lw=0.9, ls="--", alpha=0.7, zorder=5)
        ax2.axhline(0, color="#9E9E9E", lw=0.6, alpha=0.4)
        groups    = list(group_means.keys())
        ref_label = group_labels.get(groups[0], groups[0])
        for g in groups[1:]:
            color = _get_color(mode, g)
            label = group_labels.get(g, g)
            diff  = group_means[g] - group_means[groups[0]]
            ax2.plot(pos, diff, color=color, lw=1.5,
                     label=f"{label} − {ref_label}")
            ax2.fill_between(pos, 0, diff, color=color, alpha=0.12)
        ax2.set_ylabel("Δ WPS", fontsize=10)
        ax2.legend(fontsize=9, framealpha=0.85)
        ax2.spines[["top","right"]].set_visible(False)

    axes[-1].set_xlabel("Distance from marker center (bp)", fontsize=10)
    axes[-1].xaxis.set_major_formatter(
        plt.FuncFormatter(lambda x, _: f"{int(x):+d}" if x != 0 else "0")
    )
    axes[-1].set_xlim(-flank, flank)
    fig.tight_layout()
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info("저장: %s", output_path)


# ── L+S 합친 plot ────────────────────────────────────────────────────
def plot_ls_combined(
    l_group: dict[str, list[np.ndarray]],
    s_group: dict[str, list[np.ndarray]],
    group_labels: dict[str, str],
    marker_type:  str,
    flank:        int,
    output_path:  str,
    smooth_window: int = 21,
    smooth_poly:   int = 3,
) -> None:
    """
    L-WPS(위) / S-WPS(아래) 2행 × 그룹수 열 subplot.
    """
    all_groups = sorted(set(list(l_group.keys()) + list(s_group.keys())))
    n_g = len(all_groups)
    if n_g == 0:
        return

    fig, axes = plt.subplots(2, n_g, figsize=(6 * n_g, 8), sharex=True)
    if n_g == 1:
        axes = axes.reshape(2, 1)

    pos = np.arange(-flank, flank)
    fig.suptitle(f"{marker_type} — L/S WPS per group  (±{flank}bp)",
                 fontsize=12, fontweight="bold", y=1.01)

    for col, group in enumerate(all_groups):
        label = group_labels.get(group, group)

        for row, (mode, gp_dict) in enumerate([("L", l_group), ("S", s_group)]):
            ax = axes[row, col]
            ax.set_facecolor("white")
            ax.axvline(0, color="#F44336", lw=0.9, ls="--", alpha=0.7, zorder=5)
            ax.axhline(0, color="#9E9E9E", lw=0.6, alpha=0.4)

            profiles = gp_dict.get(group, [])
            if profiles:
                mat  = np.stack(profiles, axis=0)
                mean = np.nanmean(mat, axis=0)
                ci   = 1.96 * stats.sem(mat, axis=0, nan_policy="omit")
                sw = min(smooth_window, len(mean) // 4 * 2 + 1)
                ms = savgol_filter(mean, sw, smooth_poly) if sw >= 3 else mean
                color = _get_color(mode, group)
                ax.fill_between(pos, ms - ci, ms + ci, color=color, alpha=0.2)
                ax.plot(pos, ms, color=color, lw=1.8)
            else:
                ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                        ha="center", va="center", color="#9E9E9E")

            ax.set_ylabel(f"{mode}-WPS", fontsize=9)
            if row == 0:
                ax.set_title(f"{label}\n(n={len(profiles)})", fontsize=10)
            ax.spines[["top","right"]].set_visible(False)
            ax.tick_params(labelsize=8)

        axes[-1, col].set_xlabel("Distance from center (bp)", fontsize=9)
        axes[-1, col].xaxis.set_major_formatter(
            plt.FuncFormatter(lambda x, _: f"{int(x):+d}" if x != 0 else "0")
        )
        axes[-1, col].set_xlim(-flank, flank)

    fig.tight_layout()
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info("저장: %s", output_path)


# ── 메인 ─────────────────────────────────────────────────────────────
def main() -> None:
    ap = argparse.ArgumentParser(
        description="WPS NPY 기반 marker aggregate L/S 그룹 비교",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--samples",     required=True,
                    help="샘플 TSV (컬럼: sample_id, group, wps_dir)")
    ap.add_argument("--marker-bed",  required=True, dest="marker_bed")
    ap.add_argument("--outdir",      default="compare_plots")
    ap.add_argument("--metric-l",    default="norm", dest="metric_l",
                    choices=["raw","norm","frag_cov","cov"],
                    help="L mode 비교 metric")
    ap.add_argument("--metric-s",    default="norm", dest="metric_s",
                    choices=["raw","norm","frag_cov","cov"],
                    help="S mode 비교 metric")
    ap.add_argument("--flank",       type=int, default=1000)
    ap.add_argument("--flank-map",   default=None, dest="flank_map",
                    help="marker_type별 flank JSON (예: '{\"TSS\":500}')")
    ap.add_argument("--marker-flank-col", action="store_true", dest="marker_flank_col",
                    help="BED 6번째 컬럼을 marker별 flank로 사용")
    ap.add_argument("--group-labels", default=None, dest="group_labels",
                    help="그룹 표시명 JSON")
    ap.add_argument("--smooth-window", type=int, default=21, dest="smooth_window")
    ap.add_argument("--min-valid",   type=float, default=0.7, dest="min_valid")
    ap.add_argument("--no-combined", action="store_true", dest="no_combined",
                    help="L/S 합친 subplot 생략")
    args = ap.parse_args()

    # 샘플 목록
    samples_df = pd.read_csv(args.samples, sep="\t")
    for col in ("sample_id","group","wps_dir"):
        if col not in samples_df.columns:
            log.error("samples.tsv 필수 컬럼 없음: %s", col)
            sys.exit(1)

    # 표시명
    group_labels: dict[str, str] = {}
    if args.group_labels:
        try:
            group_labels = json.loads(args.group_labels)
        except json.JSONDecodeError as e:
            log.warning("group-labels 파싱 실패: %s", e)

    # flank_map
    flank_map: dict[str, int] = {}
    if args.flank_map:
        try:
            flank_map = {k: int(v) for k, v in json.loads(args.flank_map).items()}
        except Exception as e:
            log.warning("flank-map 파싱 실패: %s", e)

    # BED
    markers = load_marker_bed(args.marker_bed, use_flank_col=args.marker_flank_col)
    marker_types = markers["marker_type"].unique().tolist()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    for mtype in marker_types:
        eff_flank   = flank_map.get(mtype, args.flank)
        sub_markers = markers[markers["marker_type"] == mtype].reset_index(drop=True)
        log.info("=== %s (n=%d markers, flank=%d) ===", mtype, len(sub_markers), eff_flank)

        # 샘플별 L/S aggregate 계산
        per_sample_rows = []
        l_group: dict[str, list[np.ndarray]] = {}
        s_group: dict[str, list[np.ndarray]] = {}

        for _, row in samples_df.iterrows():
            sid     = str(row["sample_id"])
            group   = str(row["group"])
            wps_dir = str(row["wps_dir"])
            group_labels.setdefault(group, group)

            log.info("  [%s / %s]", sid, group)
            row_data = {
                "sample_id": sid,
                "group":     group,
                "marker_type": mtype,
                "flank":     eff_flank,
            }

            for mode, metric in [("L", args.metric_l), ("S", args.metric_s)]:
                profile = aggregate_sample_mode(
                    wps_dir      = wps_dir,
                    prefix       = sid,
                    mode         = mode,
                    metric       = metric,
                    markers      = sub_markers,
                    global_flank = eff_flank,
                    min_valid    = args.min_valid,
                )
                if profile is None:
                    log.warning("    [%s / %s] profile 없음", mode, sid)
                    for k in ("peak_amp","trough_depth","oscillation","period_score"):
                        row_data[f"{mode}_{k}"] = np.nan
                    continue

                # 그룹 profiles 수집
                (l_group if mode == "L" else s_group).setdefault(group, []).append(profile)

                # 통계
                st = compute_stats(profile)
                for k, v in st.items():
                    row_data[f"{mode}_{k}"] = v

            per_sample_rows.append(row_data)

        if not per_sample_rows:
            log.warning("[%s] 데이터 없음", mtype)
            continue

        # per_sample TSV
        df_per = pd.DataFrame(per_sample_rows)
        per_tsv = str(outdir / f"{mtype}.per_sample.tsv")
        df_per.to_csv(per_tsv, sep="\t", index=False, float_format="%.5f")
        log.info("샘플별 통계: %s", per_tsv)

        # group_summary TSV
        stat_cols = [c for c in df_per.columns
                     if any(c.startswith(f"{m}_") for m in ("L","S"))]
        summary = df_per.groupby("group")[stat_cols].agg(["mean","std"]).round(4)
        sum_tsv = str(outdir / f"{mtype}.group_summary.tsv")
        summary.to_csv(sum_tsv, sep="\t")
        log.info("그룹 요약: %s", sum_tsv)

        print(f"\n=== {mtype} 샘플별 L/S WPS 통계 ===")
        print(df_per[["sample_id","group"] + stat_cols].to_string(index=False))

        # L plot
        if l_group:
            plot_mode_compare(
                group_profiles = l_group,
                group_labels   = group_labels,
                mode           = "L",
                marker_type    = mtype,
                flank          = eff_flank,
                metric         = args.metric_l,
                output_path    = str(outdir / f"{mtype}.L.wps_compare.png"),
                smooth_window  = args.smooth_window,
            )

        # S plot
        if s_group:
            plot_mode_compare(
                group_profiles = s_group,
                group_labels   = group_labels,
                mode           = "S",
                marker_type    = mtype,
                flank          = eff_flank,
                metric         = args.metric_s,
                output_path    = str(outdir / f"{mtype}.S.wps_compare.png"),
                smooth_window  = args.smooth_window,
            )

        # L+S 합친 subplot
        if not args.no_combined and (l_group or s_group):
            plot_ls_combined(
                l_group      = l_group,
                s_group      = s_group,
                group_labels = group_labels,
                marker_type  = mtype,
                flank        = eff_flank,
                output_path  = str(outdir / f"{mtype}.LS.wps_compare.png"),
                smooth_window = args.smooth_window,
            )

    log.info("완료 → %s", outdir)


if __name__ == "__main__":
    main()