"""
wps_diagnostic.py
=================
bins_wps_raw.parquet 를 읽어서 WPS 값을 직접 확인하고
normalization 을 단계별로 적용하면서 결과를 PNG 로 저장합니다.

사용법
------
python wps_diagnostic.py \
    --raw  /path/to/bins_wps_raw.parquet \
    --out  ./wps_diag \
    --chrom chr9

선택 옵션
---------
--norm  /path/to/bins_wps_norm.parquet   파이프라인 결과와 비교
--chrom chr9                              특정 염색체만 (기본: 전체 4트랙 + chr9 확대)
"""

from __future__ import annotations

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd


CHROM_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
WPS_COLS = [
    "short_wps_L", "long_wps_L",
    "short_wps_S", "long_wps_S",
]
COLORS = {
    "short_wps_L": "#dc5050",
    "long_wps_L":  "#3050c8",
    "short_wps_S": "#c06000",
    "long_wps_S":  "#2090a0",
}


# ─────────────────────────────────────────────────────────────────────
# normalization
# ─────────────────────────────────────────────────────────────────────
def zscore_chrom(arr: np.ndarray, quality_mask: np.ndarray) -> np.ndarray:
    ref = arr[quality_mask & np.isfinite(arr)]
    if len(ref) < 5:
        return np.full(len(arr), np.nan, dtype=np.float32)
    med   = float(np.median(ref))
    mad   = float(np.median(np.abs(ref - med)))
    scale = mad * 1.4826 if mad > 1e-6 else float(ref.std() or 1.0)
    out   = np.full(len(arr), np.nan, dtype=np.float32)
    fin   = np.isfinite(arr)
    out[fin] = ((arr[fin] - med) / scale).astype(np.float32)
    return out


def normalize_df(df: pd.DataFrame) -> pd.DataFrame:
    """bins_wps_raw DataFrame → normalized WPS 컬럼 추가."""
    df = df.copy()

    # quality mask
    is_sex = df["chrom"].isin({"chrX", "chrY"}).values
    mp = df["mappability"].fillna(np.nan).values if "mappability" in df.columns \
         else np.ones(len(df))
    gc = df["gc"].values.astype(float) if "gc" in df.columns \
         else np.ones(len(df))

    if "is_filtered" in df.columns:
        bad = df["is_filtered"].astype(bool).values & ~is_sex
    elif "is_low_mappability" in df.columns:
        bad = df["is_low_mappability"].astype(bool).values & ~is_sex
    else:
        bad = (~np.isnan(mp) & (mp < 0.75) & ~is_sex)
    bad |= (~np.isfinite(gc) | (gc <= 0))
    quality_mask = ~bad & ~is_sex

    for col in [c for c in WPS_COLS if c in df.columns]:
        raw = df[col].values.astype(np.float32)
        raw[bad] = np.nan
        norm = np.full(len(df), np.nan, dtype=np.float32)

        for chrom, cdf in df.groupby("chrom"):
            idx     = cdf.index
            wps_c   = raw[idx]
            qmask_c = quality_mask[idx]
            norm[idx] = zscore_chrom(wps_c, qmask_c)

        # winsorization
        fin = norm[np.isfinite(norm)]
        if len(fin) > 100:
            lo = float(np.quantile(fin, 0.005))
            hi = float(np.quantile(fin, 0.995))
            norm = np.clip(norm, lo, hi)

        df[f"{col}_norm"] = norm

    return df


# ─────────────────────────────────────────────────────────────────────
# 플롯 함수
# ─────────────────────────────────────────────────────────────────────
def plot_raw_stats(df: pd.DataFrame, out_dir: str) -> None:
    """WPS raw 값 분포 및 기초 통계."""
    cols = [c for c in WPS_COLS if c in df.columns]
    if not cols:
        print("WPS 컬럼 없음")
        return

    fig, axes = plt.subplots(len(cols), 2, figsize=(14, 3 * len(cols)),
                             facecolor="white")
    if len(cols) == 1:
        axes = [axes]

    for ax_row, col in zip(axes, cols):
        vals = df[col].dropna().values
        color = COLORS.get(col, "gray")

        # 분포
        ax_row[0].hist(vals, bins=100, color=color, alpha=0.7)
        ax_row[0].axvline(float(np.median(vals)), color="black",
                          ls="--", lw=1.2, label=f"median={np.median(vals):.1f}")
        ax_row[0].set_title(f"{col} 분포")
        ax_row[0].set_xlabel("raw WPS")
        ax_row[0].legend(fontsize=9)

        # 염색체별 median
        chrom_med = df.groupby("chrom")[col].median().reindex(CHROM_ORDER).dropna()
        ax_row[1].bar(range(len(chrom_med)), chrom_med.values, color=color, alpha=0.7)
        ax_row[1].set_xticks(range(len(chrom_med)))
        ax_row[1].set_xticklabels(
            [c.replace("chr","") for c in chrom_med.index], fontsize=8)
        ax_row[1].set_title(f"{col} 염색체별 median")
        ax_row[1].axhline(0, color="black", lw=0.7, ls="--")

    fig.suptitle("WPS Raw 기초 통계", fontsize=13)
    plt.tight_layout()
    out = os.path.join(out_dir, "01_raw_stats.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"저장: {out}")


def plot_genome_track(df: pd.DataFrame, col_suffix: str,
                      title: str, out_path: str) -> None:
    """전장 유전체 WPS track (4트랙)."""
    cols = [c for c in [f"short_wps_L{col_suffix}", f"long_wps_L{col_suffix}",
                         f"short_wps_S{col_suffix}", f"long_wps_S{col_suffix}"]
            if c in df.columns]
    if not cols:
        print(f"컬럼 없음: {col_suffix}")
        return

    df2 = df[df["chrom"].isin(CHROM_ORDER)].copy()
    df2["chrom"] = pd.Categorical(df2["chrom"], categories=CHROM_ORDER, ordered=True)
    df2 = df2.sort_values(["chrom", "start"]).reset_index(drop=True)

    # x 좌표
    offsets: dict[str, int] = {}
    cursor = 0
    for chrom in CHROM_ORDER:
        if chrom not in df2["chrom"].values: continue
        offsets[chrom] = cursor
        cursor += int(df2[df2["chrom"] == chrom]["end"].max())
    x = np.array([offsets.get(r.chrom, 0) + (r.start + r.end) // 2
                  for r in df2.itertuples()])

    # y 공유 범위
    all_vals = np.concatenate([df2[c].dropna().values for c in cols])
    if len(all_vals):
        lo = float(np.quantile(all_vals, 0.01))
        hi = float(np.quantile(all_vals, 0.99))
        pad = (hi - lo) * 0.15
        bound = max(abs(lo - pad), abs(hi + pad))
        y_min, y_max = -bound, bound
    else:
        y_min, y_max = -1, 1

    labels = {
        "short_wps_L": "Short L-WPS", "short_wps_L_norm": "Short L-WPS (z)",
        "long_wps_L":  "Long  L-WPS", "long_wps_L_norm":  "Long  L-WPS (z)",
        "short_wps_S": "Short S-WPS", "short_wps_S_norm": "Short S-WPS (z)",
        "long_wps_S":  "Long  S-WPS", "long_wps_S_norm":  "Long  S-WPS (z)",
    }

    fig, axes = plt.subplots(len(cols), 1, figsize=(20, 3 * len(cols)),
                             sharex=True, facecolor="white")
    if len(cols) == 1: axes = [axes]

    for ax, col in zip(axes, cols):
        base = col.replace("_norm", "")
        color = COLORS.get(base, "gray")
        y = df2[col].values.astype(float)
        ax.fill_between(x, 0, y, where=y > 0, color=color, alpha=0.7)
        ax.fill_between(x, 0, y, where=y < 0, color=color, alpha=0.4)
        ax.axhline(0, color="black", lw=0.5, ls="--")
        ax.set_ylim(y_min, y_max)
        ax.set_ylabel(labels.get(col, col), fontsize=9)
        ax.set_facecolor("white")
        for sp in ax.spines.values(): sp.set_linewidth(0.4)

        for chrom, off in offsets.items():
            ax.axvline(off, color="lightgray", lw=0.4, ls=":")
        if ax is axes[-1]:
            for chrom, off in offsets.items():
                sub = df2[df2["chrom"] == chrom]
                if not sub.empty:
                    mid = off + int(sub["end"].max()) // 2
                    ax.text(mid, y_min * 1.15, chrom.replace("chr", ""),
                            ha="center", va="top", fontsize=7)

    ax.set_xlabel("Genomic position", fontsize=10)
    fig.suptitle(title, fontsize=13)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"저장: {out_path}")


def plot_chrom_zoom(df: pd.DataFrame, chrom: str, out_dir: str) -> None:
    """특정 염색체 확대 (raw vs norm 비교)."""
    cdf = df[df["chrom"] == chrom].sort_values("start").reset_index(drop=True)
    if cdf.empty:
        print(f"{chrom} 데이터 없음")
        return

    x = ((cdf["start"] + cdf["end"]) / 2).values

    raw_cols  = [c for c in WPS_COLS if c in cdf.columns]
    norm_cols = [f"{c}_norm" for c in raw_cols if f"{c}_norm" in cdf.columns]
    pairs     = [(c, f"{c}_norm") for c in raw_cols if f"{c}_norm" in cdf.columns]

    if not pairs:
        print(f"norm 컬럼 없음 — raw만 출력")
        pairs = [(c, None) for c in raw_cols]

    n_rows = len(pairs)
    fig, axes = plt.subplots(n_rows, 2 if pairs[0][1] else 1,
                             figsize=(18, 3 * n_rows),
                             sharex=True, facecolor="white")
    if n_rows == 1: axes = [axes] if pairs[0][1] is None else [[axes[0], axes[1]]]

    for row_axes, (raw_col, norm_col) in zip(axes, pairs):
        color = COLORS.get(raw_col, "gray")

        if isinstance(row_axes, (list, np.ndarray)):
            ax_raw, ax_norm = row_axes[0], row_axes[1]
        else:
            ax_raw, ax_norm = row_axes, None

        # raw
        y_raw = cdf[raw_col].values.astype(float)
        ax_raw.fill_between(x, 0, y_raw, where=y_raw > 0, color=color, alpha=0.7)
        ax_raw.fill_between(x, 0, y_raw, where=y_raw < 0, color=color, alpha=0.4)
        ax_raw.axhline(0, color="black", lw=0.5, ls="--")
        ax_raw.set_ylabel(f"{raw_col}\n(raw)", fontsize=9)
        ax_raw.set_facecolor("white")
        med_v = float(np.nanmedian(y_raw))
        ax_raw.set_title(f"raw  median={med_v:.2f}", fontsize=9)

        # norm
        if ax_norm and norm_col:
            y_norm = cdf[norm_col].values.astype(float)
            ax_norm.fill_between(x, 0, y_norm, where=y_norm > 0, color=color, alpha=0.7)
            ax_norm.fill_between(x, 0, y_norm, where=y_norm < 0, color=color, alpha=0.4)
            ax_norm.axhline(0, color="black", lw=0.5, ls="--")
            ax_norm.set_ylabel(f"{raw_col}\n(z-score)", fontsize=9)
            ax_norm.set_facecolor("white")
            ax_norm.set_title(f"z-score  range=[{np.nanmin(y_norm):.2f}, {np.nanmax(y_norm):.2f}]", fontsize=9)

    fig.suptitle(f"{chrom} WPS: Raw vs Normalized", fontsize=13)
    plt.tight_layout()
    out = os.path.join(out_dir, f"03_{chrom}_zoom.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"저장: {out}")


def print_summary(df: pd.DataFrame) -> None:
    """터미널 요약 출력."""
    print("\n" + "=" * 60)
    print("WPS Raw 값 요약")
    print("=" * 60)
    for col in [c for c in WPS_COLS if c in df.columns]:
        vals = df[col].dropna().values
        norm_col = f"{col}_norm"
        norm_vals = df[norm_col].dropna().values if norm_col in df.columns else None

        print(f"\n[{col}]")
        print(f"  raw : n={len(vals):,}  median={np.median(vals):.3f}"
              f"  mean={np.mean(vals):.3f}"
              f"  std={np.std(vals):.3f}"
              f"  range=[{np.min(vals):.1f}, {np.max(vals):.1f}]")
        if norm_vals is not None:
            print(f"  norm: n={len(norm_vals):,}  median={np.median(norm_vals):.3f}"
                  f"  std={np.std(norm_vals):.3f}"
                  f"  range=[{np.min(norm_vals):.2f}, {np.max(norm_vals):.2f}]")

    if "mappability" in df.columns:
        mp = df["mappability"].dropna()
        print(f"\n[mappability]  median={mp.median():.3f}"
              f"  < 0.5: {(mp<0.5).sum()}  < 0.75: {(mp<0.75).sum()}")
    if "gc" in df.columns:
        gc = df["gc"].dropna()
        print(f"[gc]           median={gc.median():.3f}"
              f"  <=0: {(gc<=0).sum()}")
    print("=" * 60 + "\n")


# ─────────────────────────────────────────────────────────────────────
# main
# ─────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="WPS raw parquet 진단 및 normalization 시각화"
    )
    parser.add_argument("--raw",   required=True,
                        help="bins_wps_raw.parquet 경로")
    parser.add_argument("--norm",  default=None,
                        help="bins_wps_norm.parquet 경로 (파이프라인 결과와 비교용)")
    parser.add_argument("--out",   default="./wps_diag",
                        help="출력 디렉토리")
    parser.add_argument("--chrom", default="chr9",
                        help="확대 볼 염색체 (기본: chr9)")
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    print(f"로드 중: {args.raw}")
    df = pd.read_parquet(args.raw)
    print(f"  → {len(df):,} bins, 컬럼: {list(df.columns)}")

    # 1. raw 통계
    print_summary(df)
    plot_raw_stats(df, args.out)

    # 2. 직접 normalization 적용
    print("z-score normalization 적용 중...")
    df_norm = normalize_df(df)

    # 3. genome-wide track (raw)
    plot_genome_track(df_norm, "",
                      "WPS Genome Track — Raw",
                      os.path.join(args.out, "02_genome_raw.png"))

    # 4. genome-wide track (norm)
    plot_genome_track(df_norm, "_norm",
                      "WPS Genome Track — Z-score",
                      os.path.join(args.out, "02_genome_norm.png"))

    # 5. 파이프라인 norm 결과 비교 (있을 때)
    if args.norm and os.path.exists(args.norm):
        print(f"파이프라인 norm 로드: {args.norm}")
        df_pipe = pd.read_parquet(args.norm)
        plot_genome_track(df_pipe, "_norm",
                          "WPS Genome Track — Pipeline Normalized",
                          os.path.join(args.out, "02_genome_pipeline_norm.png"))

    # 6. 특정 염색체 확대
    plot_chrom_zoom(df_norm, args.chrom, args.out)

    print(f"\n완료. 결과: {args.out}/")
    print("  01_raw_stats.png        : raw WPS 분포 + 염색체별 median")
    print("  02_genome_raw.png       : 전장 유전체 raw WPS")
    print("  02_genome_norm.png      : 전장 유전체 z-score WPS")
    print(f"  03_{args.chrom}_zoom.png : {args.chrom} raw vs z-score 비교")


if __name__ == "__main__":
    main()