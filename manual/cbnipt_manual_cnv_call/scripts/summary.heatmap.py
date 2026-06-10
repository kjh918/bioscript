"""
Single-Sample Chromosome Abnormality Detector
=============================================
샘플 그룹 비교 없이, 단일 샘플만으로 염색체 이상을 탐지합니다.

핵심 전략:
1. Robust LOO  - MAD 기반 정규화로 이상 염색체가 기준값에 영향 주지 않음
2. Bin-level 분포 검정 - 각 염색체의 bin 분포 vs 전체 분포 (Wilcoxon/KS test)
3. BAF 패턴 분석 - hetero/homo/imbalance rate로 3n/LOH/정상 구분
4. TER 이상 - Trans Error Rate 상승 = 구조이상 신호
5. 종합 판정 - 위 신호들을 합산해 최종 call

출력:
  - 01_chromosome_overview.png     : 염색체별 4개 지표 라인 플롯
  - 02_bin_distribution.png        : 이상 염색체 bin-level 분포
  - 03_anomaly_summary_heatmap.png : 지표별 이상 스코어 히트맵
  - 04_final_call.png              : 최종 판정 결과
  - console 요약 출력
"""

import sys
import os
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from scipy import stats
from pathlib import Path


# ─────────────────────────────────────────────────────────────
# 0. 설정값 (필요시 수정)
# ─────────────────────────────────────────────────────────────

CFG = dict(
    # QC 필터
    min_depth        = 1,
    min_coverage     = 0.5,
    min_sites        = 3,        # BAF 계산용 최소 SNP site 수

    # 이상 탐지 임계값
    robust_z_thresh  = 2.5,      # Robust Z-score (MAD 기반) 절댓값 기준
    pval_thresh      = 0.05,     # Wilcoxon/KS p-value 기준
    baf_dev_thresh   = 0.08,     # |median BAF - 0.5| 기준 (trisomy ≈ 0.17)
    homo_rate_thresh = 0.6,      # homo_like_rate 기준 (LOH 신호)
    ter_z_thresh     = 2.0,      # TER Robust Z 기준

    # 최종 판정 가중치 (합이 1이 되도록)
    w_copy     = 0.40,           # copy number (Log2 / OE_Ratio)
    w_baf      = 0.30,           # BAF 패턴
    w_ter      = 0.15,           # TER 이상
    w_pval     = 0.15,           # 분포 검정

    # 판정 임계값
    call_thresh_high = 0.50,     # 이 이상이면 Abnormal
    call_thresh_low  = 0.25,     # 이 이상이면 Suspicious
)


# ─────────────────────────────────────────────────────────────
# 1. 유틸리티
# ─────────────────────────────────────────────────────────────

def chrom_key(c):
    s = str(c).replace("chr", "")
    if s == "X": return 23
    if s == "Y": return 24
    try: return int(s)
    except: return 99


def sort_chroms(chroms):
    return sorted(chroms, key=chrom_key)


def mad_zscore(series):
    """MAD 기반 Robust Z-score. 이상값에 강건함."""
    med = series.median()
    mad = (series - med).abs().median()
    if mad < 1e-9:
        return pd.Series(np.zeros(len(series)), index=series.index)
    return (series - med) / (1.4826 * mad)  # 1.4826 = 정규분포 일관성 보정


def safe_log2fc(val, ref):
    #return np.log2((val + 1e-9) / (ref + 1e-9))
    return float((val + 1e-9) / (ref + 1e-9))


# ─────────────────────────────────────────────────────────────
# 2. QC 필터
# ─────────────────────────────────────────────────────────────

def apply_qc_filter(df):
    sex_mask = df["chrom"].isin(["chrX", "chrY"])
    auto = df[~sex_mask].copy()
    sex  = df[sex_mask].copy()

    auto_f = auto[
        (auto["raw_count"]     >= CFG["min_depth"]) &
        (auto["breadth_ratio"] >= CFG["min_coverage"])
    ].copy()

    out = pd.concat([auto_f, sex], ignore_index=True)
    out["_k"] = out["chrom"].apply(chrom_key)
    out = out.sort_values(["_k", "start"]).drop(columns=["_k"]).reset_index(drop=True)
    return out


# ─────────────────────────────────────────────────────────────
# 3. 염색체별 요약 통계
# ─────────────────────────────────────────────────────────────

def compute_chrom_summary(df):
    """
    각 염색체에 대해 핵심 지표의 중앙값/평균을 계산합니다.
    """
    agg = {}
    base_cols = [
        "raw_count", "breadth_ratio", "hetero_like_rate",
        "homo_like_rate", "imbalance_rate", "raw_bin_TER",
        "raw_bin_CER", "bin_BAF", "MAD_BAF",
    ]
    extra_cols = ["log2_chrom_norm", "Informative_OE_Ratio"]

    all_cols = base_cols + [c for c in extra_cols if c in df.columns]

    for chrom, grp in df.groupby("chrom"):
        row = {"chrom": chrom, "n_bins": len(grp)}
        for col in all_cols:
            if col in grp.columns:
                row[f"{col}_median"] = grp[col].median()
                row[f"{col}_mean"]   = grp[col].mean()
                row[f"{col}_std"]    = grp[col].std()
        agg[chrom] = row

    summary = pd.DataFrame(list(agg.values()))
    summary["_k"] = summary["chrom"].apply(chrom_key)
    summary = summary.sort_values("_k").drop(columns=["_k"]).reset_index(drop=True)
    return summary


# ─────────────────────────────────────────────────────────────
# 4. Robust LOO 스코어 (Copy Number 신호)
# ─────────────────────────────────────────────────────────────

def compute_robust_loo(summary, value_col):
    """
    각 염색체를 제외한 나머지의 MAD 기반 중앙값으로 Log2FC 계산.
    이상 염색체가 기준값에 영향을 주지 않음.

    반환: {chrom: (log2fc, robust_z)}
    """
    if value_col not in summary.columns:
        return {}

    chroms = summary["chrom"].tolist()
    vals   = summary[value_col].values
    result = {}

    for i, c in enumerate(chroms):
        others = np.array([v for j, v in enumerate(vals) if j != i and not np.isnan(v)])
        if len(others) < 3:
            result[c] = (np.nan, np.nan)
            continue

        ref_med = np.median(others)
        ref_mad = np.median(np.abs(others - ref_med))
        ref_mad = max(ref_mad, 1e-9)

        my_val  = vals[i]
        log2fc  = safe_log2fc(my_val, ref_med)
        robust_z = (my_val - ref_med) / (1.4826 * ref_mad)

        result[c] = (log2fc, robust_z)

    return result


# ─────────────────────────────────────────────────────────────
# 5. Bin-level 분포 검정 (각 염색체 vs 전체)
# ─────────────────────────────────────────────────────────────

def compute_bin_distribution_test(df, metric="raw_count"):
    """
    각 염색체의 bin 값 분포가 전체(자신 제외) 분포와 다른지 검정.
    Wilcoxon rank-sum test (Mann-Whitney U).
    반환: {chrom: (statistic, pvalue, direction)}
    """
    if metric not in df.columns:
        return {}

    global_vals = df[metric].dropna().values
    result = {}

    for chrom, grp in df.groupby("chrom"):
        chrom_vals = grp[metric].dropna().values
        others     = df[df["chrom"] != chrom][metric].dropna().values

        if len(chrom_vals) < 3 or len(others) < 3:
            result[chrom] = (np.nan, np.nan, 0)
            continue

        stat, pval = stats.mannwhitneyu(chrom_vals, others, alternative="two-sided")
        direction  = 1 if chrom_vals.mean() > others.mean() else -1
        result[chrom] = (stat, pval, direction)

    return result


# ─────────────────────────────────────────────────────────────
# 6. BAF 패턴 분석
# ─────────────────────────────────────────────────────────────

def compute_baf_signal(summary):
    """
    BAF / hetero / homo / imbalance rate 기반 이상 판정.

    정상 2n:  BAF ≈ 0.5,  hetero_high,  homo_low
    Trisomy:  BAF ≈ 0.33 or 0.67 (ABB/AAB),  imbalance_high
    LOH:      BAF → 0 or 1,  homo_high
    반환: {chrom: baf_score}  (0=정상, 1=최대이상)
    """
    result = {}
    for _, row in summary.iterrows():
        chrom = row["chrom"]
        score = 0.0
        n     = 0

        # BAF deviation from 0.5
        baf = row.get("bin_BAF_median", np.nan)
        if not np.isnan(baf):
            dev = abs(baf - 0.5)
            # 0.08 이상이면 이상 신호 시작, 0.17(trisomy 이론값)에서 1.0
            baf_s = min(dev / 0.17, 1.0)
            score += baf_s
            n += 1

        # homo_like_rate 상승 → LOH/monosomy
        homo = row.get("homo_like_rate_mean", np.nan)
        if not np.isnan(homo):
            homo_s = min(max(homo - CFG["homo_rate_thresh"], 0) / (1 - CFG["homo_rate_thresh"]), 1.0)
            score += homo_s
            n += 1

        # imbalance_rate 상승 → trisomy/CNV
        imbal = row.get("imbalance_rate_mean", np.nan)
        if not np.isnan(imbal):
            score += min(imbal * 2, 1.0)
            n += 1

        result[chrom] = score / n if n > 0 else 0.0

    return result


# ─────────────────────────────────────────────────────────────
# 7. 최종 종합 스코어 및 판정
# ─────────────────────────────────────────────────────────────

def make_final_call(chroms, loo_cn, bin_pvals, baf_signal, ter_loo):
    """
    4가지 신호를 가중 합산하여 최종 판정.

    반환: DataFrame (chrom, copy_score, baf_score, ter_score, pval_score, final_score, call, call_detail)
    """
    records = []
    for c in chroms:
        # 1. Copy number (Log2FC or OE Ratio 기반 LOO)
        if c in loo_cn:
            log2fc, rz = loo_cn[c]
            # |Robust Z| → 0~1 스코어 (방향 보존)
            copy_score = float(np.clip(rz / (CFG["robust_z_thresh"] * 2), -1, 1)) if not np.isnan(rz) else 0.0
        else:
            log2fc, rz, copy_score = np.nan, np.nan, 0.0

        # 2. BAF
        baf_score = baf_signal.get(c, 0.0)

        # 3. TER
        if c in ter_loo:
            _, ter_rz = ter_loo[c]
            ter_score = float(np.clip(ter_rz / (CFG["ter_z_thresh"] * 2), -1, 1)) if not np.isnan(ter_rz) else 0.0
        else:
            ter_score = 0.0

        # 4. 분포 검정 p-value → 스코어
        if c in bin_pvals:
            _, pval, direction = bin_pvals[c]
            if not np.isnan(pval) and pval < CFG["pval_thresh"]:
                # -log10(p) 를 0~1로 압축
                pval_score = float(np.clip(-np.log10(pval + 1e-10) / 5, 0, 1)) * direction
            else:
                pval_score = 0.0
        else:
            pval_score = 0.0

        # 가중 합산
        final = (
            CFG["w_copy"] * copy_score +
            CFG["w_baf"]  * baf_score  +
            CFG["w_ter"]  * ter_score  +
            CFG["w_pval"] * pval_score
        )

        # 판정
        abs_final = abs(final)
        if abs_final >= CFG["call_thresh_high"]:
            call = "ABNORMAL"
        elif abs_final >= CFG["call_thresh_low"]:
            call = "SUSPICIOUS"
        else:
            call = "NORMAL"

        # 세부 설명
        detail_parts = []
        if not np.isnan(log2fc):
            detail_parts.append(f"Log2FC={log2fc:+.2f}")
        if not np.isnan(rz):
            detail_parts.append(f"RobustZ={rz:+.2f}")
        detail_parts.append(f"BAF_s={baf_score:.2f}")

        records.append({
            "chrom": c,
            "log2fc": log2fc,
            "robust_z": rz,
            "copy_score": copy_score,
            "baf_score": baf_score,
            "ter_score": ter_score,
            "pval_score": pval_score,
            "final_score": final,
            "call": call,
            "detail": "  ".join(detail_parts),
        })

    return pd.DataFrame(records)


# ─────────────────────────────────────────────────────────────
# 8. 시각화
# ─────────────────────────────────────────────────────────────

CALL_COLORS = {"NORMAL": "#4CAF50", "SUSPICIOUS": "#FF9800", "ABNORMAL": "#F44336"}

def plot_chromosome_overview(df_filtered, summary, call_df, sample_id, out_dir):
    """
    Figure 1: 염색체별 4개 지표 프로필 (bin-level scatter + 염색체 중앙값 라인)
    """
    chroms = sort_chroms(summary["chrom"].tolist())
    chrom_x = {c: i for i, c in enumerate(chroms)}

    fig, axes = plt.subplots(4, 1, figsize=(20, 18), sharex=True)
    fig.suptitle(f"Chromosome Overview — {sample_id}", fontsize=14, fontweight="bold")

    panels = [
        ("log2_chrom_norm",  "Copy Number (Norm Log2)",    "steelblue",  (-2.0, 2.0),  [0.585, -1.0]),
        ("hetero_like_rate",       "Hetero-like Rate",           "mediumseagreen", (0, 1), []),
        ("homo_like_rate",         "Homo-like Rate (LOH signal)","tomato",     (0, 1),   []),
        ("raw_bin_TER",            "Trans Error Rate",           "darkorange", (0, None), []),
    ]
    # log2_chrom_norm 없으면 raw_count로 대체
    if "log2_chrom_norm" not in df_filtered.columns:
        panels[0] = ("raw_count", "Raw Count (bin-level)", "steelblue", (0, None), [])

    for ax, (col, title, color, ylim, ref_lines) in zip(axes, panels):
        if col not in df_filtered.columns:
            ax.set_title(f"{title} (data not available)", fontsize=10)
            continue

        # bin-level scatter
        for chrom, grp in df_filtered.groupby("chrom"):
            xi = chrom_x[chrom]
            jitter = np.random.uniform(-0.3, 0.3, len(grp))
            ax.scatter(xi + jitter, grp[col], alpha=0.3, s=8, color=color)

        # 염색체 중앙값 라인
        med_col = f"{col}_median"
        if med_col in summary.columns:
            xs = [chrom_x[c] for c in chroms if c in chrom_x]
            ys = [summary.loc[summary["chrom"] == c, med_col].values[0]
                  for c in chroms if c in chrom_x]
            ax.plot(xs, ys, color="black", linewidth=1.5, zorder=5, alpha=0.7)
            ax.scatter(xs, ys, color="black", s=30, zorder=6)

        # 기준선
        ax.axhline(0, color="gray", linewidth=0.8, linestyle="--", alpha=0.5)
        for ref in ref_lines:
            ax.axhline(ref, color="red", linewidth=1.0, linestyle=":", alpha=0.7,
                       label=f"Expected {ref:+.2f}")

        # 이상 염색체 배경 강조
        for _, row in call_df.iterrows():
            c = row["chrom"]
            if c not in chrom_x: continue
            xi = chrom_x[c]
            if row["call"] == "ABNORMAL":
                ax.axvspan(xi - 0.5, xi + 0.5, alpha=0.15, color="red", zorder=0)
            elif row["call"] == "SUSPICIOUS":
                ax.axvspan(xi - 0.5, xi + 0.5, alpha=0.10, color="orange", zorder=0)

        ax.set_ylabel(title, fontsize=9, fontweight="bold")
        ax.set_ylim(ylim)
        ax.grid(axis="y", linestyle="--", alpha=0.3)

    axes[-1].set_xticks(range(len(chroms)))
    axes[-1].set_xticklabels(chroms, rotation=45, ha="right", fontsize=9)

    plt.tight_layout()
    out = os.path.join(out_dir, "01_chromosome_overview.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  [+] Saved: {out}")


def plot_bin_distribution(df_filtered, call_df, sample_id, out_dir,
                          col="log2_chrom_norm"):
    """
    Figure 2: 이상 염색체의 bin-level 분포 vs 전체 분포 비교 (violinplot)
    """
    if col not in df_filtered.columns:
        col = "raw_count"
    if col not in df_filtered.columns:
        return

    abnormal_chroms = call_df[call_df["call"].isin(["ABNORMAL", "SUSPICIOUS"])]["chrom"].tolist()
    if not abnormal_chroms:
        print("  [i] No abnormal chromosomes — Figure 2 skipped.")
        return

    n = len(abnormal_chroms)
    fig, axes = plt.subplots(1, n, figsize=(max(6, 5 * n), 6), squeeze=False)
    fig.suptitle(f"Bin Distribution — Abnormal/Suspicious Chromosomes ({sample_id})",
                 fontsize=12, fontweight="bold")

    global_vals = df_filtered[col].dropna().values

    for ax, chrom in zip(axes[0], abnormal_chroms):
        chrom_vals = df_filtered[df_filtered["chrom"] == chrom][col].dropna().values
        others     = df_filtered[df_filtered["chrom"] != chrom][col].dropna().values

        # violin
        parts = ax.violinplot([others, chrom_vals], positions=[0, 1],
                              showmedians=True, showextrema=True)
        for i, (body, color) in enumerate(zip(parts["bodies"],
                                               ["#90CAF9", "#EF9A9A"])):
            body.set_facecolor(color)
            body.set_alpha(0.7)

        ax.scatter([0] * len(others),   others,      alpha=0.15, s=6, color="#1565C0")
        ax.scatter([1] * len(chrom_vals), chrom_vals, alpha=0.5,  s=8, color="#B71C1C")

        # 이론값 기준선
        ax.axhline(0.585, color="red", linestyle=":", linewidth=1, label="Trisomy (+0.585)")
        ax.axhline(-1.0,  color="blue", linestyle=":", linewidth=1, label="Monosomy (-1.0)")
        ax.axhline(0,     color="gray", linestyle="--", linewidth=0.8)

        _, pval, _ = (None, np.nan, None)
        row = call_df[call_df["chrom"] == chrom]
        if not row.empty:
            pval = row.iloc[0].get("pval_score", np.nan)

        call = call_df[call_df["chrom"] == chrom]["call"].values[0]
        color = CALL_COLORS.get(call, "black")

        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Others\n(all chroms)", f"{chrom}\n({call})"],
                           fontsize=10, fontweight="bold")
        ax.set_title(chrom, fontsize=12, fontweight="bold", color=color)
        ax.set_ylabel(col, fontsize=9)
        ax.legend(fontsize=7, loc="upper right")
        ax.grid(axis="y", linestyle="--", alpha=0.3)

    plt.tight_layout()
    out = os.path.join(out_dir, "02_bin_distribution.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  [+] Saved: {out}")


def plot_score_heatmap(call_df, sample_id, out_dir):
    """
    Figure 3: 염색체별 스코어 성분 히트맵
    """
    score_cols = ["copy_score", "baf_score", "ter_score", "pval_score", "final_score"]
    col_labels  = ["Copy\n(LOO RZ)", "BAF\nPattern", "TER\n(Error)", "Distrib\n(pval)", "FINAL\nScore"]

    chroms = sort_chroms(call_df["chrom"].tolist())
    mat = call_df.set_index("chrom").reindex(chroms)[score_cols].values.T  # (5, n_chrom)

    fig, ax = plt.subplots(figsize=(max(14, len(chroms) * 0.7), 5))
    fig.suptitle(f"Anomaly Score Components — {sample_id}",
                 fontsize=12, fontweight="bold")

    im = ax.imshow(mat, cmap="RdBu_r", aspect="auto", vmin=-1, vmax=1)

    for i in range(len(score_cols)):
        for j, c in enumerate(chroms):
            val = mat[i, j]
            text_color = "white" if abs(val) > 0.6 else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                    color=text_color, fontsize=8, fontweight="bold")

            # 최종 판정에 따라 열 테두리
            if i == len(score_cols) - 1:
                call = call_df[call_df["chrom"] == c]["call"].values[0]
                if call == "ABNORMAL":
                    rect = plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                         linewidth=3, edgecolor="red", facecolor="none")
                    ax.add_patch(rect)
                elif call == "SUSPICIOUS":
                    rect = plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                         linewidth=2, edgecolor="orange", facecolor="none")
                    ax.add_patch(rect)

    ax.set_xticks(range(len(chroms)))
    ax.set_yticks(range(len(score_cols)))
    ax.set_xticklabels(chroms, rotation=45, ha="right", fontsize=9, fontweight="bold")
    ax.set_yticklabels(col_labels, fontsize=9, fontweight="bold")

    plt.colorbar(im, ax=ax, pad=0.02, fraction=0.025).set_label(
        "Score (-1=low/normal, +1=high/abnormal)", fontsize=9)

    ax.set_xticks(np.arange(-0.5, len(chroms), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(score_cols), 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=1.5)
    ax.tick_params(which="minor", bottom=False, left=False)

    plt.tight_layout()
    out = os.path.join(out_dir, "03_anomaly_score_heatmap.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  [+] Saved: {out}")


def plot_final_call(call_df, sample_id, out_dir):
    """
    Figure 4: 최종 판정 요약 — 염색체별 final_score 바 차트 + call 컬러
    """
    chroms = sort_chroms(call_df["chrom"].tolist())
    call_df_sorted = call_df.set_index("chrom").reindex(chroms).reset_index()

    fig, ax = plt.subplots(figsize=(max(16, len(chroms) * 0.9), 6))

    colors = [CALL_COLORS.get(c, "gray") for c in call_df_sorted["call"]]
    bars = ax.bar(range(len(chroms)), call_df_sorted["final_score"].values,
                  color=colors, edgecolor="black", linewidth=0.6, width=0.75)

    ax.axhline( CFG["call_thresh_high"], color="red",    linestyle="--", linewidth=1.5,
                label=f"Abnormal threshold ({CFG['call_thresh_high']})")
    ax.axhline( CFG["call_thresh_low"],  color="orange", linestyle="--", linewidth=1.2,
                label=f"Suspicious threshold ({CFG['call_thresh_low']})")
    ax.axhline(-CFG["call_thresh_high"], color="red",    linestyle="--", linewidth=1.5)
    ax.axhline(-CFG["call_thresh_low"],  color="orange", linestyle="--", linewidth=1.2)
    ax.axhline(0, color="black", linewidth=0.8)

    ax.set_xticks(range(len(chroms)))
    ax.set_xticklabels(chroms, rotation=45, ha="right", fontsize=10, fontweight="bold")
    ax.set_ylabel("Final Anomaly Score", fontsize=11, fontweight="bold")
    ax.set_ylim(-1.1, 1.1)
    ax.set_title(f"Final Chromosome Call — {sample_id}", fontsize=13, fontweight="bold")

    # 이상 염색체에 레이블
    for i, (_, row) in enumerate(call_df_sorted.iterrows()):
        if row["call"] != "NORMAL":
            ypos = row["final_score"]
            offset = 0.05 if ypos >= 0 else -0.10
            ax.text(i, ypos + offset, row["call"], ha="center", va="bottom",
                    fontsize=8, fontweight="bold",
                    color=CALL_COLORS[row["call"]])
            # Robust Z 표시
            if not np.isnan(row["robust_z"]):
                ax.text(i, ypos + offset - 0.12, f"Z={row['robust_z']:+.1f}",
                        ha="center", fontsize=7, color="gray")

    patches = [mpatches.Patch(color=v, label=k) for k, v in CALL_COLORS.items()]
    ax.legend(handles=patches + [
        plt.Line2D([0],[0], color="red",    linestyle="--", label=f"Abnormal ≥{CFG['call_thresh_high']}"),
        plt.Line2D([0],[0], color="orange", linestyle="--", label=f"Suspicious ≥{CFG['call_thresh_low']}"),
    ], fontsize=9, loc="upper right")

    ax.grid(axis="y", linestyle="--", alpha=0.35)
    plt.tight_layout()

    out = os.path.join(out_dir, "04_final_call.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  [+] Saved: {out}")


# ─────────────────────────────────────────────────────────────
# 9. 메인
# ─────────────────────────────────────────────────────────────

def run(tsv_path: str, out_dir: str = None):
    tsv_path = Path(tsv_path)
    sample_id = tsv_path.stem.replace(".normalized", "")

    if out_dir is None:
        out_dir = str(tsv_path.parent)
    os.makedirs(out_dir, exist_ok=True)

    print("=" * 65)
    print(f"[★] Single-Sample CNV Analysis: {sample_id}")
    print("=" * 65)

    # ── 데이터 로드 & QC ──
    df = pd.read_csv(tsv_path, sep="\t")
    df_f = apply_qc_filter(df)
    print(f"  Bins after QC: {len(df_f)}  (from {len(df)})")

    auto_df = df_f[~df_f["chrom"].isin(["chrX", "chrY"])].copy()
    chroms  = sort_chroms(auto_df["chrom"].unique().tolist())

    # ── 요약 통계 ──
    summary = compute_chrom_summary(auto_df)
    #print(summary[['chrom','log2_chrom_norm_median']])
    #exit()
    # ── 신호 계산 ──
    # Copy number: log2_chrom_norm 우선, 없으면 raw_count 사용
    cn_col = "log2_chrom_norm_median" if "log2_chrom_norm_median" in summary.columns \
             else "raw_count_median"
    loo_cn  = compute_robust_loo(summary, cn_col)
    #loo_cn  = summary / cn_col
    loo_ter = compute_robust_loo(summary, "raw_bin_TER_median")

    bin_pvals = compute_bin_distribution_test(
        auto_df,
        metric="log2_chrom_norm" if "log2_chrom_norm" in auto_df.columns else "raw_count"
    )
    baf_sig = compute_baf_signal(summary)

    # ── 최종 판정 ──
    call_df = make_final_call(chroms, loo_cn, bin_pvals, baf_sig, loo_ter)
    #print(call_df)
    #exit()
    # ── 콘솔 요약 ──
    print(f"\n{'Chrom':<8} {'Log2FC':>8} {'RobustZ':>9} {'BAF_s':>7} {'Final':>8}  Call")
    print("-" * 58)
    for _, row in call_df.iterrows():
        flag = "⚠ " if row["call"] == "ABNORMAL" else ("△ " if row["call"] == "SUSPICIOUS" else "  ")
        log2 = f"{row['log2fc']:+.3f}" if not np.isnan(row["log2fc"]) else "   N/A"
        rz   = f"{row['robust_z']:+.2f}" if not np.isnan(row["robust_z"]) else "  N/A"
        print(f"{row['chrom']:<8} {log2:>8} {rz:>9} {row['baf_score']:>7.3f} "
              f"{row['final_score']:>8.3f}  {flag}{row['call']}")

    abnormal = call_df[call_df["call"] == "ABNORMAL"]["chrom"].tolist()
    susp     = call_df[call_df["call"] == "SUSPICIOUS"]["chrom"].tolist()
    print(f"\n  ABNORMAL   : {abnormal if abnormal else 'none'}")
    print(f"  SUSPICIOUS : {susp if susp else 'none'}")

    # ── 시각화 ──
    print("\n[★] Generating figures ...")
    plot_chromosome_overview(auto_df, summary, call_df, sample_id, out_dir)
    plot_bin_distribution(auto_df, call_df, sample_id, out_dir)
    plot_score_heatmap(call_df, sample_id, out_dir)
    plot_final_call(call_df, sample_id, out_dir)

    print(f"\n[Done] Results in: {out_dir}")
    return call_df


# ─────────────────────────────────────────────────────────────
# 10. CLI 실행
# ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # ── 직접 경로 지정 방식 ──────────────────────────────────
    # 아래 TSV_PATH 를 분석할 파일 경로로 변경하고 실행
    TSV_PATH = "/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Results/cbNIPT_24_04_13/data/cbNIPT_24_04_13.normalized.tsv"
    OUT_DIR  = None   # None이면 TSV 파일과 같은 폴더에 저장

    # ── CLI 인자 방식 (선택) ─────────────────────────────────
    if len(sys.argv) > 1:
        parser = argparse.ArgumentParser(description="Single-sample CNV/Aneuploidy detector")
        parser.add_argument("tsv",     help="Path to .normalized.tsv")
        parser.add_argument("--out",   default=None, help="Output directory")
        args = parser.parse_args()
        TSV_PATH = args.tsv
        OUT_DIR  = args.out

    run(TSV_PATH, OUT_DIR)