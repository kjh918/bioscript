from operator import sub

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
import numpy as np
import pandas as pd
import os
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

from utils import log, sort_chroms
from rules import CALL_COLORS, CFG
from clinical_diagnosis import _resolve_clinical_rule


# ═════════════════════════════════════════════════════════════
# 공통 스타일 스킴 (모든 plot 함수가 여기서만 색을 가져온다)
# ═════════════════════════════════════════════════════════════
STYLE = {
    # scatter / violin 기본색 (autosome vs sex-chrom vs abnormal)
    "autosome_scatter": "#1565C0",
    "autosome_violin":  "#90CAF9",
    "sex_scatter":      "#6A1B9A",
    "sex_violin":       "#CE93D8",
    "abn_scatter":      "#B71C1C",
    "abn_violin":       "#EF9A9A",

    # 배경 음영 (axvspan / rectangle 공용)
    "sex_shade":        "purple",
    "sex_shade_alpha":  0.10,
    "abn_shade":        "red",
    "abn_shade_alpha":  0.18,
    "sus_shade":        "orange",
    "sus_shade_alpha":  0.12,

    # 성염색체 전용 call 색 (final_call bar 등에서 상염색체와 구분)
    "sex_call_colors": {"NORMAL": "#AB47BC", "SUSPICIOUS": "#FF9800", "ABNORMAL": "#F44336"},
}


def _is_sex_chrom(chrom):
    return chrom in ("chrX", "chrY")


# ═════════════════════════════════════════════════════════════
# 임상 진단(Marker DB) threshold 검증용 plot
# ═════════════════════════════════════════════════════════════
def plot_clinical_diagnosis_report(report_df, out_path, sample_id="Sample", cfg=CFG):
    """
    diagnose_clinical_markers() 의 결과(report_df)를 syndrome 단위로 시각화.

    - AMP/DEL 타입 syndrome: syndrome별 서브플롯에서 각 feature(TargetChromosome
      ~ CoreGene)의 OBSERVED_CN 을 점으로 찍고, 그 syndrome에 적용된
      pos/sus threshold(config.yaml CFG["CLINICAL_DIAGNOSIS_RULES"])를
      배경 음영 + 점선으로 표시. 점 색은 DIAGNOSIS(POSITIVE/SUSPICIOUS/NEGATIVE),
      LOW_RESOLUTION_WARNING(겹치는 bin이 너무 적은 경우)이 True 이면 'x' 마커로 표시.
    - SEX_ANEUPLOIDY 타입 syndrome(Turner/Klinefelter/Jacob/TripleX 등): 하나의
      chrX-CN vs chrY-CN 평면에 config.yaml CFG["SEX_ANEUPLOIDY"] 각 규칙의
      판정 영역을 박스로 그리고, 샘플의 실측 (X, Y) 위치를 별표로 표시.
    """
    diag_colors = {"POSITIVE": "#F44336", "SUSPICIOUS": "#FF9800", "NEGATIVE": "#4CAF50"}

    df = report_df[report_df["SYNDROME"].astype(str).str.upper() != "NIPT_SEX"].copy()

    amp_del_groups = []   # (syndrome, rule, sub_df)
    sex_x_cn, sex_y_cn = None, None

    for syn, grp in df.groupby("SYNDROME", sort=False):
        syn_key = str(syn).lower()
        grp_key = str(grp["NIPT_GROUP"].iloc[0]).lower()
        rule = _resolve_clinical_rule(syn_key, grp_key)
        rtype = rule.get("rule_type", "UNKNOWN")

        if rtype in ("AMP", "DEL"):
            amp_del_groups.append((syn, rule, grp.sort_values("FEAT_RANK")))
        elif rtype == "SEX_ANEUPLOIDY":
            xrow = grp[grp["CHROMOSOME"] == "chrX"]
            yrow = grp[grp["CHROMOSOME"] == "chrY"]
            if sex_x_cn is None and not xrow.empty and pd.notna(xrow["OBSERVED_CN"].iloc[0]):
                sex_x_cn = float(xrow["OBSERVED_CN"].iloc[0])
            if sex_y_cn is None and not yrow.empty and pd.notna(yrow["OBSERVED_CN"].iloc[0]):
                sex_y_cn = float(yrow["OBSERVED_CN"].iloc[0])

    n_amp_del = len(amp_del_groups)
    n_cols = 3 if n_amp_del > 1 else 1
    n_rows_amp = int(np.ceil(n_amp_del / n_cols)) if n_amp_del else 0
    has_sex_panel = (sex_x_cn is not None) or (sex_y_cn is not None)
    total_rows = n_rows_amp + (2 if has_sex_panel else 0)   # sex panel은 2 row 높이 사용

    if total_rows == 0:
        log("[plot_clinical_diagnosis_report] No AMP/DEL/SEX_ANEUPLOIDY syndromes to plot.")
        return

    fig = plt.figure(figsize=(6.2 * n_cols, 3.0 * n_rows_amp + (6.0 if has_sex_panel else 0)))
    fig.suptitle(f"Clinical Diagnosis Thresholds — {sample_id}", fontsize=14, fontweight="bold")
    gs = fig.add_gridspec(total_rows, n_cols, hspace=0.7, wspace=0.35)
    fig.subplots_adjust(
        top=0.95,      # 그래프 시작 위치. 0.90~0.95 사이 조절
        bottom=0.06,
        left=0.08,
        right=0.98
    )

    # ── AMP/DEL 서브플롯 ─────────────────────────────────────────
    for idx, (syn, rule, sub) in enumerate(amp_del_groups):
        r, c = divmod(idx, n_cols)
        ax = fig.add_subplot(gs[r, c])
        rtype = rule["rule_type"]

        obs_vals = sub["OBSERVED_CN"].dropna()
        x_lo = min(0.0, obs_vals.min() - 0.3) if not obs_vals.empty else 0.0
        x_hi = max(4.0, obs_vals.max() + 0.3) if not obs_vals.empty else 4.0

        if rtype == "AMP":
            pos_th, sus_th = rule["pos_min"], rule["sus_min"]
            ax.axvspan(pos_th, x_hi, color="#F44336", alpha=0.10, label="POSITIVE zone")
            ax.axvspan(sus_th, pos_th, color="#FF9800", alpha=0.10, label="SUSPICIOUS zone")
        else:
            pos_th, sus_th = rule["pos_max"], rule["sus_max"]
            ax.axvspan(x_lo, pos_th, color="#F44336", alpha=0.10, label="POSITIVE zone")
            ax.axvspan(pos_th, sus_th, color="#FF9800", alpha=0.10, label="SUSPICIOUS zone")

        ax.axvline(pos_th, color="#F44336", linestyle="--", linewidth=1.2)
        ax.axvline(sus_th, color="#FF9800", linestyle="--", linewidth=1.2)
        ax.axvline(2.0, color="gray", linestyle=":", linewidth=1.0, label="Diploid (2.0)")

        y_labels, y_pos = [], []
        
        n_rows = len(sub)

        for i, (_, row) in enumerate(sub.iterrows()):
            label = row["FEATURE_NAME"]
            y_labels.append(f"{label}\n({row['FEATURE_TYPE']})")
            y_pos.append(i)

            cn = row["OBSERVED_CN"]
            color = diag_colors.get(row["DIAGNOSIS"], "gray")
            low_res = bool(row.get("LOW_RESOLUTION_WARNING", False))

            # NA인 경우: x=2 위치에 네모로 표시
            if pd.isna(cn):
                na_x = 2
                ax.text(
                    na_x, i, "NA",           # x=2 기준 가운데
                    fontsize=6.5,
                    ha="center",
                    va="center",
                    color="dimgray",
                    zorder=10
                )
                continue

            # low_res면 세모, 아니면 원
            marker = "^" if low_res else "o"

            ax.scatter(
                [cn], [i],
                color=color,
                s=90,
                marker=marker,
                edgecolor="black",
                linewidth=1,
                zorder=5
            )

            ax.text(
                cn, i + 0.22,
                f"n={int(row['OVERLAP_BINS'])}",
                fontsize=6.5,
                ha="center",
                color="dimgray",
                zorder=10
            )

        ax.set_yticks(y_pos)
        ax.set_yticklabels(y_labels, fontsize=7)
        ax.invert_yaxis()

        if n_rows == 1:
            ax.set_ylim(0.8, -0.8)
        elif n_rows == 2:
            ax.set_ylim(1.4, -0.4)
        else:
            ax.set_ylim(n_rows - 0.4, -0.6)

        ax.set_xlim(x_lo, x_hi)
        ax.set_xlabel("Observed CN", fontsize=8, fontweight="bold")
        ax.set_title(f"{syn}  [{rtype}]", fontsize=9, fontweight="bold")
        ax.grid(axis="x", linestyle="--", alpha=0.3)

        # legend 직접 생성
        legend_handles = [
            Line2D([0], [0], marker='o', color='w', label='Pass',
                markerfacecolor='gray', markeredgecolor='black', markersize=8),
            Line2D([0], [0], marker='^', color='w', label='Low resolution',
                markerfacecolor='gray', markeredgecolor='black', markersize=8),
            Line2D([0], [0], marker='s', color='w', label='NA',
                markerfacecolor='white', markeredgecolor='white', markersize=8),
        ]

        ax.legend(handles=legend_handles, fontsize=5, loc="upper right")

    # ── 성염색체 이수성 서브플롯 (chrX CN vs chrY CN 평면) ──────────
    if has_sex_panel:
        ax = fig.add_subplot(gs[n_rows_amp:n_rows_amp + 2, :])

        call_color_map = {"ABNORMAL": "#F44336", "SUSPICIOUS": "#FF9800"}
        for key, rule in cfg["SEX_ANEUPLOIDY"].items():
            x_min = rule.get("x_min", -0.3)
            x_max = rule.get("x_max", 4.0)
            y_min = rule.get("y_min", -0.3)
            y_max = rule.get("y_max", 3.0)
            face = call_color_map.get(rule["call"], "gray")
            rect = plt.Rectangle((x_min, y_min), x_max - x_min, y_max - y_min,
                                  facecolor=face, alpha=0.15, edgecolor=face, linewidth=1.3)
            ax.add_patch(rect)
            ax.text((x_min + x_max) / 2, (y_min + y_max) / 2, rule["pattern_name"],
                    ha="center", va="center", fontsize=8, fontweight="bold", color=face)

        if sex_x_cn is not None and sex_y_cn is not None:
            ax.scatter([sex_x_cn], [sex_y_cn], color="royalblue", s=220, marker="*",
                       edgecolor="black", linewidth=1.2, zorder=10,
                       label=f"Sample (X={sex_x_cn:.2f}, Y={sex_y_cn:.2f})")
            ax.legend(fontsize=8, loc="upper right")

        ax.axvline(2.0, color="gray", linestyle=":", linewidth=0.8)
        ax.axhline(1.0, color="gray", linestyle=":", linewidth=0.8)
        ax.set_xlim(-0.3, 4.0)
        ax.set_ylim(-0.3, 3.0)
        ax.set_xlabel("chrX Copy Number", fontsize=9, fontweight="bold")
        ax.set_ylabel("chrY Copy Number", fontsize=9, fontweight="bold")
        ax.set_title("Sex Chromosome Aneuploidy Thresholds",
                     fontsize=10, fontweight="bold")
        ax.grid(linestyle="--", alpha=0.3)

    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    log(f"Clinical diagnosis threshold plot saved to: {out_path}")


def plot_gc_correction(gc_stats, out_path):
    if gc_stats is None:
        return
    x, y, fit = gc_stats
    plt.figure(figsize=(8, 5))
    plt.scatter(x, y, s=1, color='gray', alpha=0.2, label="Raw")
    idx = x.argsort()
    plt.plot(x[idx], fit[idx], color='red', lw=2, label="Lowess Fit")
    plt.title("GC Bias Correction (LOWESS)")
    plt.legend()
    plt.savefig(out_path)
    plt.close()


def _draw_bin_distribution_zoom(fig, gs, row_idx, df_f, call_df, abn_chroms, cn_col, use_cn_scale):
    auto_vals = df_f[~df_f["chrom"].isin(["chrX", "chrY"])][cn_col].dropna().values

    zoom_axes = []
    for j, chrom in enumerate(abn_chroms):
        ax = fig.add_subplot(gs[row_idx, j])
        zoom_axes.append(ax)
        chrom_vals = df_f[df_f["chrom"] == chrom][cn_col].dropna().values
        is_sex = _is_sex_chrom(chrom)

        parts = ax.violinplot([auto_vals, chrom_vals], positions=[0, 1], showmedians=True, showextrema=True)
        colors = [STYLE["autosome_violin"], STYLE["sex_violin"] if is_sex else STYLE["abn_violin"]]
        for body, vc in zip(parts["bodies"], colors):
            body.set_facecolor(vc)
            body.set_alpha(0.7)

        ax.scatter([0] * len(auto_vals), auto_vals, alpha=0.1, s=5, color=STYLE["autosome_scatter"])
        ax.scatter([1] * len(chrom_vals), chrom_vals, alpha=0.4, s=8,
                    color=STYLE["sex_scatter"] if is_sex else STYLE["abn_scatter"])

        row = call_df[call_df["chrom"] == chrom].iloc[0]

        if use_cn_scale:
            ax.axhline(2.0, color="blue", linestyle=":", linewidth=1.5, label="Expected (2.0)")
            if not is_sex:
                ax.axhline(3.0, color="red", linestyle=":", linewidth=1, label="Trisomy (3.0)")
                ax.axhline(1.0, color="gray", linestyle=":", linewidth=1, label="Monosomy (1.0)")
        else:
            ax.axhline(0.0, color="blue", linestyle=":", linewidth=1.0, label="Expected (0.0)")

        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Autosomes", f"{chrom}\n({row['call']})"], fontsize=9, fontweight="bold")
        ax.set_title(f"{chrom}  {row.get('detail', '')}", fontsize=8,
                     color=CALL_COLORS.get(row['call'], "black"), fontweight="bold")
        ax.legend(fontsize=6, loc="upper right")
        ax.grid(axis="y", linestyle="--", alpha=0.3)

    if zoom_axes:
        pos = zoom_axes[0].get_position()
        fig.text(0.5, pos.y1 + 0.015, "Bin-level Distribution — Abnormal / Suspicious Chromosomes",
                  ha="center", fontsize=11, fontweight="bold")


def plot_score_heatmap(call_df, sample_id, sex, out_dir):
    score_cols = ["copy_score", "baf_score", "ter_score", "cer_score", "mosaic_score", "final_score"]
    col_labels = ["Copy\n(Signal)", "BAF\nPattern", "TER\n(Trans)", "CER\n(Cis)", "MOSAIC\nIndex", "FINAL"]
    chroms = sort_chroms(call_df["chrom"].tolist())
    mat = call_df.set_index("chrom").reindex(chroms)[score_cols].values.T

    fig, ax = plt.subplots(figsize=(max(14, len(chroms) * 0.8), 6))
    fig.suptitle(f"Anomaly & Mosaic Score Components — {sample_id}  [Sex: {sex}]", fontsize=12, fontweight="bold")
    im = ax.imshow(mat, cmap="RdBu_r", aspect="auto", vmin=-1, vmax=1)

    for i in range(len(score_cols)):
        for j, c in enumerate(chroms):
            val = mat[i, j]
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", color="white" if abs(val) > 0.6 else "black",
                    fontsize=8, fontweight="bold")
            if i == len(score_cols) - 1:
                call = call_df[call_df["chrom"] == c]["call"].values[0]
                lw = 3 if call == "ABNORMAL" else (2 if call == "SUSPICIOUS" else 0)
                if lw:
                    ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, linewidth=lw,
                                                edgecolor=STYLE["abn_shade"] if call == "ABNORMAL" else STYLE["sus_shade"],
                                                facecolor="none"))

    for j, c in enumerate(chroms):
        if _is_sex_chrom(c):
            ax.add_patch(plt.Rectangle((j - 0.5, -0.5), 1, len(score_cols), linewidth=0,
                                        facecolor=STYLE["sex_shade"], alpha=STYLE["sex_shade_alpha"], zorder=0))

    ax.set_xticks(range(len(chroms)))
    ax.set_yticks(range(len(score_cols)))
    ax.set_xticklabels(chroms, rotation=45, ha="right", fontsize=9, fontweight="bold")
    ax.set_yticklabels(col_labels, fontsize=9, fontweight="bold")
    plt.colorbar(im, ax=ax, pad=0.02).set_label("Score Metric Range", fontsize=9)
    plt.tight_layout()
    fig.savefig(os.path.join(out_dir, "03_anomaly_score_heatmap.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_chromosome_overview(df_f, summary, call_df, sex, sample_id, out_dir):
    all_chroms = sort_chroms(df_f["chrom"].unique().tolist())
    chrom_x = {c: i for i, c in enumerate(all_chroms)}

    use_cn_scale = "copy_number_signal" in df_f.columns
    if use_cn_scale:
        cn_col, cn_title, cn_ylim, base_val = "copy_number_signal", "Copy Number\n(Baseline=2.0)", (0.0, 5.0), 2.0
    else:
        cn_col, cn_title, cn_ylim, base_val = "log2_chrom_norm", "Copy Number Signal\n(Log2FC)", (-2.0, 2.0), 0.0

    panels = [
        (cn_col,             cn_title,             "steelblue",       cn_ylim, base_val),
        ("hetero_like_rate", "Hetero-like Rate",   "mediumseagreen",  (0, 1),  0.0),
        ("homo_like_rate",   "Homo-like Rate (LOH)", "tomato",        (0, 1),  0.0),
    ]
    n_top = len(panels)

    # 모든 염색체의 중앙값을 한 번에 미리 계산 (순서 보장)
    chrom_medians = {col: [] for col, _, _, _, _ in panels}
    for chrom in all_chroms:
        grp = df_f[df_f["chrom"] == chrom]
        for col, _, _, _, _ in panels:
            chrom_medians[col].append(grp[col].median())

    abn_chroms = call_df[call_df["call"].isin(["ABNORMAL", "SUSPICIOUS"])]["chrom"].tolist()
    n_zoom = len(abn_chroms)

    fig_h = 3.0 * n_top + (3.4 if n_zoom else 0)
    fig = plt.figure(figsize=(max(18, len(all_chroms) * 0.8), fig_h))
    fig.suptitle(f"Chromosome Overview — {sample_id}  [Sex: {sex}]", fontsize=15, fontweight="bold")

    n_rows = n_top + (1 if n_zoom else 0)
    height_ratios = [1] * n_top + ([1.15] if n_zoom else [])
    gs = fig.add_gridspec(n_rows, max(n_zoom, 1), height_ratios=height_ratios, hspace=0.45, wspace=0.35)

    top_axes = [fig.add_subplot(gs[i, :]) for i in range(n_top)]
    for ax in top_axes[1:]:
        ax.sharex(top_axes[0])

    for i, (ax, (col, title, color, ylim, b_val)) in enumerate(zip(top_axes, panels)):
        if col not in df_f.columns:
            continue

        # 1. 배경 Shade
        for sc in ["chrX", "chrY"]:
            if sc in chrom_x:
                ax.axvspan(chrom_x[sc] - 0.5, chrom_x[sc] + 0.5, alpha=STYLE["sex_shade_alpha"], color=STYLE["sex_shade"], zorder=0)

        # 2. Scatter 데이터 분포 그리기
        for chrom, grp in df_f.groupby("chrom"):
            xi = chrom_x[chrom]
            ax.scatter(xi + np.random.uniform(-0.3, 0.3, len(grp)), grp[col], alpha=0.3, s=10,
                       color=STYLE["sex_scatter"] if _is_sex_chrom(chrom) else color)

        # 3. 중앙값 라인 및 포인트는 루프 외부에서 계산된 리스트로 단 한 번만 그리기
        xs = range(len(all_chroms))
        ys = chrom_medians[col]
        ax.plot(xs, ys, color="black", linewidth=2.0, alpha=0.7, zorder=5)
        ax.scatter(xs, ys, color="black", s=40, zorder=6)

        ax.axhline(b_val, color="gray", linewidth=1.0, linestyle="--", alpha=0.6)

        # 4. Call 결과 Shade
        for _, row in call_df.iterrows():
            if row["chrom"] not in chrom_x: continue
            xi = chrom_x[row["chrom"]]
            if row["call"] == "ABNORMAL":
                ax.axvspan(xi - 0.5, xi + 0.5, alpha=STYLE["abn_shade_alpha"], color=STYLE["abn_shade"], zorder=0)
            elif row["call"] == "SUSPICIOUS":
                ax.axvspan(xi - 0.5, xi + 0.5, alpha=STYLE["sus_shade_alpha"], color=STYLE["sus_shade"], zorder=0)

        ax.set_ylabel(title, fontsize=13, fontweight="bold", labelpad=12)
        ax.tick_params(axis='y', labelsize=11)
        if ylim[1] is not None: ax.set_ylim(ylim)
        ax.grid(axis="y", linestyle="--", alpha=0.4)

    top_axes[-1].set_xticks(range(len(all_chroms)))
    top_axes[-1].set_xticklabels(all_chroms, rotation=45, ha="right", fontsize=13, fontweight="bold")
    for ax in top_axes[:-1]:
        plt.setp(ax.get_xticklabels(), visible=False)

    fig.subplots_adjust(top=0.93)
    fig.savefig(os.path.join(out_dir, "01_chromosome_overview.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)

def plot_final_call(call_df, sex, sample_id, out_dir):
    chroms = sort_chroms(call_df["chrom"].tolist())
    cdf = call_df.set_index("chrom").reindex(chroms).reset_index()

    fig, ax = plt.subplots(figsize=(max(16, len(chroms)), 7))
    sex_colors = STYLE["sex_call_colors"]
    bar_colors = [sex_colors.get(r["call"], "gray") if _is_sex_chrom(r["chrom"]) else CALL_COLORS.get(r["call"], "gray")
                  for _, r in cdf.iterrows()]

    ax.bar(range(len(chroms)), cdf["final_score"].values, color=bar_colors, edgecolor="black", linewidth=0.6, width=0.75)
    ax.axhline(CFG.get("call_thresh_high", 0.7), color="red", linestyle="--", linewidth=1.5)
    ax.axhline(CFG.get("call_thresh_low", 0.5), color="orange", linestyle="--", linewidth=1.2)
    ax.axhline(-CFG.get("call_thresh_high", 0.7), color="red", linestyle="--", linewidth=1.5)
    ax.axhline(-CFG.get("call_thresh_low", 0.5), color="orange", linestyle="--", linewidth=1.2)
    ax.axhline(0, color="black", linewidth=0.8)

    for i, (_, row) in enumerate(cdf.iterrows()):
        if row["call"] != "NORMAL":
            ypos = row["final_score"]
            offset = 0.05 if ypos >= 0 else -0.12
            ax.text(i, ypos + offset, row["call"], ha="center", fontsize=8, fontweight="bold",
                    color=CALL_COLORS.get(row["call"], "black"))
            detail = row.get("detail", "")
            if detail:
                ax.text(i, ypos + offset - 0.13, detail.split("  ")[0], ha="center", fontsize=6.5, color="gray")

    ax.set_xticks(range(len(chroms)))
    ax.set_xticklabels(chroms, rotation=45, ha="right", fontsize=10, fontweight="bold")
    ax.set_ylabel("Final Anomaly Score", fontsize=11, fontweight="bold")
    ax.set_ylim(-1.15, 1.15)
    ax.set_title(f"Final Chromosome Call — {sample_id}  [Sex: {sex}]", fontsize=13, fontweight="bold")
    ax.grid(axis="y", linestyle="--", alpha=0.35)
    plt.tight_layout()
    fig.savefig(os.path.join(out_dir, "04_final_call.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_comprehensive_baf_logr_qc(bins_df, segments_df, output_path, run_id):
    if bins_df is None or segments_df is None or output_path is None or run_id is None:
        raise ValueError("All inputs are mandatory.")

    ALL_CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ['chrX', 'chrY']
    CHROM_MAX_SIZES = {
        'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555,
        'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636,
        'chr9': 138394717, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
        'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345,
        'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
        'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415
    }

    fig, axes = plt.subplots(3, 1, figsize=(28, 14), sharex=True, gridspec_kw={'height_ratios': [2, 2, 1]})
    ax_logr, ax_baf, ax_qc = axes

    use_cn_scale = "copy_number_signal" in bins_df.columns
    cn_col = "copy_number_signal" if use_cn_scale else "log2_chrom_norm"

    current_genome_offset = 0
    chrom_ticks, chrom_names = [], []

    for chrom in ALL_CHROMOSOMES:
        chrom_max_end = CHROM_MAX_SIZES.get(chrom, 100000000)
        group = bins_df[bins_df["chrom"] == chrom].sort_values("start")

        if not group.empty:
            bin_mids = (group["start"].values + group["end"].values) / 2
            abs_x_positions = current_genome_offset + bin_mids

            ax_logr.scatter(abs_x_positions, group[cn_col],
                             s=1.5, c='lightgrey', alpha=0.6, edgecolors='none', zorder=1)

            ax_baf.scatter(abs_x_positions, group["bin_BAF"],
                            s=2.0, c='purple', alpha=0.5, edgecolors='none')

            qc_valid = group[group['total_sites'] >= 5]
            if not qc_valid.empty:
                abs_x_qc = current_genome_offset + (qc_valid["start"].values + qc_valid["end"].values) / 2
                ax_qc.plot(abs_x_qc, qc_valid["imbalance_rate"], color='red', alpha=0.7, linewidth=1.5,
                           label='Imbalance Rate' if chrom == 'chr1' else "")
                ax_qc.plot(abs_x_qc, qc_valid["hetero_like_rate"], color='green', alpha=0.7, linewidth=1.5,
                           label='Hetero Rate' if chrom == 'chr1' else "")

            chrom_max_end = max(chrom_max_end, group["end"].max())

        chrom_segs = segments_df[segments_df["chrom"] == chrom]
        for _, seg in chrom_segs.iterrows():
            seg_start_x = current_genome_offset + seg["start"]
            seg_end_x = current_genome_offset + seg["end"]

            seg_val_log2 = seg.get("seg_median", seg.get("seg_mean", 0.0))

            if use_cn_scale:
                y_val = 2.0 * (2 ** seg_val_log2)
            else:
                y_val = seg_val_log2

            color = 'black'
            if seg.get("copy_number", 2) < 2:
                color = 'blue'
            elif seg.get("copy_number", 2) > 2:
                color = 'red'

            ax_logr.hlines(y=y_val, xmin=seg_start_x, xmax=seg_end_x, colors=color, linewidth=4, zorder=4)

        chrom_ticks.append(current_genome_offset + (chrom_max_end / 2))
        chrom_names.append(chrom.replace("chr", ""))
        current_genome_offset += chrom_max_end

        for ax in axes:
            ax.axvline(x=current_genome_offset, color='black', linestyle='--', linewidth=0.5, alpha=0.6)

    if use_cn_scale:
        ax_logr.set_ylim(-0.5, 5.5)
        ax_logr.axhline(2.0, color='gray', linestyle='-', linewidth=1.0, alpha=0.5)
        ax_logr.axhline(3.0, color='red', linestyle=':', linewidth=1.0, alpha=0.5)
        ax_logr.axhline(1.0, color='blue', linestyle=':', linewidth=1.0, alpha=0.5)
        ax_logr.set_ylabel("Copy Number\n(Baseline=2.0)", fontsize=12)
        ax_logr.set_yticks([0, 1, 2, 3, 4, 5])
    else:
        ax_logr.set_ylim(-2.5, 3.0)
        ax_logr.axhline(0.0, color='gray', linestyle='-', linewidth=1.0, alpha=0.5)
        ax_logr.set_ylabel("Log2 Ratio\n(Normalized Depth)", fontsize=12)

    ax_logr.set_title(f"Comprehensive Multi-omic CNV Profile: {run_id}", fontsize=16, fontweight='bold')

    ax_baf.set_ylim(-0.05, 1.05)
    ax_baf.axhline(0.5, color='gray', linestyle='-', linewidth=1.0, alpha=0.5)
    ax_baf.axhline(0.33, color='orange', linestyle=':', linewidth=1.0, alpha=0.5)
    ax_baf.axhline(0.66, color='orange', linestyle=':', linewidth=1.0, alpha=0.5)
    ax_baf.set_ylabel("B-Allele Frequency\n(Variant Space)", fontsize=12)

    ax_qc.set_ylim(-0.05, 1.05)
    ax_qc.set_ylabel("AI / LOH Rates\n(0 to 1)", fontsize=12)
    ax_qc.set_xlabel("Chromosomes (Physical Position)", fontsize=14)
    ax_qc.legend(loc='upper right')

    ax_qc.set_xticks(chrom_ticks)
    ax_qc.set_xticklabels(chrom_names, fontsize=11, fontweight='bold')
    ax_qc.set_xlim(0, current_genome_offset)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    log(f"Comprehensive 3-panel CNV plot saved to: {output_path}")


# ═════════════════════════════════════════════════════════════
# 다중 샘플 코호트 종합 뷰
# ═════════════════════════════════════════════════════════════
def load_cohort_call_dfs(run_dirs):
    """
    여러 샘플의 파이프라인 실행 결과를 모아서 로드한다.

    run_dirs : dict {sample_id: OutDir}
        OutDir은 CnvPipeline(args).run() 에서 쓰는 args.OutDir 이며,
        run() 안에서 저장한 "<OutDir>/data/<SeqID>.chrom_calls.tsv" 를 읽는다.

    반환값 : dict {sample_id: call_df(DataFrame)}  (파일이 없는 샘플은 제외되고 warning 로그만 남긴다)
    """
    result = {}
    for sample_id, out_dir in run_dirs.items():
        f = os.path.join(out_dir, "data", f"{sample_id}.chrom_calls.tsv")
        if not os.path.exists(f):
            log(f"[Cohort] Skip {sample_id}: call_df not found at {f}")
            continue
        result[sample_id] = pd.read_csv(f, sep="\t")
    return result


def plot_cohort_overview(call_df_dict, out_path, score_col="final_score", sex_dict=None):
    """
    여러 샘플의 call_df를 하나의 heatmap + summary bar로 종합 시각화.

    call_df_dict : dict {sample_id: call_df}   (각 call_df는 최소 chrom/call/<score_col> 컬럼 필요)
    out_path     : 저장 경로 (.png)
    score_col    : heatmap에 표시할 점수 컬럼 (기본 final_score)
    sex_dict     : dict {sample_id: 'XX'/'XY'/'Unknown'} (선택, y축 라벨에 성별 표기)
    """
    if not call_df_dict:
        log("[Cohort] No call_df provided, skip cohort plot.")
        return

    sample_ids = list(call_df_dict.keys())
    all_chroms = sort_chroms(sorted(set(c for df in call_df_dict.values() for c in df["chrom"])))

    mat = np.full((len(sample_ids), len(all_chroms)), np.nan)
    call_mat = np.full((len(sample_ids), len(all_chroms)), "", dtype=object)
    abn_counts = {"ABNORMAL": [], "SUSPICIOUS": []}

    for i, sid in enumerate(sample_ids):
        df = call_df_dict[sid].set_index("chrom")
        n_abn, n_sus = 0, 0
        for j, c in enumerate(all_chroms):
            if c in df.index:
                mat[i, j] = df.loc[c, score_col]
                call = df.loc[c, "call"]
                call_mat[i, j] = call
                if call == "ABNORMAL":
                    n_abn += 1
                elif call == "SUSPICIOUS":
                    n_sus += 1
        abn_counts["ABNORMAL"].append(n_abn)
        abn_counts["SUSPICIOUS"].append(n_sus)

    fig_w = max(14, len(all_chroms) * 0.8 + 3)
    fig_h = max(4, len(sample_ids) * 0.55) + 2.5
    fig = plt.figure(figsize=(fig_w, fig_h))
    fig.suptitle(f"Cohort Overview — {len(sample_ids)} samples  ({score_col})", fontsize=15, fontweight="bold")

    gs = fig.add_gridspec(2, 1, height_ratios=[max(2, len(sample_ids) * 0.55), 1.4], hspace=0.55)
    ax_hm = fig.add_subplot(gs[0])
    ax_bar = fig.add_subplot(gs[1])

    im = ax_hm.imshow(mat, cmap="RdBu_r", aspect="auto", vmin=-1, vmax=1)

    for i in range(len(sample_ids)):
        for j in range(len(all_chroms)):
            val = mat[i, j]
            if np.isnan(val):
                continue
            ax_hm.text(j, i, f"{val:.2f}", ha="center", va="center",
                       color="white" if abs(val) > 0.6 else "black", fontsize=7)
            call = call_mat[i, j]
            if call == "ABNORMAL":
                ax_hm.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                               edgecolor=STYLE["abn_shade"], linewidth=2.2))
            elif call == "SUSPICIOUS":
                ax_hm.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                               edgecolor=STYLE["sus_shade"], linewidth=1.6))

    for j, c in enumerate(all_chroms):
        if _is_sex_chrom(c):
            ax_hm.add_patch(plt.Rectangle((j - 0.5, -0.5), 1, len(sample_ids), linewidth=0,
                                           facecolor=STYLE["sex_shade"], alpha=STYLE["sex_shade_alpha"], zorder=0))

    ax_hm.set_xticks(range(len(all_chroms)))
    ax_hm.set_xticklabels(all_chroms, rotation=45, ha="right", fontsize=9, fontweight="bold")
    ax_hm.set_yticks(range(len(sample_ids)))
    y_labels = [f"{sid}  [{sex_dict.get(sid, '?')}]" for sid in sample_ids] if sex_dict else sample_ids
    ax_hm.set_yticklabels(y_labels, fontsize=9, fontweight="bold")
    plt.colorbar(im, ax=ax_hm, pad=0.02).set_label(score_col, fontsize=9)

    # 샘플별 ABNORMAL/SUSPICIOUS 개수 요약 bar
    y_pos = np.arange(len(sample_ids))
    ax_bar.barh(y_pos, abn_counts["ABNORMAL"], color=STYLE["abn_shade"], alpha=0.85, label="ABNORMAL")
    ax_bar.barh(y_pos, abn_counts["SUSPICIOUS"], left=abn_counts["ABNORMAL"],
                color=STYLE["sus_shade"], alpha=0.85, label="SUSPICIOUS")
    ax_bar.set_yticks(y_pos)
    ax_bar.set_yticklabels(y_labels, fontsize=9, fontweight="bold")
    ax_bar.invert_yaxis()
    ax_bar.set_xlabel("Chromosome Call Count", fontsize=10, fontweight="bold")
    ax_bar.set_title("Abnormal / Suspicious Calls per Sample", fontsize=11, fontweight="bold")
    ax_bar.legend(fontsize=8, loc="lower right")
    ax_bar.grid(axis="x", linestyle="--", alpha=0.35)

    fig.subplots_adjust(top=0.93)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log(f"Cohort overview plot saved to: {out_path}")