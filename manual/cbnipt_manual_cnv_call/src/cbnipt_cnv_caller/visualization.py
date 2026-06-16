import matplotlib.pyplot as plt
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



def plot_gc_correction(gc_stats, out_path):
    if gc_stats is None: return
    x, y, fit = gc_stats
    plt.figure(figsize=(8, 5))
    plt.scatter(x, y, s=1, color='gray', alpha=0.2, label="Raw")
    idx = x.argsort()
    plt.plot(x[idx], fit[idx], color='red', lw=2, label="Lowess Fit")
    plt.title("GC Bias Correction (LOWESS)"); plt.legend(); plt.savefig(out_path); plt.close()


def plot_bin_distribution(df_f, call_df, sample_id, out_dir):
    # Copy Number (2.0) 기준으로 스케일 변경
    use_cn_scale = "copy_number_signal" in df_f.columns
    cn_col = "copy_number_signal" if use_cn_scale else "log2_chrom_norm"
    
    if cn_col not in df_f.columns: return
    abn_chroms = call_df[call_df["call"].isin(["ABNORMAL", "SUSPICIOUS"])]["chrom"].tolist()
    if not abn_chroms: return

    n = len(abn_chroms)
    fig, axes = plt.subplots(1, n, figsize=(max(6, 5 * n), 6), squeeze=False)
    fig.suptitle(f"Bin Distribution — Abnormal/Suspicious ({sample_id})", fontsize=12, fontweight="bold")
    auto_vals = df_f[~df_f["chrom"].isin(["chrX", "chrY"])][cn_col].dropna().values

    for ax, chrom in zip(axes[0], abn_chroms):
        chrom_vals = df_f[df_f["chrom"] == chrom][cn_col].dropna().values
        parts = ax.violinplot([auto_vals, chrom_vals], positions=[0, 1], showmedians=True, showextrema=True)
        for body, vc in zip(parts["bodies"], ["#90CAF9", "#CE93D8" if chrom in ("chrX","chrY") else "#EF9A9A"]):
            body.set_facecolor(vc); body.set_alpha(0.7)

        ax.scatter([0]*len(auto_vals), auto_vals, alpha=0.1, s=5, color="#1565C0")
        ax.scatter([1]*len(chrom_vals), chrom_vals, alpha=0.4, s=8, color="#6A1B9A" if chrom in ("chrX","chrY") else "#B71C1C")

        row = call_df[call_df["chrom"] == chrom].iloc[0]
        
        # 2.0 Baseline 가이드라인
        if use_cn_scale:
            ax.axhline(2.0, color="blue", linestyle=":", linewidth=1.5, label="Expected (2.0)")
            if chrom not in ("chrX", "chrY"):
                ax.axhline(3.0, color="red", linestyle=":", linewidth=1, label="Trisomy (3.0)")
                ax.axhline(1.0, color="gray", linestyle=":", linewidth=1, label="Monosomy (1.0)")
        else:
            ax.axhline(0.0, color="blue", linestyle=":", linewidth=1.0, label="Expected (0.0)")

        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Autosomes", f"{chrom}\n({row['call']})"], fontsize=10, fontweight="bold")
        ax.set_title(f"{chrom}\n{row.get('detail','')}", fontsize=9, color=CALL_COLORS.get(row['call'], "black"), fontweight="bold")
        ax.legend(fontsize=7, loc="upper right")
        ax.grid(axis="y", linestyle="--", alpha=0.3)

    plt.tight_layout()
    fig.savefig(os.path.join(out_dir, "02_bin_distribution.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


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
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", color="white" if abs(val)>0.6 else "black", fontsize=8, fontweight="bold")
            if i == len(score_cols) - 1:
                call = call_df[call_df["chrom"] == c]["call"].values[0]
                lw = 3 if call == "ABNORMAL" else (2 if call == "SUSPICIOUS" else 0)
                if lw: ax.add_patch(plt.Rectangle((j-0.5, i-0.5), 1, 1, linewidth=lw, edgecolor="red" if call=="ABNORMAL" else "orange", facecolor="none"))

    for j, c in enumerate(chroms):
        if c in ("chrX", "chrY"): ax.add_patch(plt.Rectangle((j-0.5, -0.5), 1, len(score_cols), linewidth=0, facecolor="purple", alpha=0.08, zorder=0))

    ax.set_xticks(range(len(chroms))); ax.set_yticks(range(len(score_cols)))
    ax.set_xticklabels(chroms, rotation=45, ha="right", fontsize=9, fontweight="bold")
    ax.set_yticklabels(col_labels, fontsize=9, fontweight="bold")
    plt.colorbar(im, ax=ax, pad=0.02).set_label("Score Metric Range", fontsize=9)
    plt.tight_layout()
    fig.savefig(os.path.join(out_dir, "03_anomaly_score_heatmap.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_chromosome_overview(df_f, summary, call_df, sex, sample_id, out_dir):
    all_chroms = sort_chroms(df_f["chrom"].unique().tolist())
    chrom_x    = {c: i for i, c in enumerate(all_chroms)}
    
    use_cn_scale = "copy_number_signal" in df_f.columns
    if use_cn_scale:
        cn_col, cn_title, cn_ylim, base_val = "copy_number_signal", "Copy Number\n(Baseline=2.0)", (0.0, 5.0), 2.0
    else:
        cn_col, cn_title, cn_ylim, base_val = "log2_chrom_norm", "Copy Number Signal\n(Log2FC)", (-2.0, 2.0), 0.0
    
    # [수정 사항] Trans, Cis 패널을 날리고 핵심 3개 패널만 직관적으로 표출
    panels = [
        (cn_col,             cn_title,               "steelblue",      cn_ylim, base_val),
        ("hetero_like_rate", "Hetero-like Rate",       "mediumseagreen", (0, 1), 0.0),
        ("homo_like_rate",   "Homo-like Rate (LOH)",   "tomato",         (0, 1), 0.0),
    ]

    fig, axes = plt.subplots(len(panels), 1, figsize=(max(18, len(all_chroms) * 0.8), 3.0 * len(panels)), sharex=True)
    fig.suptitle(f"Chromosome Overview — {sample_id}  [Sex: {sex}]", fontsize=15, fontweight="bold")

    for ax, (col, title, color, ylim, b_val) in zip(axes, panels):
        if col not in df_f.columns: continue
        for sc in ["chrX", "chrY"]:
            if sc in chrom_x: ax.axvspan(chrom_x[sc] - 0.5, chrom_x[sc] + 0.5, alpha=0.10, color="purple", zorder=0)

        xs, ys = [], []
        for chrom, grp in df_f.groupby("chrom"):
            xi = chrom_x[chrom]
            ax.scatter(xi + np.random.uniform(-0.3, 0.3, len(grp)), grp[col], alpha=0.3, s=10, color="purple" if chrom in ("chrX", "chrY") else color)
            # 동적으로 중앙값 선을 그리기 위한 데이터 적재
            xs.append(xi)
            ys.append(grp[col].median())

        # 중앙값 라인 그리기
        ax.plot(xs, ys, color="black", linewidth=2.0, alpha=0.7, zorder=5)
        ax.scatter(xs, ys, color="black", s=40, zorder=6)
        
        # 중심 가이드 라인
        ax.axhline(b_val, color="gray", linewidth=1.0, linestyle="--", alpha=0.6)

        for _, row in call_df.iterrows():
            if row["chrom"] not in chrom_x: continue
            xi = chrom_x[row["chrom"]]
            if row["call"] == "ABNORMAL": ax.axvspan(xi - 0.5, xi + 0.5, alpha=0.18, color="red", zorder=0)
            elif row["call"] == "SUSPICIOUS": ax.axvspan(xi - 0.5, xi + 0.5, alpha=0.12, color="orange", zorder=0)

        ax.set_ylabel(title, fontsize=13, fontweight="bold", labelpad=12)
        ax.tick_params(axis='y', labelsize=11)
        if ylim[1] is not None: ax.set_ylim(ylim)
        ax.grid(axis="y", linestyle="--", alpha=0.4)

    axes[-1].set_xticks(range(len(all_chroms)))
    axes[-1].set_xticklabels(all_chroms, rotation=45, ha="right", fontsize=13, fontweight="bold")
    
    plt.tight_layout()
    fig.subplots_adjust(top=0.92) 
    fig.savefig(os.path.join(out_dir, "01_chromosome_overview.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_final_call(call_df, sex, sample_id, out_dir):
    chroms = sort_chroms(call_df["chrom"].tolist())
    cdf    = call_df.set_index("chrom").reindex(chroms).reset_index()

    fig, ax = plt.subplots(figsize=(max(16, len(chroms)), 7))
    bar_colors = [{"NORMAL":"#AB47BC","SUSPICIOUS":"#FF9800","ABNORMAL":"#F44336"}.get(r["call"], "gray") if r["chrom"] in ("chrX", "chrY") else CALL_COLORS.get(r["call"], "gray") for _, r in cdf.iterrows()]
    
    ax.bar(range(len(chroms)), cdf["final_score"].values, color=bar_colors, edgecolor="black", linewidth=0.6, width=0.75)
    ax.axhline(CFG.get("call_thresh_high", 0.7), color="red", linestyle="--", linewidth=1.5)
    ax.axhline(CFG.get("call_thresh_low", 0.5),  color="orange", linestyle="--", linewidth=1.2)
    ax.axhline(-CFG.get("call_thresh_high", 0.7), color="red", linestyle="--", linewidth=1.5)
    ax.axhline(-CFG.get("call_thresh_low", 0.5),  color="orange", linestyle="--", linewidth=1.2)
    ax.axhline(0, color="black", linewidth=0.8)

    for i, (_, row) in enumerate(cdf.iterrows()):
        if row["call"] != "NORMAL":
            ypos = row["final_score"]
            offset = 0.05 if ypos >= 0 else -0.12
            ax.text(i, ypos + offset, row["call"], ha="center", fontsize=8, fontweight="bold", color=CALL_COLORS.get(row["call"], "black"))
            detail = row.get("detail", "")
            if detail: ax.text(i, ypos + offset - 0.13, detail.split("  ")[0], ha="center", fontsize=6.5, color="gray")

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

    # Copy Number (2.0) 기준으로 스케일 변경
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
            
            # [Panel 1] 실제 Copy Number 점찍기
            ax_logr.scatter(abs_x_positions, group[cn_col], 
                            s=1.5, c='lightgrey', alpha=0.6, edgecolors='none', zorder=1)
            
            # [Panel 2] BAF 
            ax_baf.scatter(abs_x_positions, group["bin_BAF"], 
                           s=2.0, c='purple', alpha=0.5, edgecolors='none')
            
            # [Panel 3] QC & Imbalance
            qc_valid = group[group['total_sites'] >= 5]
            if not qc_valid.empty:
                abs_x_qc = current_genome_offset + (qc_valid["start"].values + qc_valid["end"].values) / 2
                ax_qc.plot(abs_x_qc, qc_valid["imbalance_rate"], color='red', alpha=0.7, linewidth=1.5, label='Imbalance Rate' if chrom=='chr1' else "")
                ax_qc.plot(abs_x_qc, qc_valid["hetero_like_rate"], color='green', alpha=0.7, linewidth=1.5, label='Hetero Rate' if chrom=='chr1' else "")

            chrom_max_end = max(chrom_max_end, group["end"].max())

        # [Panel 1] 세그먼트 오버레이 (안전한 get 사용)
        chrom_segs = segments_df[segments_df["chrom"] == chrom]
        for _, seg in chrom_segs.iterrows():
            seg_start_x = current_genome_offset + seg["start"]
            seg_end_x = current_genome_offset + seg["end"]
            
            # seg_mean이나 seg_median 이름 꼬임 방지
            seg_val_log2 = seg.get("seg_median", seg.get("seg_mean", 0.0))
            
            if use_cn_scale:
                # Log2 공간에 있는 세그먼트 값을 실제 Copy Number로 복원
                y_val = 2.0 * (2 ** seg_val_log2)
            else:
                y_val = seg_val_log2
            
            color = 'black'
            if seg.get("copy_number", 2) < 2: color = 'blue'
            elif seg.get("copy_number", 2) > 2: color = 'red'
            
            ax_logr.hlines(y=y_val, xmin=seg_start_x, xmax=seg_end_x, colors=color, linewidth=4, zorder=4)

        chrom_ticks.append(current_genome_offset + (chrom_max_end / 2))
        chrom_names.append(chrom.replace("chr", ""))
        current_genome_offset += chrom_max_end
        
        for ax in axes:
            ax.axvline(x=current_genome_offset, color='black', linestyle='--', linewidth=0.5, alpha=0.6)

    # --- [Panel 1] 데코레이션 (Copy Number Scale) ---
    if use_cn_scale:
        ax_logr.set_ylim(-0.5, 5.5)
        ax_logr.axhline(2.0, color='gray', linestyle='-', linewidth=1.0, alpha=0.5)
        ax_logr.axhline(3.0, color='red', linestyle=':', linewidth=1.0, alpha=0.5) # Trisomy
        ax_logr.axhline(1.0, color='blue', linestyle=':', linewidth=1.0, alpha=0.5) # Monosomy
        ax_logr.set_ylabel("Copy Number\n(Baseline=2.0)", fontsize=12)
        ax_logr.set_yticks([0, 1, 2, 3, 4, 5])
    else:
        ax_logr.set_ylim(-2.5, 3.0)
        ax_logr.axhline(0.0, color='gray', linestyle='-', linewidth=1.0, alpha=0.5)
        ax_logr.set_ylabel("Log2 Ratio\n(Normalized Depth)", fontsize=12)

    ax_logr.set_title(f"Comprehensive Multi-omic CNV Profile: {run_id}", fontsize=16, fontweight='bold')
    
    # --- [Panel 2] 데코레이션 (BAF) ---
    ax_baf.set_ylim(-0.05, 1.05)
    ax_baf.axhline(0.5, color='gray', linestyle='-', linewidth=1.0, alpha=0.5) 
    ax_baf.axhline(0.33, color='orange', linestyle=':', linewidth=1.0, alpha=0.5) 
    ax_baf.axhline(0.66, color='orange', linestyle=':', linewidth=1.0, alpha=0.5) 
    ax_baf.set_ylabel("B-Allele Frequency\n(Variant Space)", fontsize=12)
    
    # --- [Panel 3] 데코레이션 (Rates) ---
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