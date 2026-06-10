"""
Single-Sample Chromosome Abnormality Detector (v2.2 - Sex Chromosome Ratio Hard-Cutoffs)
=============================================================================
핵심 업데이트:
  1. 임상 표준 Ratio 기반 성염색체 이상(SCA) 하드 컷오프 적용 (Turner, XXX, XXY, XYY)
  2. chrY 존재 여부(y_ratio)로 1차 성별 확정 후, 성별에 따른 교차 비율 검증
"""

import sys, os, argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
from pathlib import Path

# ═══════════════════════════════════════════════════════════════
# 0. 전역 구성 매개변수
# ═══════════════════════════════════════════════════════════════

CFG = dict(
    min_depth            = 1,
    min_coverage         = 0.5,

    # [핵심 업데이트] 성별 및 성염색체 질환 진단 컷오프
    y_presence_threshold = 0.15,   # XY 판단 기준 (chrY/autosome)
    
    sex_mono_x_ratio     = 0.75,   # XX에서 chrX/auto <= 0.75 -> Turner(X0) 의심
    sex_xxx_ratio        = 1.25,   # XX에서 chrX/auto >= 1.25 -> XXX 의심
    sex_xxy_ratio        = 0.75,   # XY에서 chrX/auto >= 0.75 -> Klinefelter(XXY) 의심
    sex_mono_y_ratio     = 0.25,   # XY에서 chrY/auto <= 0.25 -> chrY loss 의심
    sex_xyy_ratio        = 0.75,   # XY에서 chrY/auto >= 0.75 -> XYY 의심

    robust_z_thresh      = 2.5,
    pval_thresh          = 0.05,
    homo_rate_thresh     = 0.6,
    ter_z_thresh         = 2.0,

    sex_log2_abnormal    = 0.5,   
    sex_log2_suspicious  = 0.28,   
    y_noise_threshold    = 0.10,   

    mosaic_log2_min      = 0.10,
    mosaic_log2_max      = 0.45,
    mosaic_bin_shift_min = 0.55,
    mosaic_score_thresh  = 0.45,

    w_copy  = 0.35,
    w_baf   = 0.25,
    w_ter   = 0.15,
    w_cer   = 0.15,
    w_pval  = 0.10,

    call_thresh_high = 0.50,
    call_thresh_low  = 0.25,
)

CALL_COLORS = {"NORMAL": "#4CAF50", "SUSPICIOUS": "#FF9800", "ABNORMAL": "#F44336"}

# ═══════════════════════════════════════════════════════════════
# 1. 유틸리티 및 전처리
# ═══════════════════════════════════════════════════════════════

def chrom_key(c):
    s = str(c).replace("chr", "")
    if s == "X": return 23
    if s == "Y": return 24
    try: return int(s)
    except: return 99

def sort_chroms(chroms):
    return sorted(chroms, key=chrom_key)

def safe_log2fc(val, ref):
    return float(np.log2((val + 1e-9) / (ref + 1e-9)))

def robust_z_single(val, others_arr):
    others = np.array([v for v in others_arr if not np.isnan(v)])
    if len(others) < 2: return np.nan, np.nan
    med = np.median(others)
    mad = max(np.median(np.abs(others - med)), 1e-9)
    return safe_log2fc(val, med), (val - med) / (1.4826 * mad)

def apply_qc_filter(df):
    df_f = df[(df["raw_count"] >= CFG["min_depth"]) & (df["breadth_ratio"] >= CFG["min_coverage"])].copy()
    df_f["_k"] = df_f["chrom"].apply(chrom_key)
    return df_f.sort_values(["_k", "start"]).drop(columns=["_k"]).reset_index(drop=True)

# ═══════════════════════════════════════════════════════════════
# 3. [수정] 성별 추정 및 X/Y Ratio 추출
# ═══════════════════════════════════════════════════════════════

def estimate_sex(df_f):
    """chrY 비율로 성별을 확정하고, X와 Y의 비율(Ratio)을 명시적으로 반환합니다."""
    auto_df = df_f[~df_f["chrom"].isin(["chrX", "chrY"])]
    x_df    = df_f[df_f["chrom"] == "chrX"]
    y_df    = df_f[df_f["chrom"] == "chrY"]
    
    cn_col = "log2_chrom_norm" if "log2_chrom_norm" in df_f.columns else "raw_count"
    auto_med = auto_df[cn_col].median() if len(auto_df) > 0 else np.nan

    if np.isnan(auto_med) or auto_med < 1e-6:
        return "XX", 0.0, 0.0, auto_med

    x_med = x_df[cn_col].median() if len(x_df) > 0 else 0.0
    y_med = y_df[cn_col].median() if len(y_df) > 0 else 0.0
    
    x_ratio = x_med / auto_med
    y_ratio = y_med / auto_med

    # chrY가 threshold 이상이면 남성(XY)으로 확정
    sex = "XY" if y_ratio >= CFG["y_presence_threshold"] else "XX"
    return sex, x_ratio, y_ratio, auto_med

def compute_chrom_summary(df):
    cols = ["raw_count", "breadth_ratio", "hetero_like_rate", "homo_like_rate", "imbalance_rate", "raw_bin_TER", "raw_bin_CER", "bin_BAF", "MAD_BAF", "log2_chrom_norm", "Informative_OE_Ratio"]
    cols = [c for c in cols if c in df.columns]
    rows = []
    for chrom, grp in df.groupby("chrom"):
        row = {"chrom": chrom, "n_bins": len(grp)}
        for c in cols:
            row[f"{c}_median"], row[f"{c}_mean"], row[f"{c}_std"] = grp[c].median(), grp[c].mean(), grp[c].std()
        rows.append(row)
    summary = pd.DataFrame(rows)
    summary["_k"] = summary["chrom"].apply(chrom_key)
    return summary.sort_values("_k").drop(columns=["_k"]).reset_index(drop=True)

# ═══════════════════════════════════════════════════════════════
# 4. 통계 연산 (상염색체)
# ═══════════════════════════════════════════════════════════════

def compute_robust_loo(summary, value_col):
    if value_col not in summary.columns: return {}
    auto_summary = summary[~summary["chrom"].isin(["chrX", "chrY"])]
    chroms, vals = auto_summary["chrom"].tolist(), auto_summary[value_col].values
    result = {}
    for i, c in enumerate(chroms):
        others = np.array([vals[j] for j in range(len(vals)) if j != i and not np.isnan(vals[j])])
        result[c] = (np.nan, np.nan) if len(others) < 3 else robust_z_single(vals[i], others)
    return result

def compute_bin_distribution_test(df, metric):
    auto_df = df[~df["chrom"].isin(["chrX", "chrY"])]
    if metric not in auto_df.columns: return {}
    result = {}
    for chrom, grp in auto_df.groupby("chrom"):
        c_vals = grp[metric].dropna().values
        o_vals = auto_df[auto_df["chrom"] != chrom][metric].dropna().values
        if len(c_vals) < 3 or len(o_vals) < 3:
            result[chrom] = (np.nan, np.nan, 0)
        else:
            stat, pval = stats.mannwhitneyu(c_vals, o_vals, alternative="two-sided")
            result[chrom] = (stat, pval, 1 if c_vals.mean() > o_vals.mean() else -1)
    return result

def compute_baf_signal(summary):
    result = {}
    for _, row in summary.iterrows():
        c, scores, n = row["chrom"], 0.0, 0
        baf = row.get("bin_BAF_median", np.nan)
        if not np.isnan(baf) and c not in ("chrX", "chrY"):
            scores += min(abs(baf - 0.5) / 0.17, 1.0); n += 1
        homo = row.get("homo_like_rate_mean", np.nan)
        if not np.isnan(homo):
            scores += min(max(homo - CFG["homo_rate_thresh"], 0) / (1 - CFG["homo_rate_thresh"]), 1.0); n += 1
        imbal = row.get("imbalance_rate_mean", np.nan)
        if not np.isnan(imbal) and c not in ("chrX", "chrY"):
            scores += min(imbal * 2, 1.0); n += 1
        result[c] = scores / n if n > 0 else 0.0
    return result

def compute_mosaic_score(log2fc, robust_z, baf_score, ter_score, cer_score, pval_score):
    if np.isnan(log2fc): return 0.0
    abs_l2 = abs(log2fc)
    c_comp = 0.0 if abs_l2 < CFG["mosaic_log2_min"] else (
        (abs_l2 - CFG["mosaic_log2_min"]) / (CFG["mosaic_log2_max"] - CFG["mosaic_log2_min"]) if abs_l2 <= CFG["mosaic_log2_max"] 
        else max(0.0, 1.0 - (abs_l2 - CFG["mosaic_log2_max"]) / 0.40)
    )
    rz_comp = 0.0 if np.isnan(robust_z) else min(abs(robust_z) / CFG["robust_z_thresh"], 1.0)
    support = 0.40*c_comp + 0.20*rz_comp + 0.20*abs(baf_score) + 0.05*abs(ter_score) + 0.05*abs(cer_score) + 0.10*abs(pval_score)
    return float(np.clip(support, 0, 1))

def make_autosome_calls(chroms, loo_cn, bin_pvals, baf_sig, loo_ter, loo_cer):
    records = []
    for c in chroms:
        log2fc, rz = loo_cn.get(c, (np.nan, np.nan))
        copy_s = float(np.clip(rz / (CFG["robust_z_thresh"] * 2), -1, 1)) if not np.isnan(rz) else 0.0
        baf_s = baf_sig.get(c, 0.0)
        ter_s = float(np.clip(loo_ter.get(c, (0, np.nan))[1] / (CFG["ter_z_thresh"] * 2), -1, 1)) if not np.isnan(loo_ter.get(c, (0, np.nan))[1]) else 0.0
        cer_s = float(np.clip(loo_cer.get(c, (0, np.nan))[1] / (CFG["ter_z_thresh"] * 2), -1, 1)) if not np.isnan(loo_cer.get(c, (0, np.nan))[1]) else 0.0
        pval_s = 0.0
        if c in bin_pvals:
            _, pval, d = bin_pvals[c]
            if not np.isnan(pval) and pval < CFG["pval_thresh"]:
                pval_s = float(np.clip(-np.log10(pval + 1e-10) / 5, 0, 1)) * d

        final = CFG["w_copy"]*copy_s + CFG["w_baf"]*baf_s + CFG["w_ter"]*ter_s + CFG["w_cer"]*cer_s + CFG["w_pval"]*pval_s
        mosaic_s = compute_mosaic_score(log2fc, rz, baf_s, ter_s, cer_s, pval_s)
        call = "ABNORMAL" if abs(final) >= CFG["call_thresh_high"] else ("SUSPICIOUS" if abs(final) >= CFG["call_thresh_low"] else "NORMAL")

        records.append(dict(
            chrom=c, log2fc=log2fc, robust_z=rz, expected_copy="2 copy",
            copy_score=copy_s, baf_score=baf_s, ter_score=ter_s, cer_score=cer_s, pval_score=pval_s,
            mosaic_score=mosaic_s, final_score=final, call=call,
            detail=(f"Log2FC={log2fc:+.3f} RZ={rz:+.2f} Mos={mosaic_s:.2f}" if not np.isnan(log2fc) else ""), sex_note="",
        ))
    return records

# ═══════════════════════════════════════════════════════════════
# 7. [핵심 수정] 성염색체 질환 규칙 강제 적용 엔진
# ═══════════════════════════════════════════════════════════════

def analyze_sex_chromosomes(df_f, summary, sex, x_ratio, y_ratio, auto_med):
    """
    연구원님의 Ratio 규칙을 최우선으로 검사하여 질환명(Type)과 Call을 강제 할당합니다.
    """
    results = []
    
    # 내부 통계 함수들
    def get_auto_mad(metric):
        vals = summary[~summary["chrom"].isin(["chrX", "chrY"])][metric].dropna().values
        if len(vals) < 2: return 1e-9
        return max(np.median(np.abs(vals - np.median(vals))), 1e-9)

    def get_rz(chrom, metric):
        val = summary.loc[summary["chrom"]==chrom, metric].median() if not summary[summary["chrom"]==chrom].empty else np.nan
        if np.isnan(val): return 0.0
        vals = summary[~summary["chrom"].isin(["chrX", "chrY"])][metric].dropna().values
        if len(vals) < 2: return 0.0
        med = np.median(vals)
        mad = max(np.median(np.abs(vals - med)), 1e-9)
        return float(np.clip(((val - med) / (1.4826 * mad)) / (CFG["ter_z_thresh"] * 2), -1, 1))

    # 기본 스코어 연산 (X)
    cn_col = "log2_chrom_norm" if "log2_chrom_norm" in df_f.columns else "raw_count"
    x_val = summary.loc[summary["chrom"]=="chrX", f"{cn_col}_median"].values[0] if not summary[summary["chrom"]=="chrX"].empty else np.nan
    
    ref_x = auto_med if sex == "XX" else auto_med * 0.5
    expected_x = "2 copy" if sex == "XX" else "1 copy"
    note_x = f"Expected ~{ref_x:.2f}"
    
    log2fc_x = safe_log2fc(x_val, ref_x) if not np.isnan(x_val) else np.nan
    rz_x = (x_val - ref_x) / (1.4826 * get_auto_mad(f"{cn_col}_median")) if not np.isnan(x_val) else np.nan
    copy_s_x = float(np.clip(log2fc_x / 1.0, -1, 1)) if not np.isnan(log2fc_x) else 0.0
    ter_s_x = get_rz("chrX", "raw_bin_TER_median")
    cer_s_x = get_rz("chrX", "raw_bin_CER_median")
    final_x = CFG["w_copy"] * copy_s_x + CFG["w_ter"] * ter_s_x + CFG["w_cer"] * cer_s_x

    # ────────────────────────────────────────────────────────
    # [연구원님 룰 적용] Ratio 기반 질환 추정 (chrX)
    # ────────────────────────────────────────────────────────
    call_x = "NORMAL"
    type_x = ""
    
    if sex == "XX":
        if x_ratio <= CFG["sex_mono_x_ratio"]:
            call_x, type_x = "ABNORMAL", "Turner(X0) 의심"
            final_x = min(final_x - 0.5, -0.6) # 시각화를 위해 점수 강제 조정
        elif x_ratio >= CFG["sex_xxx_ratio"]:
            call_x, type_x = "ABNORMAL", "XXX 의심"
            final_x = max(final_x + 0.5, 0.6)
        else:
            if abs(log2fc_x) >= CFG["sex_log2_abnormal"]: call_x = "ABNORMAL"
            elif abs(log2fc_x) >= CFG["sex_log2_suspicious"]: call_x = "SUSPICIOUS"
            
    elif sex == "XY":
        if x_ratio >= CFG["sex_xxy_ratio"]:
            call_x, type_x = "ABNORMAL", "Klinefelter(XXY) 의심"
            final_x = max(final_x + 0.5, 0.6)
        elif x_ratio <= 0.25:
            call_x, type_x = "ABNORMAL", "chrX loss 의심"
            final_x = min(final_x - 0.5, -0.6)
        else:
            if abs(log2fc_x) >= CFG["sex_log2_abnormal"]: call_x = "ABNORMAL"
            elif abs(log2fc_x) >= CFG["sex_log2_suspicious"]: call_x = "SUSPICIOUS"

    results.append(dict(
        chrom="chrX", log2fc=log2fc_x, robust_z=rz_x, expected_copy=expected_x,
        copy_score=copy_s_x, baf_score=0.0, ter_score=ter_s_x, cer_score=cer_s_x, pval_score=0.0,
        mosaic_score=0.0, final_score=final_x, call=call_x,
        detail=f"X_Ratio={x_ratio:.2f} {type_x}", sex_note=note_x,
    ))

    # 기본 스코어 연산 (Y)
    y_val = summary.loc[summary["chrom"]=="chrY", f"{cn_col}_median"].values[0] if not summary[summary["chrom"]=="chrY"].empty else np.nan
    ref_y = auto_med * 0.5 if sex == "XY" else np.nan
    expected_y = "1 copy" if sex == "XY" else "0 copy"
    note_y = f"Expected ~{ref_y:.2f}" if sex == "XY" else "Should be absent"

    log2fc_y = safe_log2fc(y_val, ref_y) if sex == "XY" and not np.isnan(y_val) else np.nan
    rz_y = (y_val - ref_y) / (1.4826 * get_auto_mad(f"{cn_col}_median")) if sex == "XY" and not np.isnan(y_val) else np.nan
    copy_s_y = float(np.clip(log2fc_y / 1.0, -1, 1)) if not np.isnan(log2fc_y) else 0.0
    ter_s_y = get_rz("chrY", "raw_bin_TER_median")
    cer_s_y = get_rz("chrY", "raw_bin_CER_median")
    final_y = CFG["w_copy"] * copy_s_y + CFG["w_ter"] * ter_s_y + CFG["w_cer"] * cer_s_y

    # ────────────────────────────────────────────────────────
    # [연구원님 룰 적용] Ratio 기반 질환 추정 (chrY)
    # ────────────────────────────────────────────────────────
    call_y = "NORMAL"
    type_y = ""

    if sex == "XY":
        if np.isnan(y_val):
            call_y, type_y, final_y = "SUSPICIOUS", "Y data missing", -0.5
        elif y_ratio >= CFG["sex_xyy_ratio"]:
            call_y, type_y = "ABNORMAL", "XYY 의심"
            final_y = max(final_y + 0.5, 0.6)
        elif y_ratio <= CFG["sex_mono_y_ratio"]:
            call_y, type_y = "ABNORMAL", "mono Y / chrY loss 의심"
            final_y = min(final_y - 0.5, -0.6)
        else:
            if abs(log2fc_y) >= CFG["sex_log2_abnormal"]: call_y = "ABNORMAL"
            elif abs(log2fc_y) >= CFG["sex_log2_suspicious"]: call_y = "SUSPICIOUS"

    elif sex == "XX":
        if y_ratio >= CFG["y_noise_threshold"]:
            call_y, type_y, final_y = "SUSPICIOUS", "SRY contamination / mismatch", min(y_ratio / 0.5, 1.0)
        else:
            call_y, final_y = "NORMAL", 0.0

    results.append(dict(
        chrom="chrY", log2fc=log2fc_y, robust_z=rz_y, expected_copy=expected_y,
        copy_score=copy_s_y, baf_score=0.0, ter_score=ter_s_y, cer_score=cer_s_y, pval_score=0.0,
        mosaic_score=0.0, final_score=final_y, call=call_y,
        detail=f"Y_Ratio={y_ratio:.2f} {type_y}", sex_note=note_y,
    ))

    return results

# ═══════════════════════════════════════════════════════════════
# 8. 시각화 (축약형)
# ═══════════════════════════════════════════════════════════════

def plot_chromosome_overview(df_f, summary, call_df, sex, sample_id, out_dir):
    all_chroms = sort_chroms(df_f["chrom"].unique().tolist())
    chrom_x    = {c: i for i, c in enumerate(all_chroms)}
    cn_col = "log2_chrom_norm" if "log2_chrom_norm" in df_f.columns else "raw_count"
    
    panels = [
        (cn_col,             "Copy Number Signal",    "steelblue",     (-2.0, 2.0)),
        ("hetero_like_rate", "Hetero-like Rate",       "mediumseagreen",(0, 1)),
        ("homo_like_rate",   "Homo-like Rate (LOH)",  "tomato",        (0, 1)),
        ("raw_bin_TER",      "Trans Error Rate",       "darkorange",    (0, None)),
        ("raw_bin_CER",      "Cis Phasing Rate",       "mediumpurple",  (0, None)),
    ]

    fig, axes = plt.subplots(len(panels), 1, figsize=(max(20, len(all_chroms)), 4.5 * len(panels)), sharex=True)
    fig.suptitle(f"Chromosome Overview — {sample_id}  [Sex: {sex}]", fontsize=14, fontweight="bold")

    for ax, (col, title, color, ylim) in zip(axes, panels):
        if col not in df_f.columns: continue
        for sc in ["chrX", "chrY"]:
            if sc in chrom_x: ax.axvspan(chrom_x[sc] - 0.5, chrom_x[sc] + 0.5, alpha=0.10, color="purple", zorder=0)

        for chrom, grp in df_f.groupby("chrom"):
            xi = chrom_x[chrom]
            ax.scatter(xi + np.random.uniform(-0.3, 0.3, len(grp)), grp[col], alpha=0.25, s=7, color="purple" if chrom in ("chrX", "chrY") else color)

        med_col = f"{col}_median"
        if med_col in summary.columns:
            xs = [chrom_x[c] for c in all_chroms if c in summary["chrom"].values]
            ys = [summary.loc[summary["chrom"] == c, med_col].values[0] for c in all_chroms if c in summary["chrom"].values]
            ax.plot(xs, ys, color="black", linewidth=1.5, alpha=0.7, zorder=5)
            ax.scatter(xs, ys, color="black", s=28, zorder=6)

        ax.axhline(0, color="gray", linewidth=0.8, linestyle="--", alpha=0.5)

        for _, row in call_df.iterrows():
            if row["chrom"] not in chrom_x: continue
            xi = chrom_x[row["chrom"]]
            if row["call"] == "ABNORMAL": ax.axvspan(xi - 0.5, xi + 0.5, alpha=0.18, color="red", zorder=0)
            elif row["call"] == "SUSPICIOUS": ax.axvspan(xi - 0.5, xi + 0.5, alpha=0.12, color="orange", zorder=0)

        ax.set_ylabel(title, fontsize=9, fontweight="bold")
        if ylim[1] is not None: ax.set_ylim(ylim)
        ax.grid(axis="y", linestyle="--", alpha=0.3)

    axes[-1].set_xticks(range(len(all_chroms)))
    axes[-1].set_xticklabels(all_chroms, rotation=45, ha="right", fontsize=9)
    plt.tight_layout()
    fig.savefig(os.path.join(out_dir, "01_chromosome_overview.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)

def plot_final_call(call_df, sex, sample_id, out_dir):
    chroms = sort_chroms(call_df["chrom"].tolist())
    cdf    = call_df.set_index("chrom").reindex(chroms).reset_index()

    fig, ax = plt.subplots(figsize=(max(16, len(chroms)), 7))
    bar_colors = [{"NORMAL":"#AB47BC","SUSPICIOUS":"#FF9800","ABNORMAL":"#F44336"}.get(r["call"], "gray") if r["chrom"] in ("chrX", "chrY") else CALL_COLORS.get(r["call"], "gray") for _, r in cdf.iterrows()]
    
    ax.bar(range(len(chroms)), cdf["final_score"].values, color=bar_colors, edgecolor="black", linewidth=0.6, width=0.75)
    ax.axhline(CFG["call_thresh_high"], color="red", linestyle="--", linewidth=1.5)
    ax.axhline(CFG["call_thresh_low"],  color="orange", linestyle="--", linewidth=1.2)
    ax.axhline(-CFG["call_thresh_high"], color="red", linestyle="--", linewidth=1.5)
    ax.axhline(-CFG["call_thresh_low"],  color="orange", linestyle="--", linewidth=1.2)
    ax.axhline(0, color="black", linewidth=0.8)

    for i, (_, row) in enumerate(cdf.iterrows()):
        if row["call"] != "NORMAL":
            ypos = row["final_score"]
            offset = 0.05 if ypos >= 0 else -0.12
            ax.text(i, ypos + offset, row["call"], ha="center", fontsize=8, fontweight="bold", color=CALL_COLORS[row["call"]])
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

# ═══════════════════════════════════════════════════════════════
# 9. 메인 실행
# ═══════════════════════════════════════════════════════════════

def run(tsv_path: str, out_dir: str = None, force_sex: str = None):
    tsv_path  = Path(tsv_path)
    sample_id = tsv_path.stem.replace(".normalized", "")
    if out_dir is None: out_dir = str(tsv_path.parent)
    os.makedirs(out_dir, exist_ok=True)

    print("=" * 70)
    print(f"[★] NIPT SCA Ratio Hard-Cutoff Pipeline (v2.2): {sample_id}")
    print("=" * 70)

    df   = pd.read_csv(tsv_path, sep="\t")
    df_f = apply_qc_filter(df)

    # [수정] 성별과 함께 x_ratio, y_ratio 명시적 반환
    sex, x_ratio, y_ratio, auto_med = estimate_sex(df_f)
    if force_sex:
        sex = force_sex.upper()
        print(f"  Sex: {sex} [FORCED by user] (auto-estimated: y_ratio={y_ratio:.3f})")
    else:
        print(f"  Estimated sex: {sex} (X_ratio = {x_ratio:.3f}, Y_ratio = {y_ratio:.3f})")

    summary = compute_chrom_summary(df_f)

    auto_df     = df_f[~df_f["chrom"].isin(["chrX", "chrY"])].copy()
    auto_chroms = sort_chroms(auto_df["chrom"].unique().tolist())

    cn_col  = "log2_chrom_norm" if "log2_chrom_norm" in df_f.columns else "raw_count"
    loo_cn  = compute_robust_loo(summary, f"{cn_col}_median")
    loo_ter = compute_robust_loo(summary, "raw_bin_TER_median")
    loo_cer = compute_robust_loo(summary, "raw_bin_CER_median")
    bin_pvals = compute_bin_distribution_test(df_f, metric=cn_col)
    baf_sig   = compute_baf_signal(summary)

    auto_records = make_autosome_calls(auto_chroms, loo_cn, bin_pvals, baf_sig, loo_ter, loo_cer)
    
    # [수정] x_ratio와 y_ratio를 성염색체 분석 엔진에 직접 전달하여 룰 강제
    sex_records  = analyze_sex_chromosomes(df_f, summary, sex, x_ratio, y_ratio, auto_med)

    call_df = pd.DataFrame(auto_records + sex_records)

    print(f"\n{'Chrom':<8} {'Log2FC':>8} {'RobustZ':>8} {'Final':>7} Call")
    print("-" * 75)
    for _, row in call_df.set_index("chrom").reindex(sort_chroms(call_df["chrom"].tolist())).reset_index().iterrows():
        flag  = "⚠ " if row["call"] == "ABNORMAL" else ("△ " if row["call"] == "SUSPICIOUS" else "  ")
        log2  = f"{row['log2fc']:+.3f}" if not np.isnan(row["log2fc"]) else "   N/A"
        rz    = f"{row['robust_z']:+.2f}" if not np.isnan(row["robust_z"]) else "  N/A"
        print(f"{row['chrom']:<8} {log2:>8} {rz:>8} {row['final_score']:>7.3f}  {flag}{row['call']:<12} {row.get('detail','')}")

    plot_chromosome_overview(df_f, summary, call_df, sex, sample_id, out_dir)
    plot_final_call(call_df, sex, sample_id, out_dir)

    print(f"\n[Done] Pipeline finished. Out: {out_dir}")
    return call_df

if __name__ == "__main__":

    from glob import glob 

    file_list = glob("/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Results/*/data/*.normalized.tsv")
    for TSV_PATH in file_list:
        print(os.path.basename(TSV_PATH))
        run(TSV_PATH)