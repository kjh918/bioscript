"""
Single-Sample Chromosome Abnormality Detector (v2 - Sex Chromosome Support)
============================================================================
성염색체(chrX, chrY)를 포함한 전체 염색체 이상 탐지.

핵심 설계:
  [상염색체]  Robust LOO + Bin 분포 검정 + BAF + TER
  [성염색체]  성별 자동 추정 → 추정 성별 기준 독립 정규화 → 이상 판정

성별 추정 방법:
  chrY coverage ratio (chrY median / autosome median)
    ≥ 0.15  → XY (male)
    < 0.15  → XX (female)  [Y 없거나 노이즈 수준]

성염색체 정규화:
  XX female : chrX → 상염색체 중앙값 대비 Log2FC  (정상 ≈ 0)
              chrY → 없어야 정상 (coverage > threshold 이면 이상)
  XY male   : chrX → 상염색체 중앙값의 0.5배 대비 Log2FC  (정상 ≈ 0, 1 copy)
              chrY → 상염색체 중앙값의 0.5배 대비 Log2FC  (정상 ≈ 0, 1 copy)

이상 유형:
  XX : chrX↑ → XXX, chrX↓ → Turner(X0), chrY+ → SRY/contamination
  XY : chrX↑ → Klinefelter(XXY), chrX↓ → X0Y(?), chrY↑ → XYY, chrY↓ → 결실

출력:
  01_chromosome_overview.png    : 전체 염색체(성염색체 포함) 지표 라인 플롯
  02_bin_distribution.png       : 이상 염색체 bin-level violin 비교
  03_anomaly_score_heatmap.png  : 스코어 성분 히트맵
  04_final_call.png             : 최종 판정 바 차트
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
# 0. 설정
# ═══════════════════════════════════════════════════════════════

CFG = dict(
    # QC
    min_depth            = 1,
    min_coverage         = 0.5,

    # 성별 추정
    y_presence_threshold = 0.15,   # chrY/autosome median ratio 기준

    # 이상 탐지
    robust_z_thresh      = 2.5,
    pval_thresh          = 0.05,
    homo_rate_thresh     = 0.6,
    ter_z_thresh         = 2.0,

    # 성염색체 이상 판정 Log2FC 임계값
    sex_log2_abnormal    = 0.5,   # |log2fc| ≥ 0.50 → ABNORMAL
    sex_log2_suspicious  = 0.28,   # |log2fc| ≥ 0.28 → SUSPICIOUS
    y_noise_threshold    = 0.10,   # XX female에서 chrY ratio 이 이상이면 경고

    # 
    mosaic_log2_min = 0.10,
    mosaic_log2_max = 0.45,
    mosaic_bin_shift_min = 0.55,
    mosaic_score_thresh = 0.45,

    # 가중치
    w_copy  = 0.40,
    w_baf   = 0.30,
    w_ter   = 0.4,
    w_pval  = 0.15,

    # 최종 판정
    call_thresh_high = 0.50,     # 이 이상이면 Abnormal
    call_thresh_low  = 0.25,     # 이 이상이면 Suspicious
)

CALL_COLORS = {"NORMAL": "#4CAF50", "SUSPICIOUS": "#FF9800", "ABNORMAL": "#F44336"}


# ═══════════════════════════════════════════════════════════════
# 1. 유틸리티
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
    #return np.log2((val + 1e-9) / (ref + 1e-9))
    return float((val + 1e-9) / (ref + 1e-9))

def robust_z_single(val, others_arr):
    """단일 값 vs others 배열에 대해 MAD-기반 Robust Z 계산."""
    others = np.array([v for v in others_arr if not np.isnan(v)])
    if len(others) < 2:
        return np.nan, np.nan
    med = np.median(others)
    mad = np.median(np.abs(others - med))
    mad = max(mad, 1e-9)
    log2fc   = safe_log2fc(val, med)
    robust_z = (val - med) / (1.4826 * mad)
    return log2fc, robust_z


# ═══════════════════════════════════════════════════════════════
# 2. QC 필터 (성염색체 포함, 동일 기준 적용)
# ═══════════════════════════════════════════════════════════════

def apply_qc_filter(df):
    """성염색체에도 동일 QC를 적용한다."""
    df_f = df[
        (df["raw_count"]     >= CFG["min_depth"]) &
        (df["breadth_ratio"] >= CFG["min_coverage"])
    ].copy()
    df_f["_k"] = df_f["chrom"].apply(chrom_key)
    df_f = df_f.sort_values(["_k", "start"]).drop(columns=["_k"]).reset_index(drop=True)
    return df_f


# ═══════════════════════════════════════════════════════════════
# 3. 성별 추정
# ═══════════════════════════════════════════════════════════════

def estimate_sex(df_f):
    """
    chrY의 bin당 평균 coverage를 상염색체 중앙값과 비교해 성별 추정.
    반환: ("XX" | "XY", y_ratio, autosome_median)
    """
    auto_df = df_f[~df_f["chrom"].isin(["chrX", "chrY"])]
    y_df    = df_f[df_f["chrom"] == "chrY"]

    cn_col = "log2_chrom_norm" if "log2_chrom_norm" in df_f.columns else "raw_count"

    auto_med = auto_df[cn_col].median() if len(auto_df) > 0 else np.nan

    if len(y_df) == 0 or np.isnan(auto_med) or auto_med < 1e-6:
        return "XX", 0.0, auto_med

    y_med   = y_df[cn_col].median()
    y_ratio = y_med / auto_med

    sex = "XY" if y_ratio >= CFG["y_presence_threshold"] else "XX"
    return sex, y_ratio, auto_med


# ═══════════════════════════════════════════════════════════════
# 4. 염색체별 요약 통계
# ═══════════════════════════════════════════════════════════════

def compute_chrom_summary(df):
    base_cols = [
        "raw_count", "breadth_ratio", "hetero_like_rate",
        "homo_like_rate", "imbalance_rate", "raw_bin_TER",
        "raw_bin_CER", "bin_BAF", "MAD_BAF",
    ]
    extra_cols = ["log2_chrom_norm", "Informative_OE_Ratio"]
    all_cols = base_cols + [c for c in extra_cols if c in df.columns]

    rows = []
    for chrom, grp in df.groupby("chrom"):
        row = {"chrom": chrom, "n_bins": len(grp)}
        for col in all_cols:
            if col in grp.columns:
                row[f"{col}_median"] = grp[col].median()
                row[f"{col}_mean"]   = grp[col].mean()
                row[f"{col}_std"]    = grp[col].std()
        rows.append(row)

    summary = pd.DataFrame(rows)
    summary["_k"] = summary["chrom"].apply(chrom_key)
    summary = summary.sort_values("_k").drop(columns=["_k"]).reset_index(drop=True)
    return summary


# ═══════════════════════════════════════════════════════════════
# 5. 상염색체 신호 계산
# ═══════════════════════════════════════════════════════════════

def compute_robust_loo(summary, value_col):
    """상염색체 전용 LOO. 성염색체는 포함하지 않음."""
    if value_col not in summary.columns:
        return {}

    auto_summary = summary[~summary["chrom"].isin(["chrX", "chrY"])].copy()
    chroms = auto_summary["chrom"].tolist()
    vals   = auto_summary[value_col].values
    result = {}

    for i, c in enumerate(chroms):
        others = np.array([vals[j] for j in range(len(vals)) if j != i and not np.isnan(vals[j])])
        if len(others) < 3:
            result[c] = (np.nan, np.nan)
            continue
        log2fc, rz = robust_z_single(vals[i], others)
        result[c] = (log2fc, rz)

    return result


def compute_bin_distribution_test(df, metric):
    """각 상염색체 bin 분포 vs 나머지 상염색체 Mann-Whitney U 검정."""
    auto_df = df[~df["chrom"].isin(["chrX", "chrY"])]
    if metric not in auto_df.columns:
        return {}

    result = {}
    for chrom, grp in auto_df.groupby("chrom"):
        chrom_vals = grp[metric].dropna().values
        others     = auto_df[auto_df["chrom"] != chrom][metric].dropna().values
        if len(chrom_vals) < 3 or len(others) < 3:
            result[chrom] = (np.nan, np.nan, 0)
            continue
        stat, pval = stats.mannwhitneyu(chrom_vals, others, alternative="two-sided")
        direction  = 1 if chrom_vals.mean() > others.mean() else -1
        result[chrom] = (stat, pval, direction)
    return result


def compute_baf_signal(summary):
    """BAF 패턴 이상 스코어 (0=정상, 1=최대이상). 성염색체 포함."""
    result = {}
    for _, row in summary.iterrows():
        chrom = row["chrom"]
        scores, n = 0.0, 0

        baf = row.get("bin_BAF_median", np.nan)
        if not np.isnan(baf):
            # 성염색체(X)는 정상 BAF 기준이 달라 상염색체 기준만 적용
            if chrom not in ("chrX", "chrY"):
                dev = abs(baf - 0.5)
                scores += min(dev / 0.17, 1.0)
                n += 1

        homo = row.get("homo_like_rate_mean", np.nan)
        if not np.isnan(homo):
            scores += min(max(homo - CFG["homo_rate_thresh"], 0) / (1 - CFG["homo_rate_thresh"]), 1.0)
            n += 1

        imbal = row.get("imbalance_rate_mean", np.nan)
        if not np.isnan(imbal) and chrom not in ("chrX", "chrY"):
            scores += min(imbal * 2, 1.0)
            n += 1

        result[chrom] = scores / n if n > 0 else 0.0
    return result


# ═══════════════════════════════════════════════════════════════
# 6. 성염색체 독립 정규화 및 이상 판정
# ═══════════════════════════════════════════════════════════════

def analyze_sex_chromosomes(df_f, summary, sex, y_ratio, auto_med):
    """
    성별 기준으로 성염색체를 독립 정규화해 이상 판정.

    반환: list of dicts (chrom, log2fc, robust_z, expected_copy,
                         copy_score, baf_score, ter_score, pval_score,
                         final_score, call, detail, sex_note)
    """
    cn_col    = "log2_chrom_norm" if "log2_chrom_norm" in df_f.columns else "raw_count"
    cn_med    = f"{cn_col}_median"

    def get_median(chrom):
        row = summary[summary["chrom"] == chrom]
        if row.empty or cn_med not in row.columns:
            return np.nan
        return row.iloc[0][cn_med]

    def get_baf(chrom):
        row = summary[summary["chrom"] == chrom]
        if row.empty:
            return np.nan
        return row.iloc[0].get("homo_like_rate_mean", np.nan)

    def get_ter(chrom):
        row = summary[summary["chrom"] == chrom]
        if row.empty:
            return np.nan
        return row.iloc[0].get("raw_bin_TER_median", np.nan)

    results = []

    # ── 상염색체 중앙값 (정규화 기준) ──
    auto_summary = summary[~summary["chrom"].isin(["chrX", "chrY"])]
    auto_vals    = auto_summary[cn_med].dropna().values if cn_med in auto_summary.columns else np.array([])
    if len(auto_vals) == 0:
        return results

    auto_ref = np.median(auto_vals)  # 정상 diploid 2-copy 기준

    # ── chrX ──
    x_med = get_median("chrX")
    if not np.isnan(x_med):
        if sex == "XX":
            # XX: chrX는 diploid (2 copy) → 상염색체 중앙값과 동일해야 정상
            ref_x      = auto_ref
            expected_x = "2 copy (like autosomes)"
            note_x     = "XX female: chrX should equal autosome level"
        else:
            # XY: chrX는 1 copy → 상염색체 중앙값의 절반이 정상
            ref_x      = auto_ref * 0.5
            expected_x = "1 copy (~0.5× autosomes)"
            note_x     = "XY male: chrX should be ~0.5× autosome level"

        log2fc_x = safe_log2fc(x_med, ref_x)
        rz_x     = (x_med - ref_x) / (1.4826 * max(np.median(np.abs(auto_vals - np.median(auto_vals))), 1e-9))

        copy_s_x = float(np.clip(log2fc_x / 1.0, -1, 1))  # log2fc [-1,1] 범위 정규화
        baf_s_x  = min(max(get_baf("chrX") - CFG["homo_rate_thresh"], 0)
                       / (1 - CFG["homo_rate_thresh"]), 1.0) if not np.isnan(get_baf("chrX")) else 0.0
        ter_s_x  = float(np.clip((get_ter("chrX") - np.median([get_ter(c) for c in auto_summary["chrom"]
                                   if not np.isnan(get_ter(c))])) /
                                  (1.4826 * max(1e-9, np.median(np.abs(
                                      [get_ter(c) - np.median([get_ter(cc) for cc in auto_summary["chrom"]
                                                                if not np.isnan(get_ter(cc))])
                                       for c in auto_summary["chrom"] if not np.isnan(get_ter(c))])))) /
                                  (CFG["ter_z_thresh"] * 2), -1, 1)) if not np.isnan(get_ter("chrX")) else 0.0

        final_x = (CFG["w_copy"] * copy_s_x + CFG["w_baf"] * baf_s_x +
                   CFG["w_ter"] * ter_s_x)

        abs_log2 = abs(log2fc_x)
        if abs_log2 >= CFG["sex_log2_abnormal"]:
            call_x = "ABNORMAL"
        elif abs_log2 >= CFG["sex_log2_suspicious"]:
            call_x = "SUSPICIOUS"
        else:
            call_x = "NORMAL"

        # 이상 유형 추정
        type_x = ""
        if call_x != "NORMAL":
            if sex == "XX":
                type_x = "XXX?" if log2fc_x > 0 else "Turner(X0)?"
            else:
                type_x = "Klinefelter(XXY)?" if log2fc_x > 0 else "X0Y/partial del?"

        results.append(dict(
            chrom="chrX", log2fc=log2fc_x, robust_z=rz_x,
            expected_copy=expected_x,
            copy_score=copy_s_x, baf_score=baf_s_x,
            ter_score=ter_s_x, pval_score=0.0,
            final_score=final_x, call=call_x,
            detail=f"Log2FC={log2fc_x:+.3f} vs ref={ref_x:.4f}  {type_x}",
            sex_note=note_x,
        ))

    # ── chrY ──
    y_med = get_median("chrY")
    if sex == "XY":
        # XY: chrY는 1 copy → 상염색체의 절반이 정상
        ref_y      = auto_ref * 0.5
        expected_y = "1 copy (~0.5× autosomes)"
        note_y     = "XY male: chrY should be ~0.5× autosome level"

        if not np.isnan(y_med):
            log2fc_y = safe_log2fc(y_med, ref_y)
            rz_y     = (y_med - ref_y) / (1.4826 * max(np.median(np.abs(auto_vals - np.median(auto_vals))), 1e-9))
            copy_s_y = float(np.clip(log2fc_y / 1.0, -1, 1))
            final_y  = CFG["w_copy"] * copy_s_y

            abs_log2 = abs(log2fc_y)
            if abs_log2 >= CFG["sex_log2_abnormal"]:
                call_y = "ABNORMAL"
            elif abs_log2 >= CFG["sex_log2_suspicious"]:
                call_y = "SUSPICIOUS"
            else:
                call_y = "NORMAL"

            type_y = ""
            if call_y != "NORMAL":
                type_y = "XYY?" if log2fc_y > 0 else "chrY deletion?"

            results.append(dict(
                chrom="chrY", log2fc=log2fc_y, robust_z=rz_y,
                expected_copy=expected_y,
                copy_score=copy_s_y, baf_score=0.0,
                ter_score=0.0, pval_score=0.0,
                final_score=final_y, call=call_y,
                detail=f"Log2FC={log2fc_y:+.3f} vs ref={ref_y:.4f}  {type_y}",
                sex_note=note_y,
            ))
        else:
            # chrY 데이터 없음 (XY인데 Y 데이터 없으면 경고)
            results.append(dict(
                chrom="chrY", log2fc=np.nan, robust_z=np.nan,
                expected_copy=expected_y,
                copy_score=0.0, baf_score=0.0,
                ter_score=0.0, pval_score=0.0,
                final_score=0.0, call="SUSPICIOUS",
                detail="chrY data missing in XY sample",
                sex_note=note_y,
            ))

    else:  # XX female
        # XX: chrY는 없어야 정상
        note_y = "XX female: chrY should be absent (noise level)"
        if not np.isnan(y_med) and y_ratio >= CFG["y_noise_threshold"]:
            call_y  = "SUSPICIOUS"
            type_y  = "SRY contamination? or sex mismatch?"
            final_y = min(y_ratio / 0.5, 1.0)  # ratio 0.5이면 최대 이상
        else:
            call_y  = "NORMAL"
            type_y  = ""
            final_y = 0.0

        if "chrY" in summary["chrom"].values:
            log2fc_y = safe_log2fc(y_med, auto_ref) if not np.isnan(y_med) else np.nan
            results.append(dict(
                chrom="chrY", log2fc=log2fc_y, robust_z=np.nan,
                expected_copy="0 copy (absent)",
                copy_score=final_y, baf_score=0.0,
                ter_score=0.0, pval_score=0.0,
                final_score=final_y, call=call_y,
                detail=f"Y_ratio={y_ratio:.3f}  {type_y}",
                sex_note=note_y,
            ))

    return results


# ═══════════════════════════════════════════════════════════════
# 6.9 Mosaicism 확인
# ═══════════════════════════════════════════════════════════════

def compute_mosaic_score(log2fc, robust_z, baf_score, ter_score, pval_score):
    """
    Mosaicism 가능성 score.
    완전 aneuploidy보다는 약한 copy-number shift + 보조 신호가 있을 때 높게 부여.
    """

    if np.isnan(log2fc):
        return 0.0

    abs_l2 = abs(log2fc)

    # 1) copy-number shift가 너무 작으면 mosaic 가능성 낮음
    #    너무 크면 full aneuploidy 가능성이 더 큼
    if abs_l2 < CFG["mosaic_log2_min"]:
        copy_component = 0.0
    elif abs_l2 <= CFG["mosaic_log2_max"]:
        copy_component = (abs_l2 - CFG["mosaic_log2_min"]) / (
            CFG["mosaic_log2_max"] - CFG["mosaic_log2_min"]
        )
    else:
        copy_component = max(0.0, 1.0 - (abs_l2 - CFG["mosaic_log2_max"]) / 0.40)

    copy_component = float(np.clip(copy_component, 0, 1))

    # 2) robust Z는 약~중등도 편차일 때 mosaic 쪽으로 가중
    if np.isnan(robust_z):
        rz_component = 0.0
    else:
        arz = abs(robust_z)
        rz_component = min(arz / CFG["robust_z_thresh"], 1.0)

    # 3) BAF/TER/pval 보조 신호
    support_component = (
        0.45 * copy_component +
        0.25 * rz_component +
        0.20 * abs(baf_score) +
        0.05 * abs(ter_score) +
        0.05 * abs(pval_score)
    )

    return float(np.clip(support_component, 0, 1))

# ═══════════════════════════════════════════════════════════════
# 7. 상염색체 최종 판정
# ═══════════════════════════════════════════════════════════════

def make_autosome_calls(chroms, loo_cn, bin_pvals, baf_sig, loo_ter):
    records = []
    for c in chroms:
        log2fc, rz, copy_s = np.nan, np.nan, 0.0
        if c in loo_cn:
            log2fc, rz = loo_cn[c]
            copy_s = float(np.clip(rz / (CFG["robust_z_thresh"] * 2), -1, 1)) if not np.isnan(rz) else 0.0

        baf_s = baf_sig.get(c, 0.0)

        ter_s = 0.0
        if c in loo_ter:
            _, ter_rz = loo_ter[c]
            ter_s = float(np.clip(ter_rz / (CFG["ter_z_thresh"] * 2), -1, 1)) if not np.isnan(ter_rz) else 0.0
        print(c, ter_s)
        pval_s = 0.0
        if c in bin_pvals:
            _, pval, direction = bin_pvals[c]
            if not np.isnan(pval) and pval < CFG["pval_thresh"]:
                pval_s = float(np.clip(-np.log10(pval + 1e-10) / 5, 0, 1)) * direction

        final = (CFG["w_copy"] * copy_s + CFG["w_baf"] * baf_s +
                 CFG["w_ter"] * ter_s + CFG["w_pval"] * pval_s)

        abs_f = abs(final)
        call = ("ABNORMAL"  if abs_f >= CFG["call_thresh_high"] else
                "SUSPICIOUS" if abs_f >= CFG["call_thresh_low"]  else "NORMAL")

        records.append(dict(
            chrom=c, log2fc=log2fc, robust_z=rz,
            expected_copy="2 copy",
            copy_score=copy_s, baf_score=baf_s,
            ter_score=ter_s, pval_score=pval_s,
            final_score=final, call=call,
            detail=(f"Log2FC={log2fc:+.3f}  RZ={rz:+.2f}" if not np.isnan(log2fc) else ""),
            sex_note="",
        ))
    return records


# ═══════════════════════════════════════════════════════════════
# 8. 시각화
# ═══════════════════════════════════════════════════════════════

def plot_chromosome_overview(df_f, summary, call_df, sex, sample_id, out_dir):
    """Figure 1: 전체 염색체 bin-level scatter + 중앙값 라인 (성염색체 구분 음영)"""
    all_chroms = sort_chroms(df_f["chrom"].unique().tolist())
    chrom_x    = {c: i for i, c in enumerate(all_chroms)}

    cn_col = "log2_chrom_norm" if "log2_chrom_norm" in df_f.columns else "raw_count"
    panels = [
        (cn_col,            "Copy Number Signal",    "steelblue",     (-2.0, 2.0)),
        ("hetero_like_rate","Hetero-like Rate",       "mediumseagreen",(0, 1)),
        ("homo_like_rate",  "Homo-like Rate (LOH)",  "tomato",        (0, 1)),
        ("raw_bin_TER",     "Trans Error Rate",       "darkorange",    (0, None)),
    ]

    fig, axes = plt.subplots(len(panels), 1, figsize=(max(20, len(all_chroms)), 5 * len(panels)),
                              sharex=True)
    fig.suptitle(f"Chromosome Overview — {sample_id}  [Estimated sex: {sex}]",
                 fontsize=14, fontweight="bold")

    for ax, (col, title, color, ylim) in zip(axes, panels):
        if col not in df_f.columns:
            ax.set_title(f"{title} (not available)", fontsize=9)
            continue

        # 성염색체 구간 배경
        for sc in ["chrX", "chrY"]:
            if sc in chrom_x:
                ax.axvspan(chrom_x[sc] - 0.5, chrom_x[sc] + 0.5,
                           alpha=0.10, color="purple", zorder=0)

        # bin scatter
        for chrom, grp in df_f.groupby("chrom"):
            xi = chrom_x[chrom]
            jitter = np.random.uniform(-0.3, 0.3, len(grp))
            c_color = "purple" if chrom in ("chrX", "chrY") else color
            ax.scatter(xi + jitter, grp[col], alpha=0.25, s=7, color=c_color)

        # 중앙값 라인
        med_col = f"{col}_median"
        if med_col in summary.columns:
            xs = [chrom_x[c] for c in all_chroms if c in summary["chrom"].values]
            ys = [summary.loc[summary["chrom"] == c, med_col].values[0]
                  for c in all_chroms if c in summary["chrom"].values]
            ax.plot(xs, ys, color="black", linewidth=1.5, alpha=0.7, zorder=5)
            ax.scatter(xs, ys, color="black", s=28, zorder=6)

        ax.axhline(0, color="gray", linewidth=0.8, linestyle="--", alpha=0.5)

        # 이상 염색체 강조
        for _, row in call_df.iterrows():
            c = row["chrom"]
            if c not in chrom_x: continue
            xi = chrom_x[c]
            if row["call"] == "ABNORMAL":
                ax.axvspan(xi - 0.5, xi + 0.5, alpha=0.18, color="red", zorder=0)
            elif row["call"] == "SUSPICIOUS":
                ax.axvspan(xi - 0.5, xi + 0.5, alpha=0.12, color="orange", zorder=0)

        ax.set_ylabel(title, fontsize=9, fontweight="bold")
        if ylim[1]: ax.set_ylim(ylim)
        ax.grid(axis="y", linestyle="--", alpha=0.3)

    axes[-1].set_xticks(range(len(all_chroms)))
    axes[-1].set_xticklabels(all_chroms, rotation=45, ha="right", fontsize=9)

    # 범례
    legend_els = [
        mpatches.Patch(color="red",    alpha=0.4, label="ABNORMAL"),
        mpatches.Patch(color="orange", alpha=0.4, label="SUSPICIOUS"),
        mpatches.Patch(color="purple", alpha=0.4, label="Sex chromosome"),
    ]
    axes[0].legend(handles=legend_els, fontsize=8, loc="upper right")

    plt.tight_layout()
    out = os.path.join(out_dir, "01_chromosome_overview.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  [+] Saved: {out}")


def plot_bin_distribution(df_f, call_df, sample_id, out_dir):
    """Figure 2: 이상 염색체 bin-level violin (성염색체 포함)"""
    cn_col = "log2_chrom_norm" if "log2_chrom_norm" in df_f.columns else "raw_count"
    if cn_col not in df_f.columns:
        return

    abn_chroms = call_df[call_df["call"].isin(["ABNORMAL", "SUSPICIOUS"])]["chrom"].tolist()
    if not abn_chroms:
        print("  [i] No abnormal chromosomes — Figure 2 skipped.")
        return

    # 비교 기준: 성염색체는 상염색체 전체와 비교
    n   = len(abn_chroms)
    fig, axes = plt.subplots(1, n, figsize=(max(6, 5 * n), 6), squeeze=False)
    fig.suptitle(f"Bin Distribution — Abnormal/Suspicious ({sample_id})",
                 fontsize=12, fontweight="bold")

    auto_vals = df_f[~df_f["chrom"].isin(["chrX", "chrY"])][cn_col].dropna().values

    for ax, chrom in zip(axes[0], abn_chroms):
        chrom_vals = df_f[df_f["chrom"] == chrom][cn_col].dropna().values
        ref_vals   = auto_vals  # 항상 상염색체 전체와 비교

        parts = ax.violinplot([ref_vals, chrom_vals], positions=[0, 1],
                              showmedians=True, showextrema=True)
        colors_v = ["#90CAF9", "#CE93D8" if chrom in ("chrX","chrY") else "#EF9A9A"]
        for body, vc in zip(parts["bodies"], colors_v):
            body.set_facecolor(vc); body.set_alpha(0.7)

        ax.scatter([0]*len(ref_vals),   ref_vals,   alpha=0.1, s=5, color="#1565C0")
        ax.scatter([1]*len(chrom_vals), chrom_vals, alpha=0.4, s=8,
                   color="#6A1B9A" if chrom in ("chrX","chrY") else "#B71C1C")

        # 성염색체 기대값 기준선
        row = call_df[call_df["chrom"] == chrom].iloc[0]
        if chrom in ("chrX", "chrY"):
            ax.axhline(np.median(ref_vals) * 0.5, color="purple", linestyle=":",
                       linewidth=1.2, label="Expected 1-copy (×0.5)")
            ax.axhline(np.median(ref_vals), color="blue", linestyle=":",
                       linewidth=1.0, label="Expected 2-copy (autosome med)")
        else:
            ax.axhline(0.585, color="red",  linestyle=":", linewidth=1, label="Trisomy (+0.585)")
            ax.axhline(-1.0,  color="blue", linestyle=":", linewidth=1, label="Monosomy (-1.0)")
        ax.axhline(0, color="gray", linestyle="--", linewidth=0.7)

        call  = row["call"]
        color = CALL_COLORS.get(call, "black")
        note  = row.get("sex_note", "")

        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Autosomes\n(reference)", f"{chrom}\n({call})"],
                           fontsize=10, fontweight="bold")
        ax.set_title(f"{chrom}  [{row.get('expected_copy','')}]\n{note}",
                     fontsize=9, color=color, fontweight="bold")
        ax.set_ylabel(cn_col, fontsize=9)
        ax.legend(fontsize=7, loc="upper right")
        ax.grid(axis="y", linestyle="--", alpha=0.3)

    plt.tight_layout()
    out = os.path.join(out_dir, "02_bin_distribution.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  [+] Saved: {out}")


def plot_score_heatmap(call_df, sample_id, sex, out_dir):
    """Figure 3: 스코어 성분 히트맵 (성염색체 구분 표시)"""
    score_cols = ["copy_score", "baf_score", "ter_score", "pval_score", "final_score"]
    col_labels  = ["Copy\n(LOO/Norm)", "BAF\nPattern", "TER\n(Error)", "Distrib\n(pval)", "FINAL"]

    chroms = sort_chroms(call_df["chrom"].tolist())
    mat    = call_df.set_index("chrom").reindex(chroms)[score_cols].values.T

    fig, ax = plt.subplots(figsize=(max(14, len(chroms) * 0.75), 5))
    fig.suptitle(f"Anomaly Score Components — {sample_id}  [Sex: {sex}]",
                 fontsize=12, fontweight="bold")

    im = ax.imshow(mat, cmap="RdBu_r", aspect="auto", vmin=-1, vmax=1)

    for i in range(len(score_cols)):
        for j, c in enumerate(chroms):
            val = mat[i, j]
            tc  = "white" if abs(val) > 0.6 else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                    color=tc, fontsize=8, fontweight="bold")

            if i == len(score_cols) - 1:
                call = call_df[call_df["chrom"] == c]["call"].values[0]
                lw   = 3 if call == "ABNORMAL" else (2 if call == "SUSPICIOUS" else 0)
                ec   = "red" if call == "ABNORMAL" else "orange"
                if lw:
                    ax.add_patch(plt.Rectangle((j-0.5, i-0.5), 1, 1,
                                               linewidth=lw, edgecolor=ec, facecolor="none"))

    # 성염색체 컬럼 구분 선
    for j, c in enumerate(chroms):
        if c in ("chrX", "chrY"):
            ax.add_patch(plt.Rectangle((j-0.5, -0.5), 1, len(score_cols),
                                       linewidth=0, facecolor="purple", alpha=0.08, zorder=0))

    ax.set_xticks(range(len(chroms)))
    ax.set_yticks(range(len(score_cols)))
    ax.set_xticklabels(chroms, rotation=45, ha="right", fontsize=9, fontweight="bold")
    ax.set_yticklabels(col_labels, fontsize=9, fontweight="bold")

    plt.colorbar(im, ax=ax, pad=0.02, fraction=0.025).set_label(
        "Score  (negative=low, positive=high/abnormal)", fontsize=9)
    ax.set_xticks(np.arange(-0.5, len(chroms), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(score_cols), 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=1.5)
    ax.tick_params(which="minor", bottom=False, left=False)

    plt.tight_layout()
    out = os.path.join(out_dir, "03_anomaly_score_heatmap.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  [+] Saved: {out}")


def plot_final_call(call_df, sex, sample_id, out_dir):
    """Figure 4: 최종 판정 바 차트 (성염색체 구분 색상)"""
    chroms = sort_chroms(call_df["chrom"].tolist())
    cdf    = call_df.set_index("chrom").reindex(chroms).reset_index()

    fig, ax = plt.subplots(figsize=(max(16, len(chroms)), 7))

    bar_colors = []
    for _, row in cdf.iterrows():
        if row["chrom"] in ("chrX", "chrY"):
            bar_colors.append({"NORMAL":"#AB47BC","SUSPICIOUS":"#FF9800","ABNORMAL":"#F44336"}
                               .get(row["call"], "gray"))
        else:
            bar_colors.append(CALL_COLORS.get(row["call"], "gray"))

    ax.bar(range(len(chroms)), cdf["final_score"].values,
           color=bar_colors, edgecolor="black", linewidth=0.6, width=0.75)

    ax.axhline( CFG["call_thresh_high"], color="red",    linestyle="--", linewidth=1.5)
    ax.axhline( CFG["call_thresh_low"],  color="orange", linestyle="--", linewidth=1.2)
    ax.axhline(-CFG["call_thresh_high"], color="red",    linestyle="--", linewidth=1.5)
    ax.axhline(-CFG["call_thresh_low"],  color="orange", linestyle="--", linewidth=1.2)
    ax.axhline(0, color="black", linewidth=0.8)

    # 성염색체 구역 음영
    for j, c in enumerate(chroms):
        if c in ("chrX", "chrY"):
            ax.axvspan(j - 0.5, j + 0.5, alpha=0.08, color="purple", zorder=0)

    # 이상 레이블
    for i, (_, row) in enumerate(cdf.iterrows()):
        if row["call"] != "NORMAL":
            ypos   = row["final_score"]
            offset = 0.05 if ypos >= 0 else -0.12
            ax.text(i, ypos + offset, row["call"], ha="center",
                    fontsize=8, fontweight="bold", color=CALL_COLORS[row["call"]])
            detail = row.get("detail", "")
            if detail:
                ax.text(i, ypos + offset - 0.13, detail.split("  ")[0],
                        ha="center", fontsize=6.5, color="gray")

    ax.set_xticks(range(len(chroms)))
    ax.set_xticklabels(chroms, rotation=45, ha="right", fontsize=10, fontweight="bold")
    ax.set_ylabel("Final Anomaly Score", fontsize=11, fontweight="bold")
    ax.set_ylim(-1.15, 1.15)
    ax.set_title(f"Final Chromosome Call — {sample_id}  [Estimated sex: {sex}]",
                 fontsize=13, fontweight="bold")

    legend_els = [
        mpatches.Patch(color="#4CAF50", label="NORMAL (autosome)"),
        mpatches.Patch(color="#AB47BC", label="NORMAL (sex chrom)"),
        mpatches.Patch(color="#FF9800", label="SUSPICIOUS"),
        mpatches.Patch(color="#F44336", label="ABNORMAL"),
        plt.Line2D([0],[0], color="red",    linestyle="--", label=f"Abnormal ≥{CFG['call_thresh_high']}"),
        plt.Line2D([0],[0], color="orange", linestyle="--", label=f"Suspicious ≥{CFG['call_thresh_low']}"),
    ]
    ax.legend(handles=legend_els, fontsize=8, loc="upper right", ncol=2)
    ax.grid(axis="y", linestyle="--", alpha=0.35)
    plt.tight_layout()

    out = os.path.join(out_dir, "04_final_call.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  [+] Saved: {out}")


# ═══════════════════════════════════════════════════════════════
# 9. 메인
# ═══════════════════════════════════════════════════════════════

def run(tsv_path: str, out_dir: str = None, force_sex: str = None):
    tsv_path  = Path(tsv_path)
    sample_id = tsv_path.stem.replace(".normalized", "")
    if out_dir is None:
        out_dir = str(tsv_path.parent)
    os.makedirs(out_dir, exist_ok=True)

    print("=" * 65)
    print(f"[★] Single-Sample CNV Analysis (v2): {sample_id}")
    print("=" * 65)

    df   = pd.read_csv(tsv_path, sep="\t")
    df_f = apply_qc_filter(df)
    print(f"  Bins after QC: {len(df_f)}  (from {len(df)})")

    # 성별 추정
    sex, y_ratio, auto_med = estimate_sex(df_f)
    if force_sex:
        sex = force_sex.upper()
        print(f"  Sex: {sex} [FORCED by user]  (auto-estimated: y_ratio={y_ratio:.3f})")
    else:
        print(f"  Estimated sex: {sex}  (chrY/autosome ratio = {y_ratio:.3f})")

    # 요약 통계 (성염색체 포함)
    summary = compute_chrom_summary(df_f)

    # 상염색체 분석
    auto_df    = df_f[~df_f["chrom"].isin(["chrX", "chrY"])].copy()
    auto_chroms = sort_chroms(auto_df["chrom"].unique().tolist())

    cn_col  = "log2_chrom_norm" if "log2_chrom_norm" in df_f.columns else "raw_count"
    cn_med  = f"{cn_col}_median"
    loo_cn  = compute_robust_loo(summary, cn_med)
    loo_ter = compute_robust_loo(summary, "raw_bin_TER_median")

    bin_pvals = compute_bin_distribution_test(df_f, metric=cn_col)
    baf_sig   = compute_baf_signal(summary)

    auto_records = make_autosome_calls(auto_chroms, loo_cn, bin_pvals, baf_sig, loo_ter)

    # 성염색체 독립 분석
    sex_records = analyze_sex_chromosomes(df_f, summary, sex, y_ratio, auto_med)

    # 통합
    all_records = auto_records + sex_records
    call_df     = pd.DataFrame(all_records)

    # 콘솔 요약
    print(f"\n{'Chrom':<8} {'Log2FC':>8} {'RobustZ':>8} {'Final':>7}  Call         Detail")
    print("-" * 75)
    for _, row in call_df.set_index("chrom").reindex(sort_chroms(call_df["chrom"].tolist())).reset_index().iterrows():
        flag  = "⚠ " if row["call"] == "ABNORMAL" else ("△ " if row["call"] == "SUSPICIOUS" else "  ")
        log2  = f"{row['log2fc']:+.3f}" if not np.isnan(row["log2fc"]) else "   N/A"
        rz    = f"{row['robust_z']:+.2f}" if not np.isnan(row["robust_z"]) else "  N/A"
        mark  = " [SEX]" if row["chrom"] in ("chrX", "chrY") else ""
        print(f"{row['chrom']:<8} {log2:>8} {rz:>8} {row['final_score']:>7.3f}  "
              f"{flag}{row['call']:<12}{mark}  {row.get('detail','')}")

    abn  = call_df[call_df["call"] == "ABNORMAL"]["chrom"].tolist()
    susp = call_df[call_df["call"] == "SUSPICIOUS"]["chrom"].tolist()
    print(f"\n  ABNORMAL   : {abn  if abn  else 'none'}")
    print(f"  SUSPICIOUS : {susp if susp else 'none'}")

    # 시각화
    print("\n[★] Generating figures ...")
    plot_chromosome_overview(df_f, summary, call_df, sex, sample_id, out_dir)
    plot_bin_distribution(df_f, call_df, sample_id, out_dir)
    plot_score_heatmap(call_df, sample_id, sex, out_dir)
    plot_final_call(call_df, sex, sample_id, out_dir)

    print(f"\n[Done] Results saved to: {out_dir}")
    return call_df


# ═══════════════════════════════════════════════════════════════
# 10. CLI
# ═══════════════════════════════════════════════════════════════

if __name__ == "__main__":
    #TSV_PATH  = "/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Results/cb_data/cbNIPT_24_04_01_DS.normalized.tsv"
    TSV_PATH = "/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Results/cbNIPT_24_04_01_DS/data/cbNIPT_24_04_01_DS.normalized.tsv"
    
    OUT_DIR   = None
    FORCE_SEX = None  # "XX" 또는 "XY" 로 강제 지정 가능, None이면 자동 추정

    if len(sys.argv) > 1:
        parser = argparse.ArgumentParser(description="Single-sample CNV/Aneuploidy detector v2")
        parser.add_argument("tsv",          help="Path to .normalized.tsv")
        parser.add_argument("--out",        default=None, help="Output directory")
        parser.add_argument("--sex",        default=None, choices=["XX","XY"],
                            help="Force sex (XX/XY). Default: auto-estimate from chrY coverage.")
        args      = parser.parse_args()
        TSV_PATH  = args.tsv
        OUT_DIR   = args.out
        FORCE_SEX = args.sex

    run(TSV_PATH, OUT_DIR, FORCE_SEX)