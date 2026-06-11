import numpy as np
import pandas as pd
from .utils import log, safe_log2fc, chrom_key
from .rules import CFG


def apply_qc_filter(df):
    # [수정] 상염색체에는 QC 기준을 엄격히 적용하되, 성염색체(chrX, chrY)는 필터링에서 제외하여 보존
    qc_mask = (
        ((df["raw_count"] >= int(CFG["min_depth"])) & (df["breadth_ratio"] > float(CFG["min_coverage"]))) | 
        (df["chrom"].isin(["chrX", "chrY"]))
    )
    
    df_f = df[qc_mask].copy()
    df_f["_k"] = df_f["chrom"].apply(chrom_key)
    
    return df_f.sort_values(["_k", "start"]).drop(columns=["_k"]).reset_index(drop=True)

def compute_chrom_summary(df):
    cols = ["log2_chrom_norm", "hetero_like_rate", "homo_like_rate", "imbalance_rate", "raw_bin_TER", "raw_bin_CER", "bin_BAF"]
    cols = [c for c in cols if c in df.columns]
    rows = []
    for chrom, grp in df.groupby("chrom"):
        row = {"chrom": chrom, "n_bins": len(grp)}
        for c in cols: row[f"{c}_median"], row[f"{c}_mean"] = grp[c].median(), grp[c].mean()
        rows.append(row)
    summary = pd.DataFrame(rows)
    summary["_k"] = summary["chrom"].apply(chrom_key)
    return summary.sort_values("_k").drop(columns=["_k"]).reset_index(drop=True)

def compute_baf_signal(summary):
    result = {}
    for _, row in summary.iterrows():
        c, scores, n = row["chrom"], 0.0, 0
        baf = row.get("bin_BAF_median", np.nan)
        if not np.isnan(baf) and c not in ("chrX", "chrY"): scores += min(abs(baf - 0.5) / 0.17, 1.0); n += 1
        homo = row.get("homo_like_rate_mean", np.nan)
        if not np.isnan(homo): scores += min(max(homo - CFG["homo_rate_threshold"], 0) / (1 - CFG["homo_rate_threshold"]), 1.0); n += 1
        imbal = row.get("imbalance_rate_mean", np.nan)
        if not np.isnan(imbal) and c not in ("chrX", "chrY"): scores += min(imbal * 2, 1.0); n += 1
        result[c] = scores / n if n > 0 else 0.0
    return result

def compute_mosaic_score(log2fc, robust_z, baf_score, ter_score, cer_score):
    if np.isnan(log2fc): return 0.0
    abs_l2 = abs(log2fc)
    c_comp = 0.0 if abs_l2 < CFG["mosaic_log2_min"] else (
        (abs_l2 - CFG["mosaic_log2_min"]) / (CFG["mosaic_log2_max"] - CFG["mosaic_log2_min"]) if abs_l2 <= CFG["mosaic_log2_max"] 
        else max(0.0, 1.0 - (abs_l2 - CFG["mosaic_log2_max"]) / 0.40))
    rz_comp = 0.0 if np.isnan(robust_z) else min(abs(robust_z) / CFG["robust_z_hreshold"], 1.0)
    support = 0.40*c_comp + 0.20*rz_comp + 0.20*abs(baf_score) + 0.10*abs(ter_score) + 0.10*abs(cer_score)
    return float(np.clip(support, 0, 1))

# ═══════════════════════════════════════════════════════════════
# 3. 통합 진단 엔진 (상염색체 + 성염색체 하드 컷오프)
# ═══════════════════════════════════════════════════════════════
def analyze_all_chromosomes(summary, sex, x_ratio, y_ratio):
    results = []
    auto_summary = summary[~summary["chrom"].isin(["chrX", "chrY"])]
    baf_sig = compute_baf_signal(summary)
    
    def get_auto_mad(col):
        vals = auto_summary[col].dropna().values
        if len(vals) < 2: return 1e-9
        return max(np.median(np.abs(vals - np.median(vals))), 1e-9)

    for _, row in summary.iterrows():
        chrom = row["chrom"]
        log2_val = row.get("log2_chrom_norm_median", np.nan)
        
        expected_log2 = 0.0 # 상염색체 및 여성 XX는 0.0 (2-copy)
        if sex == "XY":
            if chrom == "chrX": expected_log2 = -1.0 # 남성 X는 -1.0 (1-copy)
            elif chrom == "chrY": expected_log2 = -1.0 # 남성 Y는 -1.0 (1-copy)
        
        # Z-score 계산 시: (내 값 - 정상기대치) / 상염색체_MAD
        rz = (log2_val - expected_log2) / (1.4826 * get_auto_mad("log2_chrom_norm_median")) if not np.isnan(log2_val) else np.nan
        copy_s = float(np.clip(rz / (CFG["robust_z_hreshold"] * 2), -1, 1)) if not np.isnan(rz) else 0.0
        
        baf_s = baf_sig.get(chrom, 0.0)
        
        def get_rz_score(col):
            val = row.get(col, np.nan)
            if np.isnan(val): return 0.0
            vals = auto_summary[col].dropna().values
            return float(np.clip(((val - np.median(vals)) / (1.4826 * get_auto_mad(col))) / (CFG["ter_z_thresh"] * 2), -1, 1))
        
        ter_s = get_rz_score("raw_bin_TER_median")
        cer_s = get_rz_score("raw_bin_CER_median")
        
        final_score = CFG["w_copy"]*copy_s + CFG["w_baf"]*baf_s + CFG["w_ter"]*ter_s + CFG["w_cer"]*cer_s
        # Mosaic score 계산 시에도 기대치를 뺀 순수 편차(abs_diff)를 넘겨주도록 수정
        mosaic_s = compute_mosaic_score(log2_val - expected_log2, rz, baf_s, ter_s, cer_s)
        
        call = "ABNORMAL" if abs(final_score) >= CFG["call_thresh_high"] else ("SUSPICIOUS" if abs(final_score) >= CFG["call_thresh_low"] else "NORMAL")
        detail = f"Log2={log2_val:+.3f}" if not np.isnan(log2_val) else ""

        # [하드 컷오프 오버라이드 로직 유지]
        if chrom == "chrX":
            if sex == "XX":
                if x_ratio <= CFG["sex_mono_x_ratio"]: call, detail, final_score = "ABNORMAL", "Turner(X0)", -0.7
                elif x_ratio >= CFG["sex_xxx_ratio"]: call, detail, final_score = "ABNORMAL", "XXX Syndrome", 0.7
            elif sex == "XY":
                if x_ratio >= CFG["sex_xxy_ratio"]: call, detail, final_score = "ABNORMAL", "XXY Syndrome", 0.7
                elif x_ratio <= 0.25: call, detail, final_score = "ABNORMAL", "chrX Loss", -0.7
        elif chrom == "chrY":
            if sex == "XY":
                if y_ratio >= CFG["sex_xyy_ratio"]: call, detail, final_score = "ABNORMAL", "XYY Syndrome", 0.7
                elif y_ratio <= CFG["sex_mono_y_ratio"]: call, detail, final_score = "ABNORMAL", "chrY Loss", -0.7
            elif sex == "XX":
                if y_ratio >= CFG["y_noise_threshold"]: call, detail, final_score = "SUSPICIOUS", "SRY Contamination", 0.4

        results.append(dict(
            chrom=chrom, log2fc=log2_val, robust_z=rz, expected_copy=f"Exp({expected_log2:.1f})",
            copy_score=copy_s, baf_score=baf_s, ter_score=ter_s, cer_score=cer_s, pval_score=0.0,
            mosaic_score=mosaic_s, final_score=final_score, call=call, detail=detail, sex_note=""
        ))
    return pd.DataFrame(results)