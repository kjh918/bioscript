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
    base_cols = ["log2_chrom_norm", "hetero_like_rate", "homo_like_rate", "imbalance_rate", "bin_BAF"]
    
    # 2. [핵심 업데이트] TER과 CER을 Hetero/Homo 비율로 보정 (Adjusted Metrics)
    if "raw_bin_TER" in df.columns and "hetero_like_rate" in df.columns:
        # Trans 에러는 Hetero 구간에서만 발견되므로 hetero_like_rate로 나누어 밀도 보정
        # (+1e-4는 분모가 0이 되어 무한대가 되는 것을 방지하는 안전장치)
        df["adj_bin_TER"] = df["raw_bin_TER"] / (df["hetero_like_rate"] + 1e-4)
        base_cols.append("adj_bin_TER")
        
    if "raw_bin_CER" in df.columns and "homo_like_rate" in df.columns:
        # Cis(ALT)는 변이가 존재하는 모든 구간(Hetero+Homo)에서 발견되므로 총 변이율로 보정
        total_variant_rate = df["hetero_like_rate"] + df["homo_like_rate"]
        df["adj_bin_CER"] = df["raw_bin_CER"] / (total_variant_rate + 1e-4)
        base_cols.append("adj_bin_CER")

    cols = [c for c in base_cols if c in df.columns]
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

def compute_mosaic_score(baf_median, hetero_mean, homo_mean, ter_score, cer_score):
    """
    [MODIFIED] Copy Number(Log2FC)와 Z-Score를 완전히 배제하고, 
    대립유전자(Allele) 및 Phasing 지표의 변동성만을 기반으로 Mosaicism을 산출합니다.
    """
    if np.isnan(baf_median) or np.isnan(hetero_mean) or np.isnan(homo_mean):
        return 0.0
        
    # 1. BAF 편차 증가: 0.5 기준에서 멀어질수록(늘어날수록) Mosaic 점수 상승
    baf_dev = abs(baf_median - 0.5)
    baf_comp = np.clip(baf_dev / 0.17, 0.0, 1.0)
    
    # 2. Hetero 비율 증가: 높아질수록 Mosaic 점수 상승
    hetero_comp = np.clip(hetero_mean / 0.50, 0.0, 1.0)
    
    # 3. Homo 비율 감소: 낮아질수록 Mosaic 점수 상승 (반비례 역산)
    homo_comp = np.clip(1.0 - (homo_mean / CFG["homo_rate_threshold"]), 0.0, 1.0)
    
    # 4. 통합 연산 (가중치 재분배)
    support = (0.35 * baf_comp) + (0.25 * hetero_comp) + (0.20 * homo_comp) + (0.10 * abs(ter_score)) + (0.10 * abs(cer_score))
    
    return float(np.clip(support, 0.0, 1.0))

# ═══════════════════════════════════════════════════════════════
# 3. 통합 진단 엔진 (상염색체 + 성염색체 하드 컷오프)
# ═══════════════════════════════════════════════════════════════
# ═══════════════════════════════════════════════════════════════
# 3. 통합 진단 엔진 (Copy, BAF, Homo, Hetero 순수 생물학적 지표 주도형)
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
            if chrom == "chrX": expected_log2 = -1.0 
            elif chrom == "chrY": expected_log2 = -1.0 
        
        # ─────────────────────────────────────────────────────────
        # 1. Copy Number Score (물리적 증감, -1.0 ~ 1.0)
        # ─────────────────────────────────────────────────────────
        rz = (log2_val - expected_log2) / (1.4826 * get_auto_mad("log2_chrom_norm_median")) if not np.isnan(log2_val) else np.nan
        copy_s = float(np.clip(rz / (CFG.get("robust_z_threshold", 2.5) * 2), -1.0, 1.0)) if not np.isnan(rz) else 0.0
        direction = np.sign(copy_s) if copy_s != 0 else 1.0

        # ─────────────────────────────────────────────────────────
        # 2. BAF Score (대립유전자 균형 붕괴, 0.0 ~ 1.0)
        # ─────────────────────────────────────────────────────────
        baf_s = baf_sig.get(chrom, 0.0)

        # ─────────────────────────────────────────────────────────
        # 3. Homo-like Score (결실 및 LOH 증거, 0.0 ~ 1.0)
        # ─────────────────────────────────────────────────────────
        homo_mean = row.get("homo_like_rate_mean", np.nan)
        homo_penalty = 0.0
        
        if not np.isnan(homo_mean) and chrom not in ["chrX", "chrY"]:
            # Homo가 1.0일 때 페널티 0.0, 작아질수록 페널티 최대 1.0까지 증가
            homo_penalty = 1.0 - homo_mean

        # ─────────────────────────────────────────────────────────
        # 4. Hetero-like Score (0에서 멀어질수록 페널티 증가)
        # ─────────────────────────────────────────────────────────
        hetero_mean = row.get("hetero_like_rate_mean", np.nan)
        hetero_penalty = 0.0
        
        if not np.isnan(hetero_mean) and chrom not in ["chrX", "chrY"]:
            # Hetero가 0.0일 때 페널티 0.0, 커질수록 페널티 최대 1.0까지 증가
            hetero_penalty = hetero_mean

        # ─────────────────────────────────────────────────────────
        # [품질 지표 추출] TER, CER은 메인 스코어에서 제외, 모자이시즘에만 활용
        # ─────────────────────────────────────────────────────────
        def get_rz_score(col):
            val = row.get(col, np.nan)
            if np.isnan(val): return 0.0
            vals = auto_summary[col].dropna().values
            if len(vals) < 2: return 0.05  # 1e-9 대신 기본 노이즈 할당
        
            calculated_mad = np.median(np.abs(vals - np.median(vals)))
            return max(calculated_mad, 0.08)
        
        ter_s = get_rz_score("adj_bin_TER_median")
        cer_s = get_rz_score("adj_bin_CER_median")
        
        # ─────────────────────────────────────────────────────────
        # 5. 최종 Abnormality Score 병합
        # (Copy: 40%, BAF: 25%, Homo: 20%, Hetero: 15%)
        # ─────────────────────────────────────────────────────────
        final_score = (0.65 * copy_s) + \
                      (0.15 * baf_s * direction) + \
                      (0.10 * homo_penalty * (-1.0 if copy_s < 0 else 1.0)) + \
                      (0.10 * hetero_penalty * direction)
        
        # Mosaic Score 계산 (TER, CER 포함)
        baf_median = row.get("bin_BAF_median", np.nan)
        mosaic_s = compute_mosaic_score(baf_median, hetero_mean, homo_mean, ter_s, cer_s)
        
        call = "ABNORMAL" if abs(final_score) >= CFG.get("call_thresh_high", 0.50) else ("SUSPICIOUS" if abs(final_score) >= CFG.get("call_thresh_low", 0.25) else "NORMAL")
        
        detail = f"Log2={log2_val:+.2f} | BAF={baf_median:.2f} | Hm={homo_mean:.2f}" if not np.isnan(log2_val) else ""

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