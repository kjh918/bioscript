"""
chrom_calling.py
─────────────────
Bin 단위 데이터를 염색체 단위로 요약하고, 이상/의심/정상(ABNORMAL/SUSPICIOUS/NORMAL)을
판정합니다.

  - final_score 가중치를 하드코딩(0.65/0.15/0.10/0.10) 대신 config.yaml 의
    w_copy/w_baf/w_ter/w_cer 를 사용하도록 수정.
  - 성염색체(chrX/chrY) 하드컷오프 로직을 sex_calling.classify_sex_chromosome()
    로 위임 (clinical_diagnosis.py 와 동일한 단일 판정 로직 공유).
"""

import numpy as np
import pandas as pd

from utils import log, safe_log2fc, chrom_key
from rules import CFG
from calling_sex import classify_sex_chromosome


# ═══════════════════════════════════════════════════════════════
# 1. 염색체별 요약 통계 생성 (Robust 버전)
# ═══════════════════════════════════════════════════════════════
def compute_chrom_summary(df):
    """
    bin 단위 DataFrame → 염색체별 중앙값/평균 요약.
    
    [CRITICAL FIX]
    단순 중앙값을 구하면 노이즈(특히 0 수렴 Bin) 때문에 대푯값이 왜곡됩니다.
    MAD 기반으로 ±3 MAD 범위를 벗어나는 극단적 Bin을 제외하고 요약 통계를 산출합니다.
    """
    candidate_cols = [
        "log2_chrom_norm",
        "bin_BAF",
        "hetero_like_rate",
        "homo_like_rate",
        "imbalance_rate",
        "adj_bin_TER",
        "adj_bin_CER",
    ]
    base_cols = [c for c in candidate_cols if c in df.columns]

    def _robust_median(series):
        vals = series.dropna().values
        if len(vals) < 3:
            return np.nanmedian(vals) if len(vals) > 0 else np.nan
        med = np.nanmedian(vals)
        mad = np.nanmedian(np.abs(vals - med))
        if mad == 0: return med
        # ±3 MAD 필터 적용
        valid_vals = vals[(vals >= med - 3 * mad) & (vals <= med + 3 * mad)]
        return np.nanmedian(valid_vals) if len(valid_vals) > 0 else med
        
    def _robust_mean(series):
        vals = series.dropna().values
        if len(vals) < 3:
            return np.nanmean(vals) if len(vals) > 0 else np.nan
        med = np.nanmedian(vals)
        mad = np.nanmedian(np.abs(vals - med))
        if mad == 0: return np.nanmean(vals)
        valid_vals = vals[(vals >= med - 3 * mad) & (vals <= med + 3 * mad)]
        return np.nanmean(valid_vals) if len(valid_vals) > 0 else np.nanmean(vals)

    rows = []
    for chrom, grp in df.groupby("chrom"):
        row = {"chrom": chrom, "n_bins": len(grp)}
        for c in base_cols:
            # log2_chrom_norm 등 중요 지표는 필터링된 대푯값 사용
            if c == "log2_chrom_norm":
                row[f"{c}_median"] = _robust_median(grp[c])
                row[f"{c}_mean"]   = _robust_mean(grp[c])
            else:
                row[f"{c}_median"] = grp[c].median()
                row[f"{c}_mean"]   = grp[c].mean()
        rows.append(row)

    summary = pd.DataFrame(rows)
    summary["_k"] = summary["chrom"].apply(chrom_key)
    return summary.sort_values("_k").drop(columns=["_k"]).reset_index(drop=True)


# ═══════════════════════════════════════════════════════════════
# 2. BAF 시그널 및 모자이시즘 점수
# ═══════════════════════════════════════════════════════════════
def compute_baf_signal(summary):
    result = {}
    for _, row in summary.iterrows():
        c, scores, n = row["chrom"], 0.0, 0

        baf = row.get("bin_BAF_median", np.nan)
        if not np.isnan(baf) and c not in ("chrX", "chrY"):
            scores += min(abs(baf - 0.5) / 0.17, 1.0)
            n += 1

        homo = row.get("homo_like_rate_mean", np.nan)
        if not np.isnan(homo):
            scores += min(
                max(homo - CFG["homo_rate_threshold"], 0) / (1 - CFG["homo_rate_threshold"]),
                1.0)
            n += 1

        imbal = row.get("imbalance_rate_mean", np.nan)
        if not np.isnan(imbal) and c not in ("chrX", "chrY"):
            scores += min(imbal * 2, 1.0)
            n += 1

        result[c] = scores / n if n > 0 else 0.0
    return result


def compute_mosaic_score(baf_median, hetero_mean, homo_mean, ter_rz, cer_rz):
    if any(np.isnan(v) for v in [baf_median, hetero_mean, homo_mean]):
        return 0.0

    baf_comp    = np.clip(abs(baf_median - 0.5) / 0.17, 0.0, 1.0)
    hetero_comp = np.clip(hetero_mean / 0.50, 0.0, 1.0)
    homo_comp   = np.clip(1.0 - (homo_mean / CFG["homo_rate_threshold"]), 0.0, 1.0)

    support = (0.35 * baf_comp
               + 0.25 * hetero_comp
               + 0.20 * homo_comp)
    return float(np.clip(support, 0.0, 1.0))


# ═══════════════════════════════════════════════════════════════
# 3. 통합 진단 엔진
# ═══════════════════════════════════════════════════════════════
def analyze_all_chromosomes(summary, sex, x_cn, y_cn):
    """
    [핵심 설계]
    normalize_by_chrom_with_sex() 에서 이미
        log2_chrom_norm = log2(obs / ref_median) - expected_log2
    를 수행했으므로 정상 염색체는 0.0 근방에 분포합니다.
    따라서 여기서는 expected_log2를 **다시 빼지 않고**,
    상염색체 분포를 기준으로 순수 robust-Z만 계산합니다.
    """
    results     = []
    auto_mask   = ~summary["chrom"].isin(["chrX", "chrY"])
    auto_summary = summary[auto_mask]
    baf_sig      = compute_baf_signal(summary)

    # ── 상염색체 log2 분포의 MAD 추정 ────────────────────────────
    def auto_mad(col):
        vals = auto_summary[col].dropna().values if col in auto_summary.columns else np.array([])
        if len(vals) < 2:
            return 1e-9
        return max(np.median(np.abs(vals - np.median(vals))), 1e-9)

    # ── 특정 메트릭의 robust-Z ────────────────────────────────────
    def metric_rz(col, val):
        if col not in auto_summary.columns or np.isnan(val):
            return 0.0
        vals = auto_summary[col].dropna().values
        if len(vals) < 2:
            return 0.0
        med = np.median(vals)
        mad = max(np.median(np.abs(vals - med)), 1e-9)
        return float((val - med) / (1.4826 * mad))

    thresh = CFG.get("robust_z_threshold", 2.5)

    # final_score 가중치 (config.yaml 에서 로드; 과거 하드코딩값을 대체)
    w_copy = CFG.get("w_copy", 0.4)
    w_baf  = CFG.get("w_baf", 0.25)
    w_ter  = CFG.get("w_ter", 0.15)
    w_cer  = CFG.get("w_cer", 0.15)

    for _, row in summary.iterrows():
        chrom    = row["chrom"]
        log2_val = row.get("log2_chrom_norm_median", np.nan)

        # ── Robust-Z (상염색체 0.0 기준) ─────────────────────────
        # [CRITICAL FIX] 정상 데이터가 모여있는 상염색체 그룹의 분산(MAD)을 기준으로 평가
        if not np.isnan(log2_val):
            base_mad = auto_mad("log2_chrom_norm_median")
            # 만약 MAD가 너무 작아서 RZ가 폭발하는 것을 막기 위해 최소 분산 방어 추가
            base_mad = max(base_mad, 0.05) 
            rz = (log2_val - 0.0) / (1.4826 * base_mad)
        else:
            rz = np.nan

        # copy_s : [-1, 1] 클립 (|rz| = thresh 일 때 ±0.5, thresh*2 일 때 ±1.0)
        copy_s    = float(np.clip(rz / (thresh * 2), -1.0, 1.0)) if not np.isnan(rz) else 0.0
        direction = np.sign(copy_s) if copy_s != 0 else 1.0

        # ── BAF 시그널 ────────────────────────────────────────────
        baf_s = baf_sig.get(chrom, 0.0)

        # ── Homo / Hetero 패널티 ──────────────────────────────────
        homo_mean   = row.get("homo_like_rate_mean",   np.nan)
        hetero_mean = row.get("hetero_like_rate_mean", np.nan)

        homo_penalty   = 0.0
        hetero_penalty = 0.0
        if not np.isnan(homo_mean) and chrom not in ("chrX", "chrY"):
            homo_penalty   = 1.0 - homo_mean
        if not np.isnan(hetero_mean) and chrom not in ("chrX", "chrY"):
            hetero_penalty = hetero_mean

        # ── TER / CER robust-Z (모자이시즘 전용) ──────────────────
        ter_rz = metric_rz("adj_bin_TER_median", row.get("adj_bin_TER_median", np.nan))
        cer_rz = metric_rz("adj_bin_CER_median", row.get("adj_bin_CER_median", np.nan))

        # ── Final Score (config 가중치 사용) ───────────────────────
        final_score = (
            w_copy * copy_s
            + w_baf * baf_s          * direction
            + w_ter * homo_penalty   * (-1.0 if copy_s < 0 else 1.0)
            + w_cer * hetero_penalty * direction
        )

        # ── Mosaic Score ──────────────────────────────────────────
        baf_median = row.get("bin_BAF_median", np.nan)
        mosaic_s   = compute_mosaic_score(
            baf_median, hetero_mean, homo_mean, ter_rz, cer_rz)

        # ── Call ──────────────────────────────────────────────────
        call = ("ABNORMAL"    if abs(final_score) >= CFG.get("call_thresh_high", 0.70)
                else "SUSPICIOUS" if abs(final_score) >= CFG.get("call_thresh_low", 0.50)
                else "NORMAL")

        detail = ""
        if not np.isnan(log2_val):
            baf_str  = f"{baf_median:.2f}"  if not np.isnan(baf_median)  else "N/A"
            homo_str = f"{homo_mean:.2f}"   if not np.isnan(homo_mean)   else "N/A"
            detail   = f"Log2={log2_val:+.2f} | RZ={rz:+.2f} | BAF={baf_str} | Hm={homo_str}"

        # ── 성염색체 하드 컷오프 (sex_calling 단일 규칙 사용) ────────
        sex_hard_call = classify_sex_chromosome(chrom, x_cn, y_cn)
        if sex_hard_call is not None:
            call        = sex_hard_call["call"]
            detail      = sex_hard_call["pattern_name"]
            final_score = sex_hard_call["score"]

        results.append(dict(
            chrom         = chrom,
            log2fc        = log2_val,
            robust_z      = rz,
            expected_copy = "Exp(0.0)",
            copy_score    = copy_s,
            baf_score     = baf_s,
            ter_score     = ter_rz,
            cer_score     = cer_rz,
            pval_score    = 0.0,
            mosaic_score  = mosaic_s,
            final_score   = final_score,
            call          = call,
            detail        = detail,
            sex_note      = "",
        ))

    return pd.DataFrame(results)