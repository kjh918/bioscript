import numpy as np
import pandas as pd

from utils import log, safe_log2fc, chrom_key
from rules import CFG


# ═══════════════════════════════════════════════════════════════
# 1. 염색체별 요약 통계 생성
# ═══════════════════════════════════════════════════════════════
def compute_chrom_summary(df):
    """
    bin 단위 DataFrame → 염색체별 중앙값/평균 요약.
    존재하는 컬럼만 동적으로 집계합니다.

    [전제] df의 log2_chrom_norm 은 normalize_by_chrom_with_sex() 에서
           observed_log2 - expected_log2 로 이미 보정되어 있습니다.
           즉 정상 염색체는 0.0 근방, 이상 염색체만 이탈하므로
           robust-Z 기준이 명확합니다.
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

    rows = []
    for chrom, grp in df.groupby("chrom"):
        row = {"chrom": chrom, "n_bins": len(grp)}
        for c in base_cols:
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
               + 0.20 * homo_comp
               + 0.10 * abs(ter_rz)
               + 0.10 * abs(cer_rz))
    return float(np.clip(support, 0.0, 1.0))


# ═══════════════════════════════════════════════════════════════
# 3. 통합 진단 엔진
# ═══════════════════════════════════════════════════════════════
def analyze_all_chromosomes(summary, sex, x_ratio, y_ratio):
    """
    Parameters
    ----------
    summary  : compute_chrom_summary() 결과 DataFrame
    sex      : normalize_by_chrom_with_sex() 가 반환한 gender_tag
               ("Male (XY)" / "Female (XX)" / "XY" / "XX")
    x_ratio  : raw_count 기준 chrX_median / autosome_median
    y_ratio  : raw_count 기준 chrY_median / autosome_median

    [핵심 설계]
    normalize_by_chrom_with_sex() 에서 이미
        log2_chrom_norm = log2(obs / ref_median) - expected_log2
    를 수행했으므로 정상 염색체는 0.0 근방에 분포합니다.
    따라서 여기서는 expected_log2를 **다시 빼지 않고**,
    상염색체 분포를 기준으로 순수 robust-Z만 계산합니다.
    → robust-Z ≈ 0 : 정상,  |rz| > threshold : 이상
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

    for _, row in summary.iterrows():
        chrom    = row["chrom"]
        log2_val = row.get("log2_chrom_norm_median", np.nan)

        # ── Robust-Z (상염색체 0.0 기준) ─────────────────────────
        # normalize_by_chrom_with_sex 가 이미 expected_log2를 뺐으므로
        # 여기서는 단순히 상염색체 중앙값(≈0) 대비 이탈 정도만 봅니다.
        if not np.isnan(log2_val):
            rz = (log2_val - 0.0) / (1.4826 * auto_mad("log2_chrom_norm_median"))
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

        # ── Final Score ───────────────────────────────────────────
        final_score = (
            0.65 * copy_s
            + 0.15 * baf_s          * direction
            + 0.10 * homo_penalty   * (-1.0 if copy_s < 0 else 1.0)
            + 0.10 * hetero_penalty * direction
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

        # ── 성염색체 하드 컷오프 (raw_count 비율 기준) ────────────
        # x_ratio / y_ratio 는 pipeline 에서 raw_count 기준으로 계산되어 전달됩니다.
        if chrom == "chrX":
            if sex in ("XX", "Female (XX)", "Female"):
                if x_ratio <= CFG["sex_mono_x_ratio"]:
                    call, detail, final_score = "ABNORMAL", "Turner(X0)", -0.7
                elif x_ratio >= CFG["sex_xxx_ratio"]:
                    call, detail, final_score = "ABNORMAL", "XXX Syndrome", 0.7
            elif sex in ("XY", "Male (XY)", "Male"):
                if x_ratio >= CFG["sex_xxy_ratio"]:
                    call, detail, final_score = "ABNORMAL", "XXY Syndrome", 0.7
                elif x_ratio <= 0.25:
                    call, detail, final_score = "ABNORMAL", "chrX Loss", -0.7

        elif chrom == "chrY":
            if sex in ("XY", "Male (XY)", "Male"):
                if y_ratio >= CFG["sex_xyy_ratio"]:
                    call, detail, final_score = "ABNORMAL", "XYY Syndrome", 0.7
                elif y_ratio <= CFG["sex_mono_y_ratio"]:
                    call, detail, final_score = "ABNORMAL", "chrY Loss", -0.7
            elif sex in ("XX", "Female (XX)", "Female"):
                if y_ratio >= CFG["y_noise_threshold"]:
                    call, detail, final_score = "SUSPICIOUS", "SRY Contamination", 0.4

        results.append(dict(
            chrom         = chrom,
            log2fc        = log2_val,
            robust_z      = rz,
            expected_copy = "Exp(0.0)",   # normalize 단계에서 이미 보정됨
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