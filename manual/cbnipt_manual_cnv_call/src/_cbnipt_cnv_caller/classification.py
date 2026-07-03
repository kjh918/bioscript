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
               + 0.20 * homo_comp)
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

import numpy as np
import pandas as pd

# =====================================================================
# [1] 질환별 임상 진단 룰 (Rule Dictionary)
# =====================================================================
NIPT_DIAGNOSIS_RULES = {
    # --- 특수 질환 (성염색체 증후군) ---
    "turner syndrome": {
        "rule_type": "SEX_ANEUPLOIDY",
        "x_min": 0.5, "x_max": 1.4,  # X가 1개 (0.5 ~ 1.4)
        "y_max": 0.4,                # Y가 없음 (0.0 ~ 0.4)
        "pattern_name": "Turner(X0)",
        "description": "Turner syndrome is a sex chromosomal disorder in females characterized by the absence of one X chromosome, leading to a 45,X karyotype. It can result in short stature, ovarian dysfunction, and various physical abnormalities."
    },
    "klinefelter syndrome": {
        "rule_type": "SEX_ANEUPLOIDY",
        "x_min": 1.6, "x_max": 2.4,  # X가 2개 (1.6 ~ 2.4)
        "y_min": 0.4, "y_max": 1.5,  # Y가 1개 (0.4 ~ 1.5)
        "pattern_name": "Klinefelter(XXY)",
        "description": "Klinefelter syndrome is a sex chromosomal disorder in males caused by the presence of an extra X chromosome, resulting in a 47,XXY karyotype. It can lead to hypogonadism, reduced fertility, and various physical and cognitive features."
    },
    "jacob syndrome": {
        "rule_type": "SEX_ANEUPLOIDY",
        "x_min": 0.5, "x_max": 1.5,  # X가 1개 (0.5 ~ 1.5)
        "y_min": 1.6, "y_max": 2.5,  # Y가 2개 (1.6 ~ 2.5)
        "pattern_name": "Jacob(XYY)",
        "description": "Jacob syndrome, also known as XYY syndrome, is a sex chromosomal disorder in males characterized by the presence of an extra Y chromosome, resulting in a 47,XYY karyotype. It can lead to tall stature, learning difficulties, and behavioral issues."
    },
    "triple x syndrome": {
        "rule_type": "SEX_ANEUPLOIDY",
        "x_min": 2.5, "x_max": 3.5,  # X가 3개 (2.5 ~ 3.5)
        "y_max": 0.4,                # Y가 없음 (0.0 ~ 0.4)
        "pattern_name": "TripleX(XXX)",
        "description": "Triple X syndrome is a sex chromosomal disorder in females characterized by the presence of an extra X chromosome, resulting in a 47,XXX karyotype. It can lead to tall stature, learning difficulties, and behavioral issues."
    },
    
    # --- [NEW] 성별 판별 전용 룰 (NIPT_SEX 마커 대응) ---
    "nipt_sex": {
        "rule_type": "FETAL_SEX"
    },
    
    # --- 일반 그룹 기본 룰 (Autosome / Micro Deletion) ---
    "autosome abnormality": {
        "rule_type": "AMP",
        "pos_min": 2.5,
        "sus_min": 2.3
    },
    "micro deletion": {
        "rule_type": "DEL",
        "pos_max": 1.5,
        "sus_max": 1.65,
        "loh_homo_min": 0.85,   # LOH 판정 Homo 비율
        "loh_hetero_max": 0.1
    }
}

# 계층 구조 (Depth) 매핑 - 출력 정렬 시 활용
FEATURE_RANK = {
    "TargetChromosome": 1,
    "PartialChromosome": 2,
    "PrimaryTargetRegion": 3,
    "CoreRegion": 4,
    "CoreGene": 5
}

# =====================================================================
# [2] 임상 진단 핵심 함수
# =====================================================================
def diagnose_clinical_markers(marker_df, bin_df, cnv_df, sex_tag="Unknown"):
    """
    NIPT Rule Dictionary 기반 임상 진단 모듈 (성별 판별 로직 통합)
    """
    results = []
    
    # 성염색체 증후군 다중 평가를 위해 환자의 전체 X, Y 중앙값을 사전 추출
    chrX_vals = bin_df[bin_df["chrom"] == "chrX"]["copy_number_signal"].dropna()
    chrY_vals = bin_df[bin_df["chrom"] == "chrY"]["copy_number_signal"].dropna()
    
    global_x_cn = float(chrX_vals.median()) if len(chrX_vals) > 0 else 2.0
    global_y_cn = float(chrY_vals.median()) if len(chrY_vals) > 0 else 0.0
    
    # -------------------------------------------------------------
    # [핵심] 연구원님의 직관적인 성별 판별 로직을 진단 모듈 베이스라인으로 편입
    # -------------------------------------------------------------
    if global_x_cn < 1.5 and global_y_cn > 0.2:
        inferred_sex = "XY"
    else:
        inferred_sex = "XX"
    
    for idx, row in marker_df.iterrows():
        syndrome = str(row.get("SYNDROME", "Unknown")).strip()
        nipt_group = str(row.get("NIPT_GROUP", "Unknown")).strip()
        chrom = str(row.get("CHROMOSOME", "")).strip()
        start = int(row.get("GENOMIC_POS_START", 0))
        end = int(row.get("GENOMIC_POS_END", 0))
        feature_name = str(row.get("FEATUER_NAME", "")).strip()
        feature_type = str(row.get("FEATURE_TYPE", "TargetChromosome")).strip()
        
        feat_rank = FEATURE_RANK.get(feature_type, 99)

        # 1. 겹치는 영역의 Bin Data 수집
        ov_bins = bin_df[(bin_df["chrom"] == chrom) & 
                         (bin_df["start"] < end) & 
                         (bin_df["end"] > start)]
        
        n_bins = len(ov_bins)
        if n_bins > 0:
            median_cn = float(ov_bins["copy_number_signal"].median())
            mean_hetero = float(ov_bins["hetero_like_rate"].mean())
            mean_homo = float(ov_bins["homo_like_rate"].mean())
            median_baf = float(ov_bins["bin_BAF"].median())
        else:
            median_cn, mean_hetero, mean_homo, median_baf = np.nan, np.nan, np.nan, np.nan

        # 2. 겹치는 CNV Segment 확인
        ov_cnvs = pd.DataFrame()
        if cnv_df is not None and not cnv_df.empty:
            ov_cnvs = cnv_df[(cnv_df["chrom"] == chrom) & 
                             (cnv_df["start"] < end) & 
                             (cnv_df["end"] > start)]
            
        seg_calls = ov_cnvs["cnv_call"].unique().tolist() if not ov_cnvs.empty else []
        has_del = "DEL" in seg_calls
        has_amp = "AMP" in seg_calls

        # 3. Rule Dictionary 기반 진단 로직 적용
        diagnosis = "NEGATIVE"
        reason = []
        
        syn_key = syndrome.lower()
        grp_key = nipt_group.lower()
        rule = NIPT_DIAGNOSIS_RULES.get(syn_key, NIPT_DIAGNOSIS_RULES.get(grp_key, {}))
        rule_type = rule.get("rule_type", "UNKNOWN")

        # [A] 상염색체 증폭 룰
        if rule_type == "AMP":
            if has_amp or (not np.isnan(median_cn) and median_cn >= rule["pos_min"]):
                diagnosis = "POSITIVE"
                reason.append(f"Gain (CN:{median_cn:.2f})")
            elif not np.isnan(median_cn) and median_cn >= rule["sus_min"]:
                diagnosis = "SUSPICIOUS"
                reason.append(f"Marginal Gain (CN:{median_cn:.2f})")

        # [B] 상염색체 미세결실 룰
        elif rule_type == "DEL":
            if has_del or (not np.isnan(median_cn) and median_cn <= rule["pos_max"]):
                diagnosis = "POSITIVE"
                reason.append(f"Loss (CN:{median_cn:.2f})")
            elif not np.isnan(median_cn) and median_cn <= rule["sus_max"]:
                diagnosis = "SUSPICIOUS"
                reason.append(f"Marginal Loss (CN:{median_cn:.2f})")
                
            if not np.isnan(mean_homo) and mean_homo >= rule.get("loh_homo_min", 0.85):
                reason.append(f"High LOH(Homo:{mean_homo:.2f})")

        # [C] 성염색체 질환 룰
        elif rule_type == "SEX_ANEUPLOIDY":
            if chrom in ["chrX", "chrY"]: 
                cond_x = True
                if "x_min" in rule: cond_x = cond_x and (global_x_cn >= rule["x_min"])
                if "x_max" in rule: cond_x = cond_x and (global_x_cn <= rule["x_max"])
                
                cond_y = True
                if "y_min" in rule: cond_y = cond_y and (global_y_cn >= rule["y_min"])
                if "y_max" in rule: cond_y = cond_y and (global_y_cn <= rule["y_max"])
                
                if cond_x and cond_y:
                    diagnosis = "POSITIVE"
                    reason.append(f"{rule['pattern_name']} detected (X:{global_x_cn:.2f}, Y:{global_y_cn:.2f})")
                else:
                    if has_amp and rule.get("x_min", 2.0) > 2.0 and chrom == "chrX":
                        diagnosis = "SUSPICIOUS"
                        reason.append(f"Partial X Gain")
                    elif has_del and rule.get("x_max", 2.0) < 2.0 and chrom == "chrX":
                        diagnosis = "SUSPICIOUS"
                        reason.append(f"Partial X Loss")

        # [D] 정상 성별 판별 룰 (NIPT_SEX 마커)
        elif rule_type == "FETAL_SEX":
            # POSITIVE/NEGATIVE 대신 직접적인 성별 판별 텍스트를 진단명으로 사용
            if inferred_sex == "XY":
                diagnosis = "XY (Male)"
            else:
                diagnosis = "XX (Female)"
            reason.append(f"Determined by Global CN (X:{global_x_cn:.2f}, Y:{global_y_cn:.2f})")

        # 4. 행 데이터 기록
        results.append({
            "SYNDROME": syndrome,
            "NIPT_GROUP": nipt_group,
            "FEAT_RANK": feat_rank,
            "FEATURE_TYPE": feature_type,
            "FEATURE_NAME": feature_name,
            "CHROMOSOME": chrom,
            "OVERLAP_BINS": n_bins,
            "OBSERVED_CN": median_cn,
            "BAF_MEDIAN": median_baf,
            "HETERO_RATE": mean_hetero,
            "HOMO_RATE": mean_homo,
            "DETECTED_CNV": ",".join(seg_calls) if seg_calls else "None",
            "DIAGNOSIS": diagnosis,
            "EVIDENCE": " & ".join(reason) if reason else "Normal"
        })

    df_results = pd.DataFrame(results)
    if not df_results.empty:
        df_results = df_results.sort_values(by=["SYNDROME", "FEAT_RANK"])
        
    return df_results

# =====================================================================
# [3] 요약(Summary) 헬퍼 함수
# =====================================================================
def get_diagnosis_summary(report_df):
    summary = []
    
    for syndrome, group in report_df.groupby("SYNDROME"):
        pos_features = group[group["DIAGNOSIS"] == "POSITIVE"]
        sus_features = group[group["DIAGNOSIS"] == "SUSPICIOUS"]
        
        # 성별 판독 행 처리 (NIPT_SEX)
        if syndrome == "NIPT_SEX":
            final_call = group["DIAGNOSIS"].iloc[0]
            summary.append({
                "SYNDROME": syndrome,
                "FINAL_CALL": final_call,
                "CRITICAL_HITS": "N/A",
                "DETAILS": group["EVIDENCE"].iloc[0]
            })
            continue
            
        final_call = "NEGATIVE"
        if len(pos_features) > 0:
            final_call = "POSITIVE"
        elif len(sus_features) > 0:
            final_call = "SUSPICIOUS"
            
        critical_hits = pos_features[pos_features["FEAT_RANK"] >= 4]["FEATURE_NAME"].tolist()
        
        summary.append({
            "SYNDROME": syndrome,
            "FINAL_CALL": final_call,
            "CRITICAL_HITS": ",".join(critical_hits) if critical_hits else "None",
            "DETAILS": f"Positive in {len(pos_features)} features, Suspicious in {len(sus_features)}"
        })
        
    return pd.DataFrame(summary)