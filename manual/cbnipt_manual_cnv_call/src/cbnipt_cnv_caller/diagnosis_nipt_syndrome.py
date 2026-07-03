"""
clinical_diagnosis.py
──────────────────────
NIPT 마커 DB(marker_df) 기반 임상 진단 (기존 classification.py 의 임상 진단
섹션을 옮김). 성염색체 증후군(Turner/Klinefelter/Jacob/TripleX) 판정 수치는
더 이상 이 파일에 두지 않고 config.yaml CFG["SEX_ANEUPLOIDY"] 를 통해
sex_calling.classify_sex_chromosome() 이 단일 판정합니다.
(chrom_calling.analyze_all_chromosomes() 와 완전히 동일한 규칙을 공유)
"""

import numpy as np
import pandas as pd

from rules import CFG
from calling_sex import determine_fetal_sex, classify_sex_chromosome

# =====================================================================
# [1] 질환별 임상 진단 룰 (Rule Dictionary)
# =====================================================================
# 성염색체 증후군은 CFG["SEX_ANEUPLOIDY"][sex_rule_key] 로 위임한다.
NIPT_DIAGNOSIS_RULES = {
    "turner syndrome":      {"rule_type": "SEX_ANEUPLOIDY", "sex_rule_key": "turner"},
    "klinefelter syndrome": {"rule_type": "SEX_ANEUPLOIDY", "sex_rule_key": "klinefelter"},
    "jacob syndrome":       {"rule_type": "SEX_ANEUPLOIDY", "sex_rule_key": "jacob"},
    "triple x syndrome":    {"rule_type": "SEX_ANEUPLOIDY", "sex_rule_key": "triple_x"},

    # --- 성별 판별 전용 룰 (NIPT_SEX 마커 대응) ---
    "nipt_sex": {"rule_type": "FETAL_SEX"},

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

    # sex_calling 의 단일 성별 판별 로직 사용 (기존 이 함수 안에 하드코딩되어 있던
    # "global_x_cn < 1.5 and global_y_cn > 0.2" 판별식을 대체)
    inferred_sex = determine_fetal_sex(global_x_cn, global_y_cn)
    if inferred_sex == "Unknown":
        inferred_sex = "XX" if global_x_cn >= global_y_cn else "XY"

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

        # [C] 성염색체 질환 룰 (sex_calling 단일 판정 로직 위임)
        elif rule_type == "SEX_ANEUPLOIDY":
            if chrom in ["chrX", "chrY"]:
                sex_call = classify_sex_chromosome(chrom, inferred_sex, global_x_cn, global_y_cn)
                expected_key = rule.get("sex_rule_key")
                expected_pattern = CFG["SEX_ANEUPLOIDY"].get(expected_key, {}).get("pattern_name")

                if sex_call is not None and sex_call["pattern_name"] == expected_pattern:
                    diagnosis = "POSITIVE"
                    reason.append(f"{sex_call['pattern_name']} detected (X:{global_x_cn:.2f}, Y:{global_y_cn:.2f})")
                else:
                    # 마커 정의 범위엔 못 미치지만, 겹치는 CNV segment 상 부분적 증가/감소가
                    # 관측된 경우의 보조 신호 (기존 로직 그대로 유지)
                    cfg_rule = CFG["SEX_ANEUPLOIDY"].get(expected_key, {})
                    if has_amp and cfg_rule.get("x_min", 2.0) > 2.0 and chrom == "chrX":
                        diagnosis = "SUSPICIOUS"
                        reason.append("Partial X Gain")
                    elif has_del and cfg_rule.get("x_max", 2.0) < 2.0 and chrom == "chrX":
                        diagnosis = "SUSPICIOUS"
                        reason.append("Partial X Loss")

        # [D] 정상 성별 판별 룰 (NIPT_SEX 마커)
        elif rule_type == "FETAL_SEX":
            diagnosis = "XY (Male)" if inferred_sex == "XY" else "XX (Female)"
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