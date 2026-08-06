"""
clinical_diagnosis.py
──────────────────────
NIPT 마커 DB(marker_df) 기반 임상 진단 (기존 classification.py 의 임상 진단
섹션을 옮김). 성염색체 증후군(Turner/Klinefelter/Jacob/TripleX) 판정 수치는
더 이상 이 파일에 두지 않고 config.yaml CFG["SEX_ANEUPLOIDY"] 를 통해
sex_calling.classify_sex_chromosome() 이 단일 판정합니다.
(chrom_calling.analyze_all_chromosomes() 와 완전히 동일한 규칙을 공유)

AMP(증폭)/DEL(결실) 판정 threshold는 syndrome 단위로 config.yaml
CFG["CLINICAL_DIAGNOSIS_RULES"] 에서 가져옵니다. syndrome 이름이 명시적으로
등록되어 있지 않으면 NIPT_GROUP 단위 기본값("_default_amp"/"_default_del")으로
폴백합니다.
"""

import numpy as np
import pandas as pd

from rules import CFG
from sex_calling import determine_fetal_sex, classify_sex_chromosome

# =====================================================================
# [1] 질환별 임상 진단 룰 (Rule Dictionary)
# =====================================================================
# 성염색체 증후군 / 성별판별 룰만 여기 남겨둔다. (수치 threshold가 없는 rule_type)
# AMP/DEL 은 CFG["CLINICAL_DIAGNOSIS_RULES"] 에서 syndrome 단위로 조회한다.
NIPT_DIAGNOSIS_RULES = {
    "turner syndrome":      {"rule_type": "SEX_ANEUPLOIDY", "sex_rule_key": "turner"},
    "klinefelter syndrome": {"rule_type": "SEX_ANEUPLOIDY", "sex_rule_key": "klinefelter"},
    "jacob syndrome":       {"rule_type": "SEX_ANEUPLOIDY", "sex_rule_key": "jacob"},
    "triple x syndrome":    {"rule_type": "SEX_ANEUPLOIDY", "sex_rule_key": "triple_x"},

    # --- 성별 판별 전용 룰 (NIPT_SEX 마커 대응) ---
    "nipt_sex": {"rule_type": "FETAL_SEX"},
}

# marker_df 의 NIPT_GROUP 값 -> CFG["CLINICAL_DIAGNOSIS_RULES"] 폴백 키
NIPT_GROUP_FALLBACK_KEY = {
    "autosome abnormality": "_default_amp",
    "micro deletion": "_default_del",
}

# 계층 구조 (Depth) 매핑 - 출력 정렬 시 활용
FEATURE_RANK = {
    "TargetChromosome": 1,
    "PartialChromosome": 2,
    "PrimaryTargetRegion": 3,
    "CoreRegion": 4,
    "CoreGene": 5
}


def _resolve_clinical_rule(syn_key, grp_key):
    """
    syndrome(marker_df SYNDROME, lowercase) 우선 조회 -> 없으면 NIPT_GROUP 기본값 폴백.
    반환값에는 항상 'rule_type' 이 포함된다 (없으면 'UNKNOWN').
    """
    # 1) 성염색체 이수성 / 성별판별 룰 (수치 threshold 없음)
    special_rule = NIPT_DIAGNOSIS_RULES.get(syn_key)
    if special_rule is not None:
        return special_rule

    # 2) syndrome 단위 AMP/DEL 룰
    clinical_rules = CFG["CLINICAL_DIAGNOSIS_RULES"]
    rule = clinical_rules.get(syn_key)
    if rule is not None:
        return rule

    # 3) NIPT_GROUP 단위 기본값 폴백 (신규 syndrome이 marker DB에만 추가된 경우 대비)
    fallback_key = NIPT_GROUP_FALLBACK_KEY.get(grp_key)
    if fallback_key is not None:
        return clinical_rules.get(fallback_key, {})

    return {}


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

    # sex_calling 의 단일 성별 판별 로직 사용
    inferred_sex = determine_fetal_sex(global_x_cn, global_y_cn)
    if inferred_sex == "Unknown":
        inferred_sex = "XX" if global_x_cn >= global_y_cn else "XY"

    for idx, row in marker_df.iterrows():
        syndrome = str(row.get("SYNDROME", "Unknown")).strip()
        nipt_group = str(row.get("NIPT_GROUP", "Unknown")).strip()
        chrom = str(row.get("CHROMOSOME", "")).strip()
        start = int(row.get("GENOMIC_POS_START", 0))
        end = int(row.get("GENOMIC_POS_END", 0))
        feature_name = str(row.get("FEATURE_NAME", "")).strip()
        feature_type = str(row.get("FEATURE_TYPE", "TargetChromosome")).strip()

        feat_rank = FEATURE_RANK.get(feature_type, 99)

        # 1. 겹치는 영역의 Bin Data 수집
        ov_bins = bin_df[(bin_df["chrom"] == chrom) &
                         (bin_df["start"] < end) &
                         (bin_df["end"] > start)]

        n_bins = len(ov_bins)
        if n_bins > 0:
            median_cn = float(ov_bins["copy_number_signal"].median()) ## 
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

        # 3. Rule Dictionary 기반 진단 로직 적용 (syndrome 단위 threshold)
        diagnosis = "NEGATIVE"
        reason = []

        syn_key = syndrome.lower()
        grp_key = nipt_group.lower()
        rule = _resolve_clinical_rule(syn_key, grp_key)
        rule_type = rule.get("rule_type", "UNKNOWN")

        if rule_type == "AMP":
            if not np.isnan(median_cn) and median_cn >= rule["pos_min"]:
                diagnosis = "POSITIVE"
                reason.append(f"Gain (CN:{median_cn:.2f})")
                if has_amp:
                    reason.append("Segment AMP overlap")
            elif not np.isnan(median_cn) and median_cn >= rule["sus_min"]:
                if has_amp:
                    diagnosis = "POSITIVE"
                    reason.append(f"Gain (segment overlap + marginal CN:{median_cn:.2f})")
                else:
                    diagnosis = "SUSPICIOUS"
                    reason.append(f"Marginal Gain (CN:{median_cn:.2f})")

        elif rule_type == "DEL":
            if not np.isnan(median_cn) and median_cn <= rule["pos_max"]:
                diagnosis = "POSITIVE"
                reason.append(f"Loss (CN:{median_cn:.2f})")
                if has_del:
                    reason.append("Segment DEL overlap")
            elif not np.isnan(median_cn) and median_cn <= rule["sus_max"]:
                if has_del:
                    diagnosis = "POSITIVE"
                    reason.append(f"Loss (segment overlap + marginal CN:{median_cn:.2f})")
                else:
                    diagnosis = "SUSPICIOUS"
                    reason.append(f"Marginal Loss (CN:{median_cn:.2f})")
            # median_cn 이 sus_max 초과면 has_del 여부와 무관하게 NEGATIVE 유지

            #if not np.isnan(mean_homo) and mean_homo >= rule.get("loh_homo_min", 0.85):
            #    reason.append(f"High LOH(Homo:{mean_homo:.2f})")

        # [C] 성염색체 질환 룰 (sex_calling 단일 판정 로직 위임, 성별 태그로 게이팅하지 않음)
        elif rule_type == "SEX_ANEUPLOIDY":
            if chrom in ["chrX", "chrY"]:
                sex_call = classify_sex_chromosome(chrom, global_x_cn, global_y_cn)
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

        low_resolution = n_bins < CFG.get("min_reliable_overlap_bins", 3)

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
            "EVIDENCE": " & ".join(reason) if reason else "Normal",
            # [신규] bin 수가 적어 median_cn 신뢰도가 낮은 경우 표시.
            # CoreGene처럼 물리적으로 작은 영역(수십kb)은 겹치는 bin이 1~2개뿐일 수
            # 있어 이 값이 True 인 결과는 확인검사(ddPCR 등)를 권장합니다.
            "LOW_RESOLUTION_WARNING": low_resolution,
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