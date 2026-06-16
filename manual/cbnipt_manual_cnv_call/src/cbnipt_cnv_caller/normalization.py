import numpy as np
import pandas as pd

from utils import log
from rules import CFG

def apply_qc_bin(df, min_depth=1, min_coverage=0.5):
    """
    [NEW] GC/Mappability 필터(is_filtered)를 통과하고,
    Coverage 조건(Depth, Breadth)을 만족하는 정상 Bin만 추출합니다.
    단, 성염색체(chrX, chrY)는 Coverage 필터에서 예외로 두어 삭제되지 않도록 보호합니다.
    """
    if df is None or df.empty: return pd.DataFrame()
    
    # 1. 뼈대 품질 조건 (GC 및 Mappability 정상)
    valid_quality = (df["is_filtered"] == False)
    
    # 2. 커버리지 조건 (최소 리드 뎁스 및 넓이 충족)
    valid_coverage = (df["raw_count"] >= min_depth) & (df["breadth_ratio"] >= min_coverage)
    
    # 3. 성염색체 보호 마스크
    is_sex_chrom = df["chrom"].isin(["chrX", "chrY"])
    
    # 4. 최종 마스크 결합
    # 해석: 퀄리티가 좋은 Bin 중에서, (커버리지가 충분하거나 OR 성염색체인 경우)만 통과
    final_mask = valid_quality & (valid_coverage | is_sex_chrom)
    
    df_filtered = df[final_mask].copy().reset_index(drop=True)
    
    dropped_count = len(df) - len(df_filtered)
    log(f"QC Bin Filter: Removed {dropped_count} low-quality bins. {len(df_filtered)} bins remain.")
    
    return df_filtered

def normalize_and_estimate_sex(df, value_col="raw_count"):
    """상염색체 기준 LOO 정규화 및 생물학적 성별 비율 계산"""
    autosome_df = df[~df["chrom"].isin(["chrX", "chrY"])]
    global_median = autosome_df[value_col].median()
    
    y_bins = df[df["chrom"] == "chrY"]
    y_median = y_bins[value_col].median() if not y_bins.empty else 0.0
    y_ratio = y_median / global_median if global_median > 0 else 0.0
    
    x_bins = df[df["chrom"] == "chrX"]
    x_median = x_bins[value_col].median() if not x_bins.empty else 0.0
    x_ratio = x_median / global_median if global_median > 0 else 0.0

    is_male = y_ratio > CFG["y_presence_threshold"]
    sex_tag = "XY" if is_male else "XX"

    norm_results = []
    for chrom, group in df.groupby("chrom", sort=False):
        ref_df = df[~df["chrom"].isin(["chrX", "chrY", chrom])]
        current_median = ref_df[value_col].median() if not ref_df.empty else global_median
        
        # [핵심 수정] 타겟 평탄화 로직 완전 제거! 
        # 관측된 생물학적 Log2 비율을 그대로 저장합니다. (남성 XY의 경우 알아서 -1.0으로 계산됨)
        observed_log2 = np.log2((group[value_col] + 1e-9) / (current_median + 1e-9))
        group["log2_chrom_norm"] = observed_log2
        
        norm_results.append(group)

    return pd.concat(norm_results).reset_index(drop=True), sex_tag, x_ratio, y_ratio

def normalize_by_chrom_with_sex(df, value_col="log2_corrected"):
    """
    [PERFECT FIX] Leave-One-Out (LOO) 정규화 적용.
    각 상염색체를 정규화할 때, 자기 자신을 제외한 나머지 상염색체들의 중앙값을 기준으로 삼아
    거대 염색체의 Aneuploidy가 전체 기저선(Baseline)을 왜곡하는 것을 방지합니다.
    """
    df_norm = df.copy()

    # 1. Log2 공간에 있는 보정값을 선형(Linear) 카운트 공간으로 복원
    if "log2" in value_col:
        linear_vals = 2 ** df_norm[value_col]
    else:
        linear_vals = df_norm[value_col]

    df_norm["loo_baseline"] = 1e-9 # 임시 기저선 컬럼
    auto_mask = ~df_norm["chrom"].isin(["chrX", "chrY"])
    
    # 성염색체 평가 및 안전장치용 전체 상염색체 중앙값
    global_auto_median = linear_vals[auto_mask].median()
    if global_auto_median == 0 or np.isnan(global_auto_median):
        global_auto_median = 1e-9

    # 2. Leave-One-Out (LOO) Baseline 동적 할당
    for chrom in df_norm["chrom"].unique():
        chrom_mask = df_norm["chrom"] == chrom
        
        if chrom not in ["chrX", "chrY"]:
            # 자기 자신을 제외한 나머지 상염색체(Rest of Autosomes)의 중앙값 추출
            loo_mask = auto_mask & (~chrom_mask)
            loo_median = linear_vals[loo_mask].median()
            
            if loo_median == 0 or np.isnan(loo_median):
                loo_median = global_auto_median
                
            df_norm.loc[chrom_mask, "loo_baseline"] = loo_median
        else:
            # 성염색체(chrX, chrY)는 전체 상염색체의 중앙값을 기준으로 평가
            df_norm.loc[chrom_mask, "loo_baseline"] = global_auto_median

    # 3. Y축 기준점을 2.0으로 하는 'Copy Number Signal' 생성 (LOO 적용 완료)
    df_norm["copy_number_signal"] = 2.0 * (linear_vals / df_norm["loo_baseline"])

    # 4. 하위 모듈(Segmentation) 호환성을 위해 log2_chrom_norm 도 생성
    with np.errstate(divide='ignore'):
        df_norm["log2_chrom_norm"] = np.log2(df_norm["copy_number_signal"] / 2.0)
    df_norm["log2_chrom_norm"] = df_norm["log2_chrom_norm"].replace([-np.inf], -3.0)

    # 임시 컬럼 삭제
    df_norm = df_norm.drop(columns=["loo_baseline"])

    # 5. 직관적인 성별 판별
    x_mask = df_norm["chrom"] == "chrX"
    y_mask = df_norm["chrom"] == "chrY"

    x_median = df_norm.loc[x_mask, "copy_number_signal"].median() if x_mask.sum() > 0 else 0
    y_median = df_norm.loc[y_mask, "copy_number_signal"].median() if y_mask.sum() > 0 else 0

    if x_median < 1.5 and y_median > 0.2:
        sex_tag = "XY"
    else:
        sex_tag = "XX"

    log(f"Normalization Mode: Leave-One-Out (LOO) | Sex: {sex_tag} (X Copy: {x_median:.2f}, Y Copy: {y_median:.2f})")

    return df_norm, sex_tag, float(x_median / 2.0), float(y_median / 2.0)