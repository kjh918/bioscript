import numpy as np
import pandas as pd
from utils import log 
import pandas as pd
import numpy as np
from utils import log

def normalize_by_chrom_with_sex(df, value_col, sex_threshold=0.001):
    """
    [MODIFIED] 
    1. Default 인자 제거: value_col, sex_threshold 미지정 시 에러 발생.
    2. Log2 Space 적용: 2n=0.0, 1n=-1.0 기준으로 타겟 설정.
    3. Fail-fast: 필수 컬럼 및 데이터 부재 시 즉시 Exception 발생.
    """
    
    # [CHECK] 필수 컬럼 및 인자 검증
    if value_col not in df.columns:
        raise KeyError(f"Column '{value_col}' not found in DataFrame.")
    if "chrom" not in df.columns:
        raise KeyError("Column 'chrom' not found in DataFrame.")
    
    # [MODIFIED] 상염색체 중앙값 계산 (Baseline)
    autosome_df = df[~df["chrom"].isin(["chrX", "chrY"])]
    if autosome_df.empty:
        raise ValueError("No autosomal data found to establish global baseline.")
        
    global_median = autosome_df[value_col].median()
    log(f"Global Baseline (Autosomal Median): {global_median:.4f}")

    # [MODIFIED] 성별 판별 로직 (기준값 예외처리 포함)
    y_bins = df[df["chrom"] == "chrY"]
    if y_bins.empty:
        log("No chrY data found. Inferring as Female (or data loss).")
        y_median = -np.inf
    else:
        y_median = y_bins[value_col].median()
        
    is_male_cell = y_median > (global_median + sex_threshold)
    gender_tag = "Male" if is_male_cell else "Female"
    log(f"Inferred Sex: {gender_tag} (chrY Median: {y_median:.4f}, Threshold: {sex_threshold})")

    norm_results = []
    
    # [MODIFIED] 염색체별 루프 진행
    for chrom, group in df.groupby("chrom", sort=False):
        # Leave-one-out 방식의 참조군 설정
        ref_df = df[~df["chrom"].isin(["chrX", "chrY", chrom])]
        
        # 만약 참조할 상염색체가 부족할 경우 글로벌 중앙값 강제 사용
        current_median = ref_df[value_col].median() if not ref_df.empty else global_median
        
        # [MODIFIED] Log2 Space Target 설정 (2n=0, 1n=-1)
        target_log2 = 0.0 # Default for Autosomes (2n)
                
        if chrom == "chrX":
            # 남성 1n (-1.0), 여성 2n (0.0)
            target_log2 = -1.0 if is_male_cell else 0.0
            
        elif chrom == "chrY":
            if is_male_cell:
                # 남성 Y 1n (-1.0)
                target_log2 = -1.0
            else:
                # [MODIFIED] 여성 Y: 조작 없이 글로벌 기준 대비 노이즈 유지
                # 2n Baseline(0.0) 상태를 유지하기 위해 global_median만 감쇄
                group["log2_chrom_norm"] = group[value_col] - global_median
                norm_results.append(group)
                log(f"Processed {chrom:5}: Noise preserved relative to global.")
                continue
        
        # 정규화 계산: (관측치 - 주변 중앙값) + 생물학적 기대값(Log2)
        group["log2_chrom_norm"] = group[value_col] - current_median + target_log2
        
        log(f"Normalized {chrom:5}: Target(Log2)={target_log2:4.1f}, Ref_Median={current_median:7.4f}")
        norm_results.append(group)

    return pd.concat(norm_results).reset_index(drop=True)

