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
def normalize_all_metrics_with_sex_log2(df, depth_col="raw_total_depth", qc_frag_col="qc_pass_fragments", sex_threshold=0.0001):
    """
    [REVISED PERFECT VERSION + Zero-Drop 방어]
    1. Trans 지표는 고품질 클린 자산인 qc_frag_col (qc_pass_fragments)을 분모로 정밀 제어합니다.
    2. Hetero/Homo 사이트 지표는 게놈 와이드 물리 뎁스 및 ADO 생존율 모니터링을 위해 depth_col (raw_total_depth)을 분모로 이원화합니다.
    3. 성별 판별 기반의 Log2 Space (2n=0.0, 1n=-1.0) 공식을 적용합니다.
    4. [NEW] 성염색체(chrX, Y)는 원래 물리량이 적어 0이 자주 나오므로, 플롯 붕괴 방지를 위해 대조군 5%의 안전 쿠션(Pseudocount)을 적용합니다.
    """
    # [CHECK] 필수 컬럼 검증
    target_cols = ["total_trans_fragments", "hetero_sites_count", "homo_sites_count"]
    for col in target_cols + [depth_col, qc_frag_col, "chrom"]:
        if col not in df.columns:
            raise KeyError(f"Required column '{col}' not found in DataFrame.")

    df = df.copy()
    
    # 0. 분모 Zero Division 방지 및 이원화 밀도(Density) 생성
    safe_raw_depth = df[depth_col].replace(0, 1)
    safe_qc_frags = df[qc_frag_col].replace(0, 1)
    
    # [핵심] 지표의 생물학적 성격에 따른 분모 스위칭 알고리즘
    df["density_trans"] = df["total_trans_fragments"] / safe_qc_frags   # QC 패스 분자 풀 기준
    df["density_hetero"] = df["hetero_sites_count"] / safe_raw_depth   # 원천 뎁스(ADO 생존 모니터링) 기준
    df["density_homo"] = df["homo_sites_count"] / safe_raw_depth       # 원천 뎁스(CNV 백본 동기화) 기준

    # [CHECK] 상염색체 중앙값 계산 (Baseline 기준점 설정)
    autosome_df = df[~df["chrom"].isin(["chrX", "chrY"])]
    if autosome_df.empty:
        raise ValueError("No autosomal data found to establish global baseline.")
        
    global_med_trans = max(autosome_df["density_trans"].median(), 1e-8)
    global_med_hetero = max(autosome_df["density_hetero"].median(), 1e-8)
    global_med_homo = max(autosome_df["density_homo"].median(), 1e-8)
    
    log(f"Global Baseline - Trans Density (QC-based):  {global_med_trans:.6f}")
    log(f"Global Baseline - Hetero Density (Raw-based): {global_med_hetero:.6f}")
    log(f"Global Baseline - Homo Density (Raw-based):   {global_med_homo:.6f}")

    # [CHECK] 성별 판별 로직
    y_bins = df[df["chrom"] == "chrY"]
    if y_bins.empty:
        log("No chrY data found. Inferring as Female (or data loss).")
        y_median_trans = -np.inf
    else:
        y_median_trans = y_bins["density_trans"].median()
        
    is_male_cell = y_median_trans > (global_med_trans * sex_threshold)
    gender_tag = "Male" if is_male_cell else "Female"
    log(f"Inferred Sex: {gender_tag} (chrY Trans Median: {y_median_trans:.6f})")

    norm_results = []
    
    # 염색체별 정규화 루프 진행 (Leave-one-out 대조군 매칭)
    for chrom, group in df.groupby("chrom", sort=False):
        ref_df = df[~df["chrom"].isin(["chrX", "chrY", chrom])]
        
        ref_med_trans = max(ref_df["density_trans"].median() if not ref_df.empty else global_med_trans, 1e-8)
        ref_med_hetero = max(ref_df["density_hetero"].median() if not ref_df.empty else global_med_hetero, 1e-8)
        ref_med_homo = max(ref_df["density_homo"].median() if not ref_df.empty else global_med_homo, 1e-8)
        
        # -----------------------------------------------------------------
        # [NEW] 염색체 타입에 따른 하한선(Pseudocount) 이원화
        # -----------------------------------------------------------------
        if chrom in ["chrX", "chrY"]:
            # 성염색체: 본래 물리량이 적어 0이 자주 나오므로 대조군의 5%를 하한선으로 방어 (그래프 붕괴 방지)
            p_trans = ref_med_trans * 0.05
            p_hetero = ref_med_hetero * 0.05
            p_homo = ref_med_homo * 0.05
        else:
            # 상염색체: 0이 나오면 진짜 결실(Deletion)이므로, 순수 수학적 에러만 막는 극소값 사용
            p_trans = 1e-10
            p_hetero = 1e-10
            p_homo = 1e-10
        # -----------------------------------------------------------------
        
        target_log2 = 0.0 
                
        if chrom == "chrX":
            target_log2 = -1.0 if is_male_cell else 0.0
            
        elif chrom == "chrY":
            if is_male_cell:
                target_log2 = -1.0
            else:
                # 여성 Y 노이즈 보존 처리 (성염색체 전용 p_value 적용)
                group["trans_log2_norm"] = np.log2((group["density_trans"] + p_trans) / ref_med_trans)
                group["hetero_log2_norm"] = np.log2((group["density_hetero"] + p_hetero) / ref_med_hetero)
                group["homo_log2_norm"] = np.log2((group["density_homo"] + p_homo) / ref_med_homo)
                norm_results.append(group)
                continue
        
        # 각 지표 고유의 대조군 비율 계산 후 생물학적 축(target_log2) 동기화 (이원화 p_value 적용)
        group["trans_log2_norm"] = np.log2((group["density_trans"] + p_trans) / ref_med_trans) + target_log2
        group["hetero_log2_norm"] = np.log2((group["density_hetero"] + p_hetero) / ref_med_hetero) + target_log2
        group["homo_log2_norm"] = np.log2((group["density_homo"] + p_homo) / ref_med_homo) + target_log2
        
        norm_results.append(group)

    return pd.concat(norm_results).reset_index(drop=True)