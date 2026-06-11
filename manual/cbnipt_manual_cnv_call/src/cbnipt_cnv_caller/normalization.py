import numpy as np
import pandas as pd
from .utils import log
from .rules import CFG


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

def normalize_by_chrom_with_sex(df, value_col="raw_count", sex_ratio_thresh=0.15):
    """
    원시 TSV 데이터를 읽어 성별 추정 및 Leave-One-Out(LOO) Log2 정규화를 수행합니다.
    """
    if value_col not in df.columns:
        raise KeyError(f"Column '{value_col}' not found in DataFrame.")
    if "chrom" not in df.columns:
        raise KeyError("Column 'chrom' not found in DataFrame.")
    
    # 1. 상염색체 중앙값 계산 (Global Baseline)
    autosome_df = df[~df["chrom"].isin(["chrX", "chrY"])]
    if autosome_df.empty:
        raise ValueError("No autosomal data found to establish global baseline.")
        
    global_median = autosome_df[value_col].median()
    log(f"Global Baseline (Autosomal Median): {global_median:.4f}")

    # 2. 성별 판별 로직 (비율 기반)
    y_bins = df[df["chrom"] == "chrY"]
    if y_bins.empty or global_median == 0:
        log("No chrY data or zero baseline found. Inferring as Female.")
        y_median = 0.0
        y_ratio = 0.0
    else:
        y_median = y_bins[value_col].median()
        y_ratio = y_median / global_median
        
    is_male_cell = y_ratio > sex_ratio_thresh
    gender_tag = "Male (XY)" if is_male_cell else "Female (XX)"
    log(f"Inferred Sex: {gender_tag} (chrY/Auto Ratio: {y_ratio:.3f}, Thresh: {sex_ratio_thresh})")

    norm_results = []
    
    # 3. 염색체별 LOO 정규화 루프
    for chrom, group in df.groupby("chrom", sort=False):
        ref_df = df[~df["chrom"].isin(["chrX", "chrY", chrom])]
        current_median = ref_df[value_col].median() if not ref_df.empty else global_median
        
        # Log2 Space Target (Expected Log2FC)
        expected_log2 = 0.0 # Default for Autosomes (2n)
                
        if chrom == "chrX":
            expected_log2 = -1.0 if is_male_cell else 0.0 # XY: 1n, XX: 2n
            
        elif chrom == "chrY":
            if is_male_cell:
                expected_log2 = -1.0 # XY: 1n
            else:
                # 여성 Y: 데이터 소실/노이즈 보존 (관측치 Log2FC 그대로 투영)
                group["log2_chrom_norm"] = np.log2((group[value_col] + 1e-9) / (global_median + 1e-9))
                norm_results.append(group)
                log(f"Processed {chrom:5}: Noise preserved relative to global.")
                continue
        
        # [수학적 교정] (관측치 / 중앙값)의 Log2 - (생물학적 기대 Log2)
        # 이렇게 하면, 남성 chrX라도 정상(1-copy)이면 log2_chrom_norm이 0.0으로 맞춰집니다.
        observed_log2 = np.log2((group[value_col] + 1e-9) / (current_median + 1e-9))
        group["log2_chrom_norm"] = observed_log2 - expected_log2
        
        log(f"Normalized {chrom:5}: Expected(Log2)={expected_log2:4.1f}, Ref_Median={current_median:7.1f}")
        norm_results.append(group)

    return pd.concat(norm_results).reset_index(drop=True), gender_tag, y_ratio



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