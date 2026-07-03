import numpy as np
import pandas as pd

from utils import log
from rules import CFG
from calling_sex import determine_fetal_sex


def apply_qc_bin(df, min_depth=None, min_coverage=None):
    """
    GC/Mappability 필터(is_filtered)를 통과하고,
    Coverage 조건(Depth, Breadth)을 만족하는 정상 Bin만 추출합니다.
    단, 성염색체(chrX, chrY)는 Coverage 필터에서 예외로 두어 삭제되지 않도록 보호합니다.

    min_depth / min_coverage 를 명시적으로 넘기지 않으면 config.yaml 의
    CFG["min_depth"] / CFG["min_coverage"] 를 기본값으로 사용합니다.
    """
    if df is None or df.empty:
        return pd.DataFrame()

    if min_depth is None:
        min_depth = CFG["min_depth"]
    if min_coverage is None:
        min_coverage = CFG["min_coverage"]

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
    """상염색체 기준 LOO 정규화 및 생물학적 성별 비율 계산 (raw_count 스케일, 레거시 경로)"""
    autosome_df = df[~df["chrom"].isin(["chrX", "chrY"])]
    global_median = autosome_df[value_col].median()

    y_bins = df[df["chrom"] == "chrY"]
    y_median = y_bins[value_col].median() if not y_bins.empty else 0.0
    y_ratio = y_median / global_median if global_median > 0 else 0.0

    x_bins = df[df["chrom"] == "chrX"]
    x_median = x_bins[value_col].median() if not x_bins.empty else 0.0
    x_ratio = x_median / global_median if global_median > 0 else 0.0

    # [주의] 이 함수는 raw_count 비율(0~1대) 스케일이라 sex_calling의 CN 스케일
    # 기준(x_female_min_cn=1.5 등)과 단위가 다릅니다. 여기서는 기존 동작 보존을
    # 위해 y_presence 자체 판정을 유지합니다. (파이프라인 본 경로는
    # normalize_by_chrom_with_sex 를 사용하며 sex_calling 기준을 그대로 씁니다.)
    is_male = y_ratio > 0.15
    sex_tag = "XY" if is_male else "XX"

    norm_results = []
    for chrom, group in df.groupby("chrom", sort=False):
        ref_df = df[~df["chrom"].isin(["chrX", "chrY", chrom])]
        current_median = ref_df[value_col].median() if not ref_df.empty else global_median

        # 관측된 생물학적 Log2 비율을 그대로 저장합니다. (남성 XY의 경우 알아서 -1.0으로 계산됨)
        observed_log2 = np.log2((group[value_col] + 1e-9) / (current_median + 1e-9))
        group["log2_chrom_norm"] = observed_log2

        norm_results.append(group)

    return pd.concat(norm_results).reset_index(drop=True), sex_tag, x_ratio, y_ratio


def get_mad(data):
    """결측치를 제외한 Median Absolute Deviation 계산"""
    med = np.nanmedian(data)
    return np.nanmedian(np.abs(data - med))

def determine_fetal_sex(x_median, y_median):
    """임시 성별 판별 헬퍼 (필요시 diagnosis 모듈 로직과 동기화)"""
    if x_median < 1.5 and y_median > 0.2:
        return "XY"
    return "XX"

def normalize_by_chrom_with_sex(df, value_col="log2_corrected"):
    """
    [MODIFIED] 변경 이유: Chr19 등 GC-rich 염색체의 WGA Dropout Bias 극복.
    전체 Bin을 사용하지 않고, MAD 기반으로 통계적 유효 Bin(Passed Regions)만 
    추출하여 LOO Baseline을 계산하는 Robust 방식으로 개선함.
    """
    if df is None or df.empty:
        raise ValueError("[Error] Input DataFrame is empty before normalization.")
        
    if value_col not in df.columns:
        raise KeyError(f"[Error] Required column '{value_col}' not found in DataFrame.")

    df_norm = df.copy()

    # 1. Log2 공간 -> 선형(Linear) 카운트 공간 복원
    if "log2" in value_col:
        linear_vals = 2 ** df_norm[value_col]
    else:
        linear_vals = df_norm[value_col]

    df_norm["loo_baseline"] = np.nan
    auto_mask = ~df_norm["chrom"].isin(["chrX", "chrY"])

    # -----------------------------------------------------------------
    # [NEW] Robust Valid Bin 색출 (QC Passed Regions)
    # 극단적인 Dropout(0에 수렴)이나 증폭 에러(Spike)를 Baseline 계산에서 배제
    # -----------------------------------------------------------------
    global_auto_vals = linear_vals[auto_mask].dropna().values
    if len(global_auto_vals) == 0:
        raise ValueError("[Error] No valid autosome bins found. Cannot compute baseline.")

    global_med = np.nanmedian(global_auto_vals)
    global_mad = get_mad(global_auto_vals)

    # 중심값 기준 ± 3 MAD 이내이며, 아예 0으로 죽어버린 Bin(Dropout) 제외
    valid_mask = (
        (linear_vals > 0.1) &  # WGA Dropout 배제
        (linear_vals >= (global_med - 3 * global_mad)) & 
        (linear_vals <= (global_med + 3 * global_mad))
    )

    valid_auto_mask = auto_mask & valid_mask
    robust_global_auto_median = linear_vals[valid_auto_mask].median()

    if np.isnan(robust_global_auto_median) or robust_global_auto_median == 0:
        raise ValueError("[Error] Robust global median is invalid (NaN or 0) after filtering.")

    # 2. Leave-One-Out (LOO) Baseline 동적 할당
    for chrom in df_norm["chrom"].unique():
        chrom_mask = df_norm["chrom"] == chrom

        if chrom not in ["chrX", "chrY"]:
            # 자기 자신을 제외한 "유효한(Valid)" 상염색체만의 중앙값 추출
            loo_mask = valid_auto_mask & (~chrom_mask)
            
            # 유효 Bin이 충분한지 확인 (기본값 세팅 불가 원칙)
            if loo_mask.sum() < 100:
                raise ValueError(f"[Error] Not enough valid baseline bins (<100) to normalize {chrom}.")
                
            loo_median = linear_vals[loo_mask].median()
            df_norm.loc[chrom_mask, "loo_baseline"] = loo_median
        else:
            # 성염색체(chrX, chrY)는 "유효한" 전체 상염색체 중앙값 적용
            df_norm.loc[chrom_mask, "loo_baseline"] = robust_global_auto_median

    # 3. Y축 기준점을 2.0으로 하는 'Copy Number Signal' 생성
    df_norm["copy_number_signal"] = 2.0 * (linear_vals / df_norm["loo_baseline"])

    # 4. 하위 모듈 호환성을 위한 log2_chrom_norm 생성
    with np.errstate(divide='ignore', invalid='ignore'):
        df_norm["log2_chrom_norm"] = np.log2(df_norm["copy_number_signal"] / 2.0)
    
    # Inf 및 NaN 값 예외처리 (연산 에러 방지)
    df_norm["log2_chrom_norm"] = df_norm["log2_chrom_norm"].replace([np.inf, -np.inf], np.nan).fillna(-3.0)

    df_norm = df_norm.drop(columns=["loo_baseline"])

    # 5. 성별 판별
    x_mask = df_norm["chrom"] == "chrX"
    y_mask = df_norm["chrom"] == "chrY"

    x_median = df_norm.loc[x_mask, "copy_number_signal"].median() if x_mask.sum() > 0 else 0
    y_median = df_norm.loc[y_mask, "copy_number_signal"].median() if y_mask.sum() > 0 else 0

    sex_tag = determine_fetal_sex(x_median, y_median)

    log(f"Normalization Mode: Robust LOO | Valid Bins Used: {valid_auto_mask.sum()} | Sex: {sex_tag} (X: {x_median:.2f}, Y: {y_median:.2f})")

    return df_norm, sex_tag, float(x_median / 2.0), float(y_median / 2.0)