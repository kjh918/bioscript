import numpy as np
import pandas as pd

from utils import log
from rules import CFG
from sex_calling import determine_fetal_sex


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


def normalize_by_chrom_with_sex(df, value_col="log2_corrected"):
    """
    Leave-One-Out (LOO) 정규화 적용.
    각 상염색체를 정규화할 때, 자기 자신을 제외한 나머지 상염색체들의 중앙값을 기준으로 삼아
    거대 염색체의 Aneuploidy가 전체 기저선(Baseline)을 왜곡하는 것을 방지합니다.
    """
    df_norm = df.copy()

    # 1. Log2 공간에 있는 보정값을 선형(Linear) 카운트 공간으로 복원
    if "log2" in value_col:
        linear_vals = 2 ** df_norm[value_col]
    else:
        linear_vals = df_norm[value_col]

    df_norm["loo_baseline"] = 1e-9  # 임시 기저선 컬럼
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

    # 3-1. 염색체별 체계적 편향(bias) 보정 (예: chr19 GC-rich 잔여 편향)
    # CFG["CHROM_BIAS_CORRECTION"] 에 등록된 염색체만 보정하며, 값은 log2 스케일
    # 오프셋입니다 (정상 코호트에서 관측된 log2_chrom_norm 중앙값을 빼주는 것과 동일).
    bias_cfg = CFG.get("CHROM_BIAS_CORRECTION", {})
    for chrom, log2_offset in bias_cfg.items():
        if log2_offset == 0.0:
            continue
        chrom_mask = df_norm["chrom"] == chrom
        if chrom_mask.sum() == 0:
            continue
        # copy_number_signal = 2^(log2 + offset) 형태로 보정 (offset 만큼 위로 밀어줌)
        df_norm.loc[chrom_mask, "copy_number_signal"] = (
            df_norm.loc[chrom_mask, "copy_number_signal"] * (2.0 ** log2_offset)
        )

    # 4. 하위 모듈(Segmentation) 호환성을 위해 log2_chrom_norm 도 생성
    with np.errstate(divide='ignore'):
        df_norm["log2_chrom_norm"] = np.log2(df_norm["copy_number_signal"] / 2.0)
    df_norm["log2_chrom_norm"] = df_norm["log2_chrom_norm"].replace([-np.inf], -3.0)

    # 임시 컬럼 삭제
    df_norm = df_norm.drop(columns=["loo_baseline"])

    # 5. 성별 판별 (sex_calling 단일 로직 사용 — 기존 하드코딩된 1.5/0.2 대체)
    x_mask = df_norm["chrom"] == "chrX"
    y_mask = df_norm["chrom"] == "chrY"

    x_median = df_norm.loc[x_mask, "copy_number_signal"].median() if x_mask.sum() > 0 else 0
    y_median = df_norm.loc[y_mask, "copy_number_signal"].median() if y_mask.sum() > 0 else 0

    sex_tag = determine_fetal_sex(x_median, y_median)
    if sex_tag == "Unknown":
        sex_tag = "XX" if x_median >= y_median else "XY"

    log(f"Normalization Mode: Leave-One-Out (LOO) | Sex: {sex_tag} (X Copy: {x_median:.2f}, Y Copy: {y_median:.2f})")

    return df_norm, sex_tag, float(x_median / 2.0), float(y_median / 2.0)