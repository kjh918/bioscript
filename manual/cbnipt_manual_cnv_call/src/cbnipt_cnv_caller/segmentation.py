import numpy as np
import pandas as pd
import ruptures as rpt
import os, sys
 
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

from utils import log
def segment_one_cell(meta, df_signals, signal_col="log2_chrom_norm", penalty=10.0):
    """
    [SIMPLIFIED & MEDIAN FIXED] 순수 Copy Number (Log2FC) 신호 기반 Change-point 탐지.
    단일 세포 WGA의 극단적 노이즈를 방어하기 위해 평균(Mean) 대신 중앙값(Median)을 사용합니다.
    """
    if penalty is None:
        raise ValueError("Penalty value must be provided for PELT algorithm.")

    segments = []
    chrom_list = meta["chrom"].unique()
    
    for chrom in chrom_list:
        idx = meta.index[meta["chrom"] == chrom]
        y = df_signals.loc[idx, signal_col].values.astype(float).flatten()
        
        valid = np.isfinite(y)
        if valid.sum() < 10: 
            continue
            
        y_v = y[valid].reshape(-1, 1)
        v_idx = np.where(valid)[0]
        
        # PELT 알고리즘 실행
        algo = rpt.Pelt(model="l2").fit(y_v)
        bkps = algo.predict(pen=penalty)
        
        start_idx = 0
        for bp in bkps:
            # 다음 세그먼트와 겹치거나 빈틈이 없도록 타일링
            if bp < len(v_idx):
                end_idx = v_idx[bp - 1]
            else:
                end_idx = len(idx) - 1
                
            full_seg_idx = idx[start_idx : end_idx + 1]
            
            if len(full_seg_idx) > 0:
                seg_vals = y[start_idx : end_idx + 1]
                
                # [CRITICAL FIX] np.nanmean -> np.nanmedian 으로 변경!
                # 노이즈 구덩이(-3, -4)에 휘둘리지 않고 진짜 베이스라인을 찾습니다.
                seg_median = 0.0 if np.isnan(seg_vals).all() else float(np.nanmedian(seg_vals))
                    
                segments.append({
                    "chrom": chrom,
                    "start": int(meta.loc[full_seg_idx[0], "start"]),
                    "end": int(meta.loc[full_seg_idx[-1], "end"]),
                    "n_bins": len(full_seg_idx),
                    "seg_median": seg_median  # 이름도 median으로 변경
                })
            
            start_idx = end_idx + 1
            
    if not segments:
        log("Warning: No segments were generated. Check input signal or penalty.")
        return pd.DataFrame()
        
    return pd.DataFrame(segments)


def assign_cn_state(seg_df, baseline_ploidy=2, sex_tag="UNKNOWN"):
    """
    [SIMPLIFIED & FIXED] 1D seg_mean 기반 상태 할당.
    성별(Sex)에 따라 성염색체(chrX, chrY)의 기본 Ploidy를 동적으로 조절합니다.
    """
    if seg_df is None or seg_df.empty:
        return seg_df
        
    df = seg_df.copy()
    
    if "seg_mean" not in df.columns:
        raise KeyError("Cannot find 'seg_mean' column in segments dataframe.")

    # 1. Log2FC -> 실제 비율 변환
    ratios = np.exp2(df["seg_mean"].fillna(0.0))
    
    # 2. 성별에 따른 염색체별 예상 Copy Number (Ploidy) 설정
    expected_ploidy = np.full(len(df), baseline_ploidy, dtype=int)
    chrom_array = df["chrom"].values
    
    if sex_tag == "XY":
        expected_ploidy[chrom_array == "chrX"] = 1
        expected_ploidy[chrom_array == "chrY"] = 1
    elif sex_tag == "XX":
        expected_ploidy[chrom_array == "chrY"] = 0
        
    # 3. 정수형 Copy Number 추정
    df["copy_number"] = np.round(ratios * expected_ploidy).astype(int)
    df["copy_number"] = df["copy_number"].clip(lower=0) # 음수 방지
    
    # 4. 직관적인 CNV Call
    conditions = [
        df["copy_number"] > expected_ploidy,
        df["copy_number"] < expected_ploidy
    ]
    choices = ["AMP", "DEL"]
    df["cnv_call"] = np.select(conditions, choices, default="NEUT")
    
    return df