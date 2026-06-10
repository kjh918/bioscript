import numpy as np
import pandas as pd
import ruptures as rpt
from utils import log
#from scipy.signal import find_peaks, gaussian_filter1d

def segment_one_cell(meta, signal, penalty):
    """
    [MODIFIED] 단일 세포 WGA 샘플의 Log2 신호를 기반으로 Change-point detection을 수행합니다.
    Reference: Killick, R. et al. (2012) Optimal Detection of Changepoints. JASA. 
    (https://doi.org/10.1080/01621459.2012.737745)
    """
    # [MODIFIED] Default 값 제거 및 필수 인자 검증
    if penalty is None:
        raise ValueError("Penalty value must be provided for PELT algorithm.")
    if meta is None or signal is None:
        raise ValueError("Input meta and signal dataframes cannot be None.")
    if len(meta) != len(signal):
        raise ValueError("Metadata and Signal length mismatch.")

    segments = []
    # [MODIFIED] 데이터 일관성을 위해 index 정렬 보장
    chrom_list = meta["chrom"].unique()
    
    for chrom in chrom_list:
        idx = meta.index[meta["chrom"] == chrom]
        # [MODIFIED] signal이 Series인 경우와 DataFrame인 경우 모두 대응
        y = signal.loc[idx].values.astype(float).flatten()
        
        valid = np.isfinite(y)
        if valid.sum() < 10: 
            log(f"Skipping {chrom}: insufficient valid bins ({valid.sum()})")
            continue
            
        y_v = y[valid].reshape(-1, 1)
        v_idx = np.where(valid)[0]
        
        # PELT 알고리즘: L2(Gaussian) 모델 사용
        algo = rpt.Pelt(model="l2").fit(y_v)
        bkps = algo.predict(pen=penalty)
        
        start = 0
        for bp in bkps:
            seg_v = v_idx[start:bp]
            if len(seg_v) == 0: continue
            
            g_idx = idx[seg_v]
            
            # [MODIFIED] 정수형 좌표 및 통계치 계산 최적화
            segments.append({
                "chrom": chrom,
                "start": int(meta.loc[g_idx[0], "start"]),
                "end": int(meta.loc[g_idx[-1], "end"]),
                "n_bins": len(g_idx),
                # [MODIFIED] 기존 signal.loc[g_idx]의 평균을 직접 계산
                "seg_mean": float(np.nanmean(y[seg_v]))
            })
            start = bp
            
    if not segments:
        raise RuntimeError("No segments were generated. Check input signal or penalty.")
        
    return pd.DataFrame(segments)

def estimate_ploidy(seg_df, max_ploidy=4):
    """
    Segmentation 결과의 seg_mean 분포를 기반으로 
    가장 가능성 높은 Ploidy(보통 2 또는 3)를 추정합니다.
    """
    # 1. seg_mean 분포 추출 (가중치: 각 세그먼트의 bin 개수)
    values = seg_df["seg_mean"].values
    weights = seg_df["n_bins"].values
    
    # 2. 히스토그램 생성 (Log2 Space에서 2n은 0.0 부근에 피크)
    hist, bin_edges = np.histogram(values, bins=50, weights=weights, density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    # 3. 노이즈 완화를 위한 가우시안 평탄화
    smoothed = gaussian_filter1d(hist, sigma=1)
    peaks, _ = find_peaks(smoothed, height=np.max(smoothed)*0.1)
    
    if len(peaks) == 0:
        return 2 # 기본값
    
    # 4. 가장 높은 피크의 위치가 0(Ratio 1.0)에 가장 가까운 Ploidy를 찾음
    best_peak_val = bin_centers[peaks[np.argmax(smoothed[peaks])]]
    
    # 2^best_peak_val 이 Ratio인데, 2n이 1.0(Log2 0)이라고 가정할 때
    # Ploidy = 2 * 2^(best_peak_val) 에 가장 가까운 정수를 찾음
    estimated = np.round(2 * (2 ** best_peak_val))
    return int(np.clip(estimated, 1, max_ploidy))

def assign_cn_state(seg_df, baseline_ploidy):
    """
    [MODIFIED] Log2 Space 기준을 반영하여 정수형 Copy Number를 할당합니다.
    Target: 2n = 0 (Log2 ratio 1.0)
    Formula: CN = baseline_ploidy * 2^(seg_mean)
    """
    # [MODIFIED] Default 값 제거
    if baseline_ploidy is None:
        raise ValueError("baseline_ploidy must be explicitly defined (e.g., 2).")
        
    if seg_df.empty:
        return seg_df

    # [MODIFIED] 사용자 요청 Log2 Space 반영: 2n=0 -> ratio=1
    # 2n * 2^(0) = 2, 2n * 2^(-1) = 1, 2n * 2^(0.58) = 3
    ratios = 2 ** seg_df["seg_mean"]
    seg_df["copy_number"] = np.clip(
        np.round(baseline_ploidy * ratios).astype(int), 
        0, 
        None
    )
    
    return seg_df