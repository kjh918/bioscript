import numpy as np
import pandas as pd
import pysam
import os 
from utils import log

from statsmodels.nonparametric.smoothers_lowess import lowess

def process_bam_to_coverage(bam_path, bins_df, min_mapq=20):
    """
    한 번의 BAM 탐색으로 리드 카운트와 실제 점유 면적(Breadth)을 동시에 계산합니다.
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    
    # 필터링 통계 초기화
    filter_stats = {
        "total_fetched_reads": 0,
        "filtered_unmapped": 0,
        "filtered_duplicate": 0,
        "filtered_low_mapq": 0,
        "filtered_secondary": 0,
        "filtered_by_midpoint": 0,
        "passed_reads": 0,
        "low_breadth_bins_count": 0  # 50% 미만 점유 빈 카운트
    }
    
    counts = []
    breadth_ratios = []

    for r in bins_df.itertuples():
        bin_len = r.end - r.start
        # 이번 빈의 영역을 나타내는 boolean 배열 (False로 초기화)
        coverage_mask = np.zeros(bin_len, dtype=bool)
        bin_count = 0
        
        # BAM fetch: 한 번의 루프에서 모든 것을 처리
        for read in bam.fetch(r.chrom, r.start, r.end):
            filter_stats["total_fetched_reads"] += 1
            
            # 1. 기본 필터링
            if read.is_unmapped:
                filter_stats["filtered_unmapped"] += 1
                continue
            if read.is_duplicate:
                filter_stats["filtered_duplicate"] += 1
                continue
            if read.mapping_quality < min_mapq:
                filter_stats["filtered_low_mapq"] += 1
                continue
            if read.is_secondary or read.is_supplementary:
                filter_stats["filtered_secondary"] += 1
                continue
            
            # 2. Breadth of Coverage 계산 (실제 리드가 덮은 염기 마킹)
            # 리드의 시작/끝이 빈의 경계를 벗어날 수 있으므로 자르기(clipping)
            ref_start = max(read.reference_start, r.start)
            ref_end = min(read.reference_end, r.end)
            
            if ref_start < ref_end:
                # 빈 내부의 상대 좌표로 변환하여 마킹
                local_start = ref_start - r.start
                local_end = ref_end - r.start
                coverage_mask[local_start:local_end] = True
            
            # 3. Midpoint Counting (CNV 분석용 리드 수 산출)
            mid = (read.reference_start + read.reference_end) // 2
            if r.start <= mid < r.end:
                bin_count += 1
                filter_stats["passed_reads"] += 1
            else:
                filter_stats["filtered_by_midpoint"] += 1
                
        # 이번 빈의 최종 Breadth 계산
        actual_covered_bases = np.sum(coverage_mask)
        breadth_ratio = actual_covered_bases / bin_len
        
        counts.append(bin_count)
        breadth_ratios.append(breadth_ratio)
        
            
    bam.close()
    
    # 데이터프레임에 결과 저장
    bins_df["raw_count"] = counts
    bins_df["breadth_ratio"] = breadth_ratios
    
    # [참고] 나중에 apply_final_filters에서 breadth_ratio < 0.5인 행을 drop하게 됩니다.
    return bins_df, filter_stats

def apply_low_quality_filter(df, filter_stats, min_depth=1, min_coverage=0.5):
    """
    Breadth of Coverage와 Depth를 기준으로 불량 빈을 제거합니다.
    순차적으로 필터링하여 각 단계에서 '추가로' 제거된 개수를 기록합니다.
    """
    initial_count = len(df)
    
    # 1. Depth 필터링 (데이터가 아예 없거나 너무 적은 경우)
    df_after_depth = df[df["raw_count"] >= min_depth].copy()
    filter_stats["low_depth_count"] = initial_count - len(df_after_depth)
    
    # 2. Breadth 필터링 (데이터는 있으나 특정 지점에 몰려 있는 경우)
    # [수정] 이전 단계에서 남은 df_after_depth를 기준으로 계산해야 합니다.
    df_filtered = df_after_depth[df_after_depth["breadth_ratio"] >= min_coverage].copy()
    filter_stats["low_breadth_bins_count"] = len(df_after_depth) - len(df_filtered)
    
    total_removed = initial_count - len(df_filtered)
    log(f"Applied low quality filter: {total_removed} bins removed "
        f"(Depth < {min_depth}: {filter_stats['low_depth_count']}, "
        f"Breadth < {min_coverage*100}%: {filter_stats['low_breadth_bins_count']})")
    
    return df_filtered, pd.DataFrame([filter_stats])

def gc_correct_lowess(df, frac=0.2):
    x, y = df["gc"].values, np.log2(df["raw_count"].values + 1)
    valid = np.isfinite(x) & (y > 0)
    fit = lowess(y[valid], x[valid], frac=frac, return_sorted=False)
    df["log2_corrected"] = np.nan
    df.loc[valid, "log2_corrected"] = y[valid] - fit + np.nanmedian(fit)
    return df, (x[valid], y[valid], fit)