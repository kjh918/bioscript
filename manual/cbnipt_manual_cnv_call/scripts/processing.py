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
    [REVISED] 
    상염색체(Autosomes)는 Breadth of Coverage와 Depth를 기준으로 불량 빈을 엄격히 제거하지만,
    성염색체(chrX, chrY)는 생물학적 성별/Plody 판정을 위해 필터링 대상에서 완전히 예외(보호) 처리합니다.
    """
    if df is None or df.empty:
        return pd.DataFrame(), pd.DataFrame([filter_stats])
        
    initial_count = len(df)
    
    # [STEP 1] 성염색체와 상염색체 영역 분리
    sex_chrom_mask = df["chrom"].isin(["chrX", "chrY"])
    df_autosomes = df[~sex_chrom_mask].copy()
    df_sex_chroms = df[sex_chrom_mask].copy()
    
    # -----------------------------------------------------------------
    # [STEP 2] 상염색체(Autosomes)에만 엄격한 퀄리티 필터 적용
    # -----------------------------------------------------------------
    autosome_initial_count = len(df_autosomes)
    
    # 1. Depth 필터링
    df_auto_after_depth = df_autosomes[df_autosomes["raw_count"] >= min_depth].copy()
    filter_stats["low_depth_count"] = autosome_initial_count - len(df_auto_after_depth)
    
    # 2. Breadth 필터링
    df_auto_filtered = df_auto_after_depth[df_auto_after_depth["breadth_ratio"] >= min_coverage].copy()
    filter_stats["low_breadth_bins_count"] = len(df_auto_after_depth) - len(df_auto_filtered)
    
    # -----------------------------------------------------------------
    # [STEP 3] 필터링된 상염색체 데이터와 '온전한 성염색체 데이터' 결합
    # -----------------------------------------------------------------
    df_filtered = pd.concat([df_auto_filtered, df_sex_chroms], ignore_index=True)
    
    # 원래 분석 파이프라인의 순서 정렬 상태(염색체 및 시작 좌표 순) 복원
    def chrom_key(c):
        c_str = str(c).replace('chr', '')
        if c_str == 'X': return 23
        if c_str == 'Y': return 24
        try: return int(c_str)
        except: return 99

    df_filtered['chrom_sort_key'] = df_filtered['chrom'].apply(chrom_key)
    df_filtered = df_filtered.sort_values(by=['chrom_sort_key', 'start']).drop(columns=['chrom_sort_key']).reset_index(drop=True)
    
    total_removed = initial_count - len(df_filtered)
    log(f"Applied localized quality filter (Sex Chromosomes Protected): {total_removed} bins removed "
        f"(Autosome Depth < {min_depth}: {filter_stats['low_depth_count']}, "
        f"Autosome Breadth < {min_coverage*100}%: {filter_stats['low_breadth_bins_count']})")
    
    return df_filtered, pd.DataFrame([filter_stats])

def gc_correct_lowess(df, frac=0.2):
    x, y = df["gc"].values, np.log2(df["raw_count"].values + 1)
    valid = np.isfinite(x) & (y > 0)
    fit = lowess(y[valid], x[valid], frac=frac, return_sorted=False)
    df["log2_corrected"] = np.nan
    df.loc[valid, "log2_corrected"] = y[valid] - fit + np.nanmedian(fit)
    return df, (x[valid], y[valid], fit)