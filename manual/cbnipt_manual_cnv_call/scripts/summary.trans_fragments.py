import pandas as pd 
from pathlib import Path
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


from normalization import normalize_all_metrics_with_sex_log2


def plot_normalized_fragments(summary_df, output_path="normalized_fragments_profile.png"):
    """
    BAF를 제외하고, Trans/Cis 프래그먼트 수치를 상염색체 중간값으로 
    정규화(Normalization)하여 세로 2단 산점도로 시각화합니다.
    """
    if summary_df.empty:
        print("Warning: Input DataFrame is empty. Skipping plot.")
        return

    df = summary_df.copy()
    
    # 1. 자연스러운 염색체 순서 정렬 (chr1~22, X, Y)
    def chrom_key(c):
        c_str = str(c).replace('chr', '')
        if c_str == 'X': return 23
        if c_str == 'Y': return 24
        try: return int(c_str)
        except: return 99
    
    df['chrom_sort'] = df['chrom'].apply(chrom_key)
    df = df.sort_values(['chrom_sort', 'start']).reset_index(drop=True)

    # 2. [CRITICAL] 상염색체(Autosome) 기준 중간값 계산 (성염색체 제외)
    autosomes = df[~df['chrom'].isin(['chrX', 'chrY'])]
    
    median_trans = autosomes['total_trans_fragments'].median()
    median_cis = autosomes['total_cis_fragments'].median()
    
    if median_trans == 0: median_trans = 1.0
    if median_cis == 0: median_cis = 1.0

    # 3. 프래그먼트 지표 정규화 실행 (Autosome Median Scaling)
    df['trans_norm'] = df['total_trans_fragments'] / median_trans
    df['cis_norm'] = df['total_cis_fragments'] / median_cis

    # 4. 연속적인 X축 게놈 좌표(Cumulative Position) 계산
    cumulative_pos = []
    current_offset = 0
    last_chrom = None
    chrom_boundaries = []
    chrom_ticks = []
    chrom_labels = []

    for row in df.itertuples():
        if row.chrom != last_chrom:
            if last_chrom is not None:
                chrom_boundaries.append(current_offset)
            chrom_labels.append(str(row.chrom).replace('chr', ''))
            last_chrom = row.chrom
            
        mid_point = (row.start + row.end) / 2
        cumulative_pos.append(current_offset + mid_point)
        
        if row.Index < len(df) - 1:
            next_row = df.iloc[row.Index + 1]
            if next_row['chrom'] != row.chrom:
                current_offset += row.end

    df['genomic_pos'] = cumulative_pos
    chrom_boundaries.append(current_offset)

    for clabel in df['chrom'].unique():
        sub = df[df['chrom'] == clabel]
        chrom_ticks.append(sub['genomic_pos'].mean())

    # 5. 2단 서브플롯 구성 (plt.figure() 우회 제약 준수)
    fig, axes = plt.subplots(2, 1, figsize=(16, 10), sharex=True)

    # Panel 1: Normalized Trans Fragments
    axes[0].scatter(df['genomic_pos'], df['trans_norm'], c='#ff7f0e', alpha=0.6, s=15, label='Normalized Trans')
    # 기준선 (1.0) 및 Trisomy 예상선 (1.15) 추가
    axes[0].axhline(1.0, color='black', linestyle='-', alpha=0.5, label='Diploid Baseline (1.0)')
    axes[0].axhline(1.15, color='red', linestyle='--', alpha=0.4, label='Potential Trisomy (1.15)')
    
    axes[0].set_ylabel('Normalized Trans Factor', fontsize=12)
    axes[0].set_title(f"Genome-wide Normalized Fragment Profiles (BAF Excluded)", fontsize=16, fontweight='bold')
    axes[0].grid(True, linestyle='--', alpha=0.3)
    axes[0].legend(loc='upper right')
    axes[0].set_ylim(0, 2.0)  # 가독성을 위한 Y축 가이드 한계 설정

    # Panel 2: Normalized Cis Fragments
    axes[1].scatter(df['genomic_pos'], df['cis_norm'], c='#1f77b4', alpha=0.6, s=15, label='Normalized Cis')
    axes[1].axhline(1.0, color='black', linestyle='-', alpha=0.5)
    axes[1].set_ylabel('Normalized Cis Factor', fontsize=12)
    axes[1].set_xlabel('Chromosomes', fontsize=12)
    axes[1].grid(True, linestyle='--', alpha=0.3)
    axes[1].legend(loc='upper right')
    axes[1].set_ylim(0, 2.5)

    # 염색체 경계선 및 라벨 세팅
    for boundary in chrom_boundaries[:-1]:
        for ax in axes:
            ax.axvline(boundary, color='grey', linestyle=':', alpha=0.5)

    axes[1].set_xticks(chrom_ticks)
    axes[1].set_xticklabels(chrom_labels, rotation=0, fontsize=10)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f"Success: Normalized fragment plot saved to -> {output_path}")
    
    # 정규화된 데이터프레임 반환 (파이프라인 다음 단계 연동용)
    return df

def plot_baf_distribution_profile(summary_df, output_path="genome_wide_baf_profile.png"):
    """
    Hetero BAF(Mean/Median)와 Homo BAF의 유전체 전역 Bin별 변화 분포를
    상하 2단 패널 산점도로 시각화합니다.
    """
    if summary_df.empty:
        print("Warning: Input DataFrame is empty. Skipping plot.")
        return

    df = summary_df.copy()
    
    # 1. 자연스러운 염색체 순서 정렬 (chr1~22, X, Y)
    def chrom_key(c):
        c_str = str(c).replace('chr', '')
        if c_str == 'X': return 23
        if c_str == 'Y': return 24
        try: return int(c_str)
        except: return 99
    
    df['chrom_sort'] = df['chrom'].apply(chrom_key)
    df = df.sort_values(['chrom_sort', 'start']).reset_index(drop=True)

    # 2. 연속적인 X축 게놈 좌표(Cumulative Position) 계산
    cumulative_pos = []
    current_offset = 0
    last_chrom = None
    chrom_boundaries = []
    chrom_ticks = []
    chrom_labels = []

    for row in df.itertuples():
        if row.chrom != last_chrom:
            if last_chrom is not None:
                chrom_boundaries.append(current_offset)
            chrom_labels.append(str(row.chrom).replace('chr', ''))
            last_chrom = row.chrom
            
        mid_point = (row.start + row.end) / 2
        cumulative_pos.append(current_offset + mid_point)
        
        if row.Index < len(df) - 1:
            next_row = df.iloc[row.Index + 1]
            if next_row['chrom'] != row.chrom:
                current_offset += row.end

    df['genomic_pos'] = cumulative_pos
    chrom_boundaries.append(current_offset)

    for clabel in df['chrom'].unique():
        sub = df[df['chrom'] == clabel]
        chrom_ticks.append(sub['genomic_pos'].mean())

    # 3. 2단 서브플롯 구성 (plt.figure() 우회 제약 준수)
    fig, axes = plt.subplots(2, 1, figsize=(16, 10), sharex=True)

    # Panel 1: Heterozygous BAF 분포 (Mean vs Median 대조)
    # 2x 저심도 특성상 데이터가 유효한 빈(sites_count > 0)만 플로팅
    valid_hetero = df[df['hetero_sites_count'] > 0]
    axes[0].scatter(valid_hetero['genomic_pos'], valid_hetero['hetero_baf_mean'], 
                    c='#2ca02c', alpha=0.5, s=15, label='Hetero BAF Mean (Trend)')
    axes[0].scatter(valid_hetero['genomic_pos'], valid_hetero['hetero_baf_median'], 
                    c='#d62728', alpha=0.5, s=15, marker='x', label='Hetero BAF Median (Robust)')
    
    axes[0].set_ylabel('Observed BAF (Hetero Sites)', fontsize=12)
    axes[0].set_title("Genome-wide Allelic Frequency (BAF) Shift Profile", fontsize=16, fontweight='bold')
    axes[0].grid(True, linestyle='--', alpha=0.3)
    axes[0].legend(loc='upper right')
    axes[0].set_ylim(-0.02, 1.02)

    # Panel 2: Homozygous BAF 분포 (Background Noise 대조군)
    valid_homo = df[df['homo_sites_count'] > 0]
    axes[1].scatter(valid_homo['genomic_pos'], valid_homo['homo_baf_mean'], 
                    c='#7f7f7f', alpha=0.4, s=12, label='Homo BAF Mean (Noise Indicator)')
    
    axes[1].set_ylabel('Observed BAF (Homo Sites)', fontsize=12)
    axes[1].set_xlabel('Chromosomes', fontsize=12)
    axes[1].grid(True, linestyle='--', alpha=0.3)
    axes[1].legend(loc='upper right')
    axes[1].set_ylim(-0.02, 1.02)

    # 염색체 경계 세로선 주입
    for boundary in chrom_boundaries[:-1]:
        for ax in axes:
            ax.axvline(boundary, color='grey', linestyle=':', alpha=0.5)

    axes[1].set_xticks(chrom_ticks)
    axes[1].set_xticklabels(chrom_labels, rotation=0, fontsize=10)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f"Success: BAF distribution profile saved to -> {output_path}")

def recompute_bin_summary(position_df, hetero_range=(0.4, 0.6)):
    """
    포지션 레벨의 상세 데이터프레임(또는 저장된 파일 로드 데이터)을 입력받아
    후속 정규화 및 Ploidy 판정에 필요한 Bin 단위 요약 통계(summary_df)를 재연산합니다.
    """
    if position_df is None or position_df.empty:
        return pd.DataFrame()
    
    # 1. 안전한 데이터 처리를 위해 카피 및 컬럼 정형화
    df = position_df.copy()
    
    # 만약 기존에 group 컬럼이 없거나 범위를 수정하고 싶다면 재할당
    h_min, h_max = hetero_range
    df['group'] = 'Other'
    df.loc[(df['pop_af'] >= h_min) & (df['pop_af'] <= h_max), 'group'] = 'Hetero'
    df.loc[(df['pop_af'] <= 0.1) | (df['pop_af'] >= 0.9), 'group'] = 'Homo'
    
    summary_list = []
    
    # 2. 각 Bin ID별로 그룹화하여 통계량 재집계
    for bin_id, group in df.groupby("bin_id", sort=False):
        # 1) [CRITICAL FOR NORMALIZATION] 필터링 전 영역 내 원천 데이터 합계 추출
        raw_ref_sum = int(group["ref_depth"].sum())
        raw_alt_sum = int(group["alt_depth"].sum())
        raw_total_depth = raw_ref_sum + raw_alt_sum
        
        # 2) 유전체 매핑 위치 복원 (bin_id 파싱: "chr13:28470000-28570000")
        chrom = group["chrom"].iloc[0]
        try:
            # bin_id 형태가 깨져있을 경우를 대비해 group의 실제 값에서 백업
            start = group["pos"].min() 
            end = group["pos"].max()
        except:
            start, end = 0, 0
            
        # 3) 퀄리티 필터링 적용 (포지션당 읽힌 파편 수가 2개 이상인 곳만 BAF 통계에 반영)
        # 2x 저심도 데이터의 통계적 왜곡을 방어하기 위한 기준선 유지
        qc_group = group[group['total_fragments'] >= 2]
        
        # Groupby 연산을 이용해 Hetero/Homo 그룹별 mean, median, count 신속 추출
        grouped_stats = qc_group.groupby('group')['observed_baf'].agg(['mean', 'median', 'count'])
        
        res = {
            "bin_id": bin_id,
            "chrom": chrom,
            "start": start,
            "end": end,
            # 원천 depth 정보 바인딩 (후속 스케일링/정규화 전용 변수)
            "raw_ref_depth_sum": raw_ref_sum,
            "raw_alt_depth_sum": raw_alt_sum,
            "raw_total_depth": raw_total_depth
        }
        
        # Hetero, Homo 각각의 BAF 분포 및 마커 개수 매핑
        for g in ['Hetero', 'Homo']:
            if g in grouped_stats.index:
                res[f"{g.lower()}_baf_mean"] = grouped_stats.loc[g, 'mean']
                res[f"{g.lower()}_baf_median"] = grouped_stats.loc[g, 'median']
                res[f"{g.lower()}_sites_count"] = int(grouped_stats.loc[g, 'count'])
            else:
                res[f"{g.lower()}_baf_mean"] = res[f"{g.lower()}_baf_median"] = res[f"{g.lower()}_sites_count"] = 0
                
        # Haplotype 검증용 프래그먼트 총합 복원 (포지션별 support 중 최대값 혹은 고유 파레트의 합산 개념)
        # 여기서는 포지션 테이블에 이미 누적된 support 기반 정보를 집계합니다.
        res["total_trans_fragments"] = int(group["trans_support"].sum())
        res["total_cis_fragments"] = int(group["cis_support"].sum())
        
        # 전체 영역 뎁스 대비 고유 위상 비 (Density normalization 뼈대 변수)
        res["raw_bin_TER"] = res["total_trans_fragments"] / (raw_total_depth + 1e-6)
        res["raw_bin_CER"] = res["total_cis_fragments"] / (raw_total_depth + 1e-6)
        
        summary_list.append(res)
        
    return pd.DataFrame(summary_list)

import numpy as np
import pandas as pd

def recompute_bin_summary_with_frag_count(position_df, hetero_range=(0.4, 0.6)):
    """
    [UPDATED]
    각 Bin 영역 내에 존재하는 순수 고유 프래그먼트 총합(raw_bin_fragments)을 추가로 집계하여,
    후속 프래그먼트 기반 정규화(Normalization)의 완벽한 분모 뼈대를 구축합니다.
    """
    if position_df is None or position_df.empty:
        return pd.DataFrame()
    
    df = position_df.copy()
    
    h_min, h_max = hetero_range
    df['group'] = 'Other'
    df.loc[(df['pop_af'] >= h_min) & (df['pop_af'] <= h_max), 'group'] = 'Hetero'
    df.loc[(df['pop_af'] <= 0.1) | (df['pop_af'] >= 0.9), 'group'] = 'Homo'
    
    summary_list = []
    
    for bin_id, group in df.groupby("bin_id", sort=False):
        # 1) 원천 데이터(Depth) 합계 추출
        raw_ref_sum = int(group["ref_depth"].sum())
        raw_alt_sum = int(group["alt_depth"].sum())
        raw_total_depth = raw_ref_sum + raw_alt_sum
        
        # 2) [NEW] 이 영역 안에 물리적으로 '프래그먼트'가 총 몇 개나 지나갔는지 캡처
        # 포지션별로 찍힌 total_fragments 중 최댓값이나 평균 스케일을 고려하여 
        # 이 영역의 고유 프래그먼트 풀(Pool) 크기를 정의합니다. (BAM 파싱 원천 개수 대용)
        raw_bin_fragments = int(group["total_fragments"].max()) if not group.empty else 0
        if raw_bin_fragments == 0:
            # 안전 장치: 만약 max가 0이면 sum 스케일 기반 보정
            raw_bin_fragments = int(group["total_fragments"].median())
        
        chrom = group["chrom"].iloc[0]
        try:
            start = group["pos"].min() 
            end = group["pos"].max()
        except:
            start, end = 0, 0
            
        # 3) QC 퀄리티 필터링 적용 (포지션당 파편 2개 이상)
        qc_group = group[group['total_fragments'] >= 2]
        grouped_stats = qc_group.groupby('group')['observed_baf'].agg(['mean', 'median', 'count'])
        
        # 4) 위상(Phasing) 지원 개수 합산
        total_trans = int(group["trans_support"].sum())
        total_cis = int(group["cis_support"].sum())
        
        res = {
            "bin_id": bin_id,
            "chrom": chrom,
            "start": start,
            "end": end,
            "raw_total_depth": raw_total_depth,
            # [DATA BINDING] 후속 정규화 분모로 정밀 조율할 원천 변수
            "raw_bin_fragments": raw_bin_fragments,
            "total_trans_fragments": total_trans,
            "total_cis_fragments": total_cis
        }
        
        # BAF 데이터 매핑
        for g in ['Hetero', 'Homo']:
            if g in grouped_stats.index:
                res[f"{g.lower()}_baf_mean"] = grouped_stats.loc[g, 'mean']
                res[f"{g.lower()}_baf_median"] = grouped_stats.loc[g, 'median']
                res[f"{g.lower()}_sites_count"] = int(grouped_stats.loc[g, 'count'])
            else:
                res[f"{g.lower()}_baf_mean"] = res[f"{g.lower()}_baf_median"] = res[f"{g.lower()}_sites_count"] = 0
        
        # [NEW FORMULA] 영역 내 순수 프래그먼트 총합 대비 위상 비율로 정형화
        # 시퀀싱 리드 깊이가 아닌 '분자 풀의 크기' 대비 발생 비로 보정됩니다.
        denom = raw_bin_fragments if raw_bin_fragments > 0 else 1
        res["raw_bin_TER"] = total_trans / denom
        res["raw_bin_CER"] = total_cis / denom
        
        summary_list.append(res)
        
    return pd.DataFrame(summary_list)

#sample_set = [
#    'cbNIPT_24_04_01_DS',
#    'cbNIPT_24_04_02_DS',
#    'cbNIPT_24_04_03_DS',
#    'cbNIPT_24_04_04',
#    'cbNIPT_24_04_05',
#    'cbNIPT_24_04_06',
#    'cbNIPT_24_04_07',
#    'cbNIPT_24_04_08',
#    'cbNIPT_24_04_09',
#]


path = '/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Results/manual'

sample_set = os.listdir(path)
s = []
#print(sample_set)
for sample_id in sample_set:
    if sample_id.find('Insilico') > 0:
        s.append(sample_id)
        
sample_set = s

for sample in sample_set:

    #path_list = Path(path).glob(f'**/data/{sample}.genetic_evidence.tsv')
    path_list = Path(path).glob(f'**/baf/{sample}.position_baf_detail.tsv')
    
    for trans_fragment_path in path_list:
        DataId = str(Path(trans_fragment_path).name).split('.')[0]
        print(DataId)
        df = pd.read_csv(trans_fragment_path,sep='\t')
        #print(df)
        summary_df = recompute_bin_summary(df)
        print(summary_df)
        summary_df.to_csv(str(trans_fragment_path).replace("tsv","summary.tsv"),sep='\t',index=False)
        #print(df[['homo_sites_count','hetero_sites_count','total_trans_fragments','total_cis_fragments','chrom']].groupby('chrom').mean())
        #exit()
        #plot_baf_distribution_profile(df)
        #plot_normalized_fragments(df)