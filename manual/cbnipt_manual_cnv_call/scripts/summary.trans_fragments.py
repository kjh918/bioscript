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


def recompute_bin_summary_with_frag_norm(position_df, hetero_range=(0.4, 0.6)):
    """
    포지션 레벨 상세 데이터에서 영역별 원천 total_fragments 및 QC 패스 fragments를 정확히 카운트하고,
    이를 분모로 삼아 프래그먼트 기반의 TER/CER 정규화 지표를 산출합니다.
    """
    if position_df is None or position_df.empty:
        return pd.DataFrame()
    
    # 1. 안전한 데이터 처리를 위해 카피 
    df = position_df.copy()
    
    summary_list = []
    
    # 2. 각 Bin ID별로 그룹화하여 통계량 재집계
    for bin_id, group in df.groupby("bin_id", sort=False):
        # 1) 필터링 전 영역 내 원천 데이터(Depth & Total Fragments) 합계 추출
        raw_ref_sum = int(group["ref_depth"].sum())
        raw_alt_sum = int(group["alt_depth"].sum())
        raw_total_depth = raw_ref_sum + raw_alt_sum
        
        # [NEW] 이 영역의 원천 total_fragments 전체 포지션 합산 값
        raw_total_fragments_sum = int(group["total_fragments"].sum())
        
        # 2) 유전체 매핑 위치 복원
        chrom = group["chrom"].iloc[0]
        try:
            start = group["pos"].min() 
            end = group["pos"].max()
        except:
            start, end = 0, 0
            
        # 3) QC 퀄리티 필터링 적용 (포지션당 읽힌 파편 수가 2개 이상인 클린 지점 격리)
        qc_mask = group['total_fragments'] >= 2
        qc_group = group[qc_mask]
        
        # [NEW] QC를 통과한 알짜배기 고품질 프래그먼트 총합 격리 카운트
        qc_pass_fragments = int(qc_group["total_fragments"].sum()) if not qc_group.empty else 0
        
        # QC 통과 그룹 내에서만 안정적인 BAF 요약치 추출
        summary_grouped = qc_group.groupby('group')['observed_baf'].agg(['mean', 'median', 'count']) if not qc_group.empty else pd.DataFrame()
        
        # 4) 위상(Phasing) 지원 프래그먼트 개수 단순 합산
        total_trans = int(group["trans_support"].sum())
        total_cis = int(group["cis_support"].sum())
        
        res = {
            "bin_id": bin_id,
            "chrom": chrom,
            "start": start,
            "end": end,
            # 원천 데이터 자산 보존
            "raw_ref_depth_sum": raw_ref_sum,
            "raw_alt_depth_sum": raw_alt_sum,
            "raw_total_depth": raw_total_depth,
            "raw_total_fragments_sum": raw_total_fragments_sum,
            # QC 패스된 청정 프래그먼트 풀 데이터 자산 격리
            "qc_pass_fragments": qc_pass_fragments,
            "total_trans_fragments": total_trans,
            "total_cis_fragments": total_cis
        }
        
        # Hetero, Homo 각각의 BAF 분포 및 마커 개수 매핑
        for g in ['Hetero', 'Homo']:
            if not summary_grouped.empty and g in summary_grouped.index:
                res[f"{g.lower()}_baf_mean"] = summary_grouped.loc[g, 'mean']
                res[f"{g.lower()}_baf_median"] = summary_grouped.loc[g, 'median']
                res[f"{g.lower()}_sites_count"] = int(summary_grouped.loc[g, 'count'])
            else:
                res[f"{g.lower()}_baf_mean"] = res[f"{g.lower()}_baf_median"] = res[f"{g.lower()}_sites_count"] = 0
                
        # -----------------------------------------------------------------
        # [NEW NORMALIZATION FORMULA] 영역별 total_fragments 기반 스케일링
        # -----------------------------------------------------------------
        # 분모 오염과 왜곡을 방어하기 위해, QC를 완전히 통과하여 살아남은 
        # 순수 유효 프래그먼트 풀(qc_pass_fragments)을 기준으로 밀도 분율을 최종 결정합니다.
        denom = qc_pass_fragments if qc_pass_fragments > 0 else 1
        res["raw_bin_TER"] = total_trans / denom
        res["raw_bin_CER"] = total_cis / denom
        
        summary_list.append(res)
        
    return pd.DataFrame(summary_list)

import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed


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


path = '/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Results/manual.new'

sample_set = os.listdir(path)
s = []
#print(sample_set)
for sample_id in sample_set:
    if sample_id.find('Insilico') > 0:
        s.append(sample_id)
        
sample_set = s
def process_single_file(file_path):
    """
    각 CPU 코어(Process)에 할당되어 독립적으로 실행될 워커 함수입니다.
    """
    DataId = str(Path(file_path).name).split('.')[0]
    
    try:
        # 1. 데이터 로드
        df = pd.read_csv(file_path, sep='\t')
        
        # 2. 무거운 병목 연산 (recompute_bin_summary_with_frag_norm)
        summary_df = recompute_bin_summary_with_frag_norm(df)
        
        # 3. 결과 저장
        out_path = str(file_path).replace(".position_baf_detail.tsv", ".summary.tsv")
        summary_df.to_csv(out_path, sep='\t', index=False)
        
        return f"[Success] {DataId} (Shape: {summary_df.shape})"
        
    except Exception as e:
        return f"[Error] Failed on {DataId}: {e}"

def main():
    path = '/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Results/manual.new'
    
    sample_set = [
        'cbNIPT_24_04_01_DS', 'cbNIPT_24_04_02_DS', 'cbNIPT_24_04_03_DS',
        'cbNIPT_24_04_04', 'cbNIPT_24_04_05', 'cbNIPT_24_04_06',
        'cbNIPT_24_04_07', 'cbNIPT_24_04_08', 'cbNIPT_24_04_09',
    ]

    # 1. 처리해야 할 전체 타겟 파일 목록을 먼저 수집합니다.
    target_files = []
    for sample in sample_set:
        path_list = Path(path).glob(f'**/baf/{sample}.position_baf_detail.tsv')
        target_files.extend(list(path_list))
        
    print(f"[*] Total files to process: {len(target_files)}")
    
    if not target_files:
        print("[-] No files found. Check your path and sample_set.")
        return

    # 2. 사용 가능한 CPU 코어 수 확인 (안정성을 위해 최대 코어-1 권장)
    max_workers = min(os.cpu_count() - 1, len(target_files))
    if max_workers < 1: max_workers = 1
    print(f"[*] Firing up {max_workers} processes...\n")

    # 3. ProcessPoolExecutor를 이용한 병렬 연산 가동
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # 각 파일을 개별 프로세스에 매핑하여 비동기 실행 예약
        futures = {executor.submit(process_single_file, file_path): file_path for file_path in target_files}
        
        # 완료되는 순서대로 결과를 반환받아 출력 (as_completed)
        for future in as_completed(futures):
            # 워커 함수의 return 값(Success/Error 문자열)을 출력합니다.
            result_msg = future.result()
            print(result_msg)
            
    print("\n[*] All parallel processing completed successfully.")

# 멀티프로세싱을 사용할 때는 이 구문이 필수입니다. (운영체제 안전성 확보)
if __name__ == '__main__':
    main()