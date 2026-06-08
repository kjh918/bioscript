import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# 커스텀 정규화 모듈 임포트
from normalization import normalize_all_metrics_with_sex_log2

path = '/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Results/manual.new'

# Insilico 등 특정 샘플을 필터링하려면 이 리스트를 수정하세요.
sample_set = [
    'cbNIPT_24_04_01_DS', 'cbNIPT_24_04_02_DS', 'cbNIPT_24_04_03_DS',
    'cbNIPT_24_04_04', 'cbNIPT_24_04_05', 'cbNIPT_24_04_06',
    'cbNIPT_24_04_07', 'cbNIPT_24_04_08', 'cbNIPT_24_04_09',
]
sample_summary_accumulator = []

for sample in sample_set:
    path_list = Path(path).glob(f'**/baf/{sample}.position_baf_detail.tsv')
    
    for trans_fragment_path in path_list:
        DataId = str(Path(trans_fragment_path).name).split('.')[0]
        summary_path = str(trans_fragment_path).replace(".position_baf_detail.tsv", ".summary.tsv")
        print(f"[Process] Loading: {summary_path}")
        
        if not os.path.exists(summary_path):
            summary_path = str(trans_fragment_path).replace(".position_baf_detail.tsv", ".position_baf_detail.summary.tsv")
            
        try:
            summary_df = pd.read_csv(summary_path, sep='\t')
        except FileNotFoundError:
            continue
        
        try:
            # 1. 정규화 수행
            norm_df = normalize_all_metrics_with_sex_log2(
                df=summary_df, 
                depth_col="raw_total_depth", 
                sex_threshold=0.0001
            )
            
            # 2. [핵심 변경] chrom 단위가 아닌 bin_id 단위로 데이터 추출
            # Bin ID 원본 순서를 유지하기 위해 복사
            bin_data = norm_df[['bin_id', 'chrom', 'hetero_sites_count', 'homo_sites_count']].copy()
            bin_data['Sample_ID'] = DataId
            
            sample_summary_accumulator.append(bin_data)
            
        except Exception as e:
            print(f"[-] Error processing {DataId}: {e}")
            continue

# -----------------------------------------------------------------
# [BIN-LEVEL EXPECTED DELTA NORMALIZATION ENGINE]
# -----------------------------------------------------------------
if sample_summary_accumulator:
    master_df = pd.concat(sample_summary_accumulator, ignore_index=True)

    print("\n" + "="*80)
    print("[★] Calculating Bin-level Hetero & Homo Expected Delta (Log2FC)...")
    print("="*80)
    
    # 1. 샘플별 앵커 염색체(chr1~5) 총합 산출
    ANCHOR_CHROMS = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5'] 
    anchor_df = master_df[master_df['chrom'].isin(ANCHOR_CHROMS)]
    
    sample_anchors = anchor_df.groupby('Sample_ID')[['hetero_sites_count', 'homo_sites_count']].sum().reset_index()
    sample_anchors.rename(columns={'hetero_sites_count': 'Anchor_Hetero', 'homo_sites_count': 'Anchor_Homo'}, inplace=True)
    master_df = master_df.merge(sample_anchors, on='Sample_ID', how='left')
    
    # 2. 개별 Bin의 앵커 대비 '상대 비율' 연산
    master_df['Ratio_Hetero'] = master_df['hetero_sites_count'] / master_df['Anchor_Hetero']
    master_df['Ratio_Homo'] = master_df['homo_sites_count'] / master_df['Anchor_Homo']
    
    # 3. [핵심] 전체 샘플 기준, 각 "Bin"별 중앙값(Median)을 Expected Ratio로 고정
    # 특정 Bin(예: chr21_1000000)이 평소에 얼마나 읽히는지를 기준점으로 잡음
    expected_df = master_df[~master_df['chrom'].isin(['chrX', 'chrY'])].groupby('bin_id')[['Ratio_Hetero', 'Ratio_Homo']].median().reset_index()
    expected_df.rename(columns={'Ratio_Hetero': 'Expected_Hetero', 'Ratio_Homo': 'Expected_Homo'}, inplace=True)
    master_df = master_df.merge(expected_df, on='bin_id', how='left')
    
    # 4. 차이값 연산 (Log2 관측치/기대치)
    master_df['Hetero_Delta_Log2'] = np.log2((master_df['Ratio_Hetero'] + 1e-9) / (master_df['Expected_Hetero'] + 1e-9))
    master_df['Homo_Delta_Log2'] = np.log2((master_df['Ratio_Homo'] + 1e-9) / (master_df['Expected_Homo'] + 1e-9))
    
    # 성염색체 제외
    master_df = master_df[~master_df['chrom'].isin(['chrX', 'chrY'])]
    
    # 5. 매트릭스 피벗 변환 (Row: bin_id, Col: Sample_ID)
    hetero_matrix = master_df.pivot(index='bin_id', columns='Sample_ID', values='Hetero_Delta_Log2')
    homo_matrix = master_df.pivot(index='bin_id', columns='Sample_ID', values='Homo_Delta_Log2')
    
    # [중요] Bin ID를 원래 유전체 물리적 순서대로 재정렬
    bin_order_map = {b: i for i, b in enumerate(master_df['bin_id'].unique())}
    hetero_matrix['order'] = hetero_matrix.index.map(bin_order_map)
    hetero_matrix = hetero_matrix.sort_values('order').drop(columns=['order'])
    
    homo_matrix['order'] = homo_matrix.index.map(bin_order_map)
    homo_matrix = homo_matrix.sort_values('order').drop(columns=['order'])

    # Y축 염색체 경계선 계산 (어디서 염색체가 바뀌는지 탐지)
    bin_chrom_map = master_df.drop_duplicates('bin_id').set_index('bin_id')['chrom']
    bin_chrom_ordered = bin_chrom_map.loc[hetero_matrix.index]
    
    # 염색체가 바뀌는 인덱스 추출
    chrom_changes = bin_chrom_ordered != bin_chrom_ordered.shift()
    boundary_indices = np.where(chrom_changes)[0]
    boundary_labels = bin_chrom_ordered.iloc[boundary_indices].values

    # -----------------------------------------------------------------
    # [VISUALIZATION] 고해상도 Bin 단위 히트맵 렌더링
    # -----------------------------------------------------------------
    print("\n[★] Generating High-Resolution Bin-level Heatmaps...")
    
    # Y축 데이터가 엄청나게 길어지므로 세로로 긴 Figure 사용
    fig, axes = plt.subplots(1, 2, figsize=(18, 14), sharey=True)
    
    def draw_high_res_heatmap(ax, matrix, title):
        # 텍스트 없이 색상으로만 매핑 (-1.5 ~ 1.5 스케일 컷오프)
        cax = ax.imshow(matrix.values, cmap='RdBu_r', aspect='auto', vmin=-1.5, vmax=1.5)
        
        ax.set_xticks(np.arange(len(matrix.columns)))
        ax.set_xticklabels(matrix.columns, rotation=45, ha='right', fontsize=11, fontweight='bold')
        
        # 염색체 경계 가이드라인 투사
        for idx in boundary_indices:
            ax.axhline(idx, color='black', linewidth=0.8, alpha=0.6)
            
        # Y축 라벨은 Bin 이름 대신 염색체 이름으로 심플하게 배치
        ax.set_yticks(boundary_indices)
        ax.set_yticklabels([str(lbl).replace('chr', '') for lbl in boundary_labels], fontsize=10, fontweight='bold')
        
        ax.set_title(title, fontsize=15, fontweight='bold', pad=15)
        return cax

    # 좌측: Hetero 차이값
    cax1 = draw_high_res_heatmap(axes[0], hetero_matrix, "Hetero Sites Difference (Log2 FC)")
    axes[0].set_ylabel("Chromosomes (Bin-Resolution)", fontsize=14, labelpad=10)
    
    # 우측: Homo 차이값
    cax2 = draw_high_res_heatmap(axes[1], homo_matrix, "Homo Sites Difference (Log2 FC)")
    
    # 공통 컬러바
    cbar = fig.colorbar(cax2, ax=axes, pad=0.02, fraction=0.03, aspect=40)
    cbar.set_label('Log2 (Observed / Expected)', rotation=270, labelpad=20, fontsize=12, fontweight='bold')
    
    plt.suptitle("Bin-level High Resolution Profiling: ADO vs Genomic Alteration", fontsize=20, fontweight='bold', y=0.97)
    plt.tight_layout(rect=[0, 0, 0.95, 0.95]) 
    
    output_heatmap_map = os.path.join(path, "BinLevel_Hetero_Homo_Heatmaps.png")
    plt.savefig(output_heatmap_map, dpi=300)
    plt.close()
    
    print(f"[+] Output High-Res Bin Heatmap Chart saved to: {output_heatmap_map}")
    print("="*80)

else:
    print("\n[-] Failed to generate sample comparison matrix. No data parsed.")