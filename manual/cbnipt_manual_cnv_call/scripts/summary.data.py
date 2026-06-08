import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt  # 시각화를 위한 필수 라이브러리

# 커스텀 정규화 모듈 임포트
from normalization import normalize_all_metrics_with_sex_log2

path = '/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Results/manual.new'

# 모든 샘플의 정규화 결과를 하나로 모으기 위한 리스트
sample_set = [
    #'cbNIPT_24_04_01_DS',
    'cbNIPT_24_04_02_DS',
    'cbNIPT_24_04_02_DS_Ratio_1-1-Insilico_cbNIPT_24_04_04',
    'cbNIPT_24_04_02_DS_Ratio_1-1-Insilico_cbNIPT_24_04_05',
    'cbNIPT_24_04_02_DS_Ratio_1-1-Insilico_cbNIPT_24_04_06',
    'cbNIPT_24_04_02_DS_Ratio_1-2-Insilico_cbNIPT_24_04_04',
    'cbNIPT_24_04_02_DS_Ratio_1-2-Insilico_cbNIPT_24_04_05',
    'cbNIPT_24_04_02_DS_Ratio_1-2-Insilico_cbNIPT_24_04_06',
    'cbNIPT_24_04_02_DS_Ratio_2-1-Insilico_cbNIPT_24_04_04',
    'cbNIPT_24_04_02_DS_Ratio_2-1-Insilico_cbNIPT_24_04_05',
    'cbNIPT_24_04_02_DS_Ratio_2-1-Insilico_cbNIPT_24_04_06',
    #'cbNIPT_24_04_03_DS',
    #'cbNIPT_24_04_04',
    #'cbNIPT_24_04_05',
    #'cbNIPT_24_04_06',
    #'cbNIPT_24_04_07',
    #'cbNIPT_24_04_08',
    #'cbNIPT_24_04_09',
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
            norm_df = normalize_all_metrics_with_sex_log2(
                df=summary_df, 
                depth_col="raw_total_depth", 
                sex_threshold=0.0001
            )
            
            # ADO 감별을 위해 이형(Hetero) 및 동형(Homo) 마커 개수만 추출
            chrom_mean = norm_df.groupby('chrom')[['hetero_sites_count', 'homo_sites_count']].mean().reset_index()
            chrom_mean['Sample_ID'] = DataId
            sample_summary_accumulator.append(chrom_mean)
            
        except Exception as e:
            print(f"[-] Error processing {DataId}: {e}")
            continue

# -----------------------------------------------------------------
# [EXPECTED DELTA NORMALIZATION ENGINE]
# -----------------------------------------------------------------
if sample_summary_accumulator:
    master_df = pd.concat(sample_summary_accumulator, ignore_index=True)

    print("\n" + "="*80)
    print("[★] Calculating Hetero & Homo Expected Delta (Log2FC)...")
    print("="*80)
    
    # 1. 뼈대가 될 상위 5개 다중 앵커 지정
    ANCHOR_CHROMS = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5'] 
    anchor_df = master_df[master_df['chrom'].isin(ANCHOR_CHROMS)]
    
    # 샘플별 앵커 염색체 총합 산출
    sample_anchors = anchor_df.groupby('Sample_ID')[['hetero_sites_count', 'homo_sites_count']].sum().reset_index()
    sample_anchors.rename(columns={'hetero_sites_count': 'Anchor_Hetero', 'homo_sites_count': 'Anchor_Homo'}, inplace=True)
    master_df = master_df.merge(sample_anchors, on='Sample_ID', how='left')
    
    # 개별 염색체의 앵커 대비 상대 비율 연산
    master_df['Ratio_Hetero'] = master_df['hetero_sites_count'] / master_df['Anchor_Hetero']
    master_df['Ratio_Homo'] = master_df['homo_sites_count'] / master_df['Anchor_Homo']
    
    # 2. 전체 Batch 기준 중앙값(Median)을 'Expected Ratio'로 고정
    expected_df = master_df[~master_df['chrom'].isin(['chrX', 'chrY'])].groupby('chrom')[['Ratio_Hetero', 'Ratio_Homo']].median().reset_index()
    expected_df.rename(columns={'Ratio_Hetero': 'Expected_Hetero', 'Ratio_Homo': 'Expected_Homo'}, inplace=True)
    master_df = master_df.merge(expected_df, on='chrom', how='left')
    
    # 3. 차이값 연산 (Log2 관측치/기대치)
    master_df['Hetero_Delta_Log2'] = np.log2((master_df['Ratio_Hetero'] + 1e-9) / (master_df['Expected_Hetero'] + 1e-9))
    master_df['Homo_Delta_Log2'] = np.log2((master_df['Ratio_Homo'] + 1e-9) / (master_df['Expected_Homo'] + 1e-9))
    
    # 성염색체 제외 및 정렬
    master_df = master_df[~master_df['chrom'].isin(['chrX', 'chrY'])]
    
    def chrom_key_auto(c):
        c_str = str(c).replace('chr', '')
        try: return int(c_str)
        except: return 99
        
    # 매트릭스 피벗 변환
    hetero_matrix = master_df.pivot(index='chrom', columns='Sample_ID', values='Hetero_Delta_Log2')
    homo_matrix = master_df.pivot(index='chrom', columns='Sample_ID', values='Homo_Delta_Log2')
    
    sorted_chroms = sorted(hetero_matrix.index, key=chrom_key_auto)
    hetero_matrix = hetero_matrix.reindex(sorted_chroms)
    homo_matrix = homo_matrix.reindex(sorted_chroms)

    # -----------------------------------------------------------------
    # [VISUALIZATION] 2개의 히트맵(Hetero vs Homo) 상하 배치 렌더링
    # -----------------------------------------------------------------
    print("\n[★] Generating Hetero & Homo Delta Heatmaps...")
    
    fig, axes = plt.subplots(2, 1, figsize=(16, 14), sharex=True)
    
    # 공통 히트맵 플로팅 함수
    def draw_heatmap(ax, matrix, title):
        # 0을 기준으로 파란색(-)과 빨간색(+)이 대칭되도록 RdBu_r 컬러맵 사용
        cax = ax.imshow(matrix.values, cmap='RdBu_r', aspect='auto', vmin=-1.0, vmax=1.0)
        
        ax.set_xticks(np.arange(len(matrix.columns)))
        ax.set_yticks(np.arange(len(matrix.index)))
        ax.set_xticklabels(matrix.columns, rotation=45, ha='right', fontsize=11, fontweight='bold')
        ax.set_yticklabels(matrix.index, fontsize=11, fontweight='bold')
        
        # 텍스트 삽입
        for i in range(len(matrix.index)):
            for j in range(len(matrix.columns)):
                val = matrix.values[i, j]
                # 글씨 가독성을 위해 색이 진한 양극단에서는 흰색 텍스트 사용
                text_color = "white" if abs(val) > 0.6 else "black"
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", color=text_color, fontsize=10)
                
        cbar = fig.colorbar(cax, ax=ax, pad=0.01, fraction=0.03)
        cbar.set_label('Log2(Observed / Expected)', rotation=270, labelpad=20, fontsize=11, fontweight='bold')
        ax.set_title(title, fontsize=16, fontweight='bold', pad=15)
        
        # 격자선
        ax.set_xticks(np.arange(-.5, len(matrix.columns), 1), minor=True)
        ax.set_yticks(np.arange(-.5, len(matrix.index), 1), minor=True)
        ax.grid(which="minor", color="gray", linestyle='-', linewidth=0.5, alpha=0.5)
        ax.tick_params(which="minor", bottom=False, left=False)

    # 상단: Hetero 차이값 히트맵
    draw_heatmap(axes[0], hetero_matrix, "Hetero Sites Difference from Expected (Log2 FC)")
    
    # 하단: Homo 차이값 히트맵
    draw_heatmap(axes[1], homo_matrix, "Homo Sites Difference from Expected (Log2 FC)")
    
    axes[1].set_xlabel("Samples", fontsize=14, labelpad=10)
    fig.text(0.04, 0.5, 'Chromosomes', va='center', rotation='vertical', fontsize=14, fontweight='bold')
    
    plt.tight_layout(rect=[0.05, 0, 1, 1]) # 왼쪽 Y축 라벨 여백 확보
    
    output_heatmap_map = os.path.join(path, "Hetero_Homo_Delta_Heatmaps.png")
    plt.savefig(output_heatmap_map, dpi=300)
    plt.close()
    
    print(f"[+] Output Dual Heatmap Chart saved to: {output_heatmap_map}")
    print("="*80)

else:
    print("\n[-] Failed to generate sample comparison matrix. No data parsed.")