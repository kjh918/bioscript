import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
def apply_low_quality_filter(df, min_depth=1, min_coverage=0.5):
    """
    Parquet 워커에서 병합되어 나온 최종 DataFrame(df)을 필터링합니다.
    성염색체는 보호하고 상염색체만 필터링합니다.
    """
    if df is None or df.empty: return pd.DataFrame()
        
    initial_count = len(df)
    sex_chrom_mask = df["chrom"].isin(["chrX", "chrY"])
    df_autosomes = df[~sex_chrom_mask].copy()
    df_sex_chroms = df[sex_chrom_mask].copy()
    
    # Depth 및 Breadth 필터링
    df_auto_filtered = df_autosomes[
        (df_autosomes["raw_count"] >= min_depth) & 
        (df_autosomes["breadth_ratio"] >= min_coverage)
    ].copy()
    
    df_filtered = pd.concat([df_auto_filtered, df_sex_chroms], ignore_index=True)
    
    def chrom_key(c):
        c_str = str(c).replace('chr', '')
        if c_str == 'X': return 23
        if c_str == 'Y': return 24
        try: return int(c_str)
        except: return 99

    df_filtered['chrom_sort_key'] = df_filtered['chrom'].apply(chrom_key)
    df_filtered = df_filtered.sort_values(by=['chrom_sort_key', 'start']).drop(columns=['chrom_sort_key']).reset_index(drop=True)
    
    return df_filtered
# 경로 설정
path = '/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Results'

sample_set = [
    'cbNIPT_24_04_01_DS', 'cbNIPT_24_04_02_DS', 'cbNIPT_24_04_03_DS',
    'cbNIPT_24_04_04', 'cbNIPT_24_04_05', 'cbNIPT_24_04_06',
    'cbNIPT_24_04_07', 'cbNIPT_24_04_08', 'cbNIPT_24_04_09',
    'cbNIPT_24_04_10', 'cbNIPT_24_04_11', 'cbNIPT_24_04_12',
    'cbNIPT_24_04_13', 'cbNIPT_24_04_14', 'cbNIPT_24_04_15',
]

sample_summary_accumulator = []

print("\n" + "="*80)
print("[★] Calculating Leave-One-Out (LOO) Internal Log2FC...")
print("="*80)

for sample in sample_set:
    file_paths = list(Path(path).glob(f'**/data/{sample}.normalized.tsv'))
    if not file_paths: continue
        
    target_path = file_paths[0]
    
    try:
        df = pd.read_csv(target_path, sep='\t')
        
        # 1. 퀄리티 필터링
        #valid_df = df[df['total_sites'] >= 5].copy()
        valid_df = apply_low_quality_filter(df)
        #print(valid_df.columns)
        #exit()
        chrom_mean = valid_df.groupby('chrom')[['hetero_like_count', 'homo_like_count']].mean().reset_index()
        
        # -------------------------------------------------------------
        # [핵심] Leave-One-Out 중앙값 연산 함수
        # -------------------------------------------------------------
        def calc_loo_median(current_chrom, col_name):
            # 성염색체 제외 및 '자기 자신(current_chrom)'도 제외한 나머지 상염색체만 추출
            pool = chrom_mean[
                (~chrom_mean['chrom'].isin(['chrX', 'chrY'])) & 
                (chrom_mean['chrom'] != current_chrom)
            ]
            return pool[col_name].median()
            
        # 2. 각 염색체별로 자신을 제외한 베이스라인(Median)을 계산하여 컬럼에 저장
        chrom_mean['LOO_Het_Median'] = chrom_mean['chrom'].apply(lambda c: calc_loo_median(c, 'hetero_like_count'))
        chrom_mean['LOO_Hommo_Median'] = chrom_mean['chrom'].apply(lambda c: calc_loo_median(c, 'homo_like_count'))
        
        # 3. 내 염색체 값을 -> '나를 뺀 나머지 중앙값'으로 나누고 Log2 씌우기
        chrom_mean['LOO_Log2FC_Hetero'] = np.log2((chrom_mean['hetero_like_count'] + 1e-4) / (chrom_mean['LOO_Het_Median'] + 1e-4))
        chrom_mean['LOO_Log2FC_Imbal'] = np.log2((chrom_mean['homo_like_count'] + 1e-4) / (chrom_mean['LOO_Hommo_Median'] + 1e-4))
        
        chrom_mean['Sample_ID'] = sample
        sample_summary_accumulator.append(chrom_mean)
        
    except Exception as e:
        print(f"[-] Error processing {sample}: {e}")

# -----------------------------------------------------------------
# [VISUALIZATION] LOO Log2FC Heatmap
# -----------------------------------------------------------------
if sample_summary_accumulator:
    master_df = pd.concat(sample_summary_accumulator, ignore_index=True)
    
    # 성염색체 제외
    master_df = master_df[~master_df['chrom'].isin(['chrX', 'chrY'])]
    
    def chrom_key_auto(c):
        try: return int(str(c).replace('chr', ''))
        except: return 99
        
    # 매트릭스 피벗
    hetero_matrix = master_df.pivot(index='chrom', columns='Sample_ID', values='LOO_Log2FC_Hetero')
    imbal_matrix = master_df.pivot(index='chrom', columns='Sample_ID', values='LOO_Log2FC_Imbal')
    
    sorted_chroms = sorted(hetero_matrix.index, key=chrom_key_auto)
    hetero_matrix = hetero_matrix.reindex(sorted_chroms)
    imbal_matrix = imbal_matrix.reindex(sorted_chroms)

    print("\n[★] Generating Heatmaps...")
    fig, axes = plt.subplots(2, 1, figsize=(16, 14), sharex=True)
    
    def draw_heatmap(ax, matrix, title, label_text):
        # -1.0 ~ 1.0 대칭 컬러 스케일
        cax = ax.imshow(matrix.values, cmap='RdBu_r', aspect='auto', vmin=-1.0, vmax=1.0)
        
        ax.set_xticks(np.arange(len(matrix.columns)))
        ax.set_yticks(np.arange(len(matrix.index)))
        ax.set_xticklabels(matrix.columns, rotation=45, ha='right', fontsize=11, fontweight='bold')
        ax.set_yticklabels(matrix.index, fontsize=11, fontweight='bold')
        
        for i in range(len(matrix.index)):
            for j in range(len(matrix.columns)):
                val = matrix.values[i, j]
                # 가독성을 위해 색상이 진해지면 텍스트를 흰색으로 변경
                text_color = "white" if abs(val) > 0.6 else "black"
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", color=text_color, fontsize=10)
                
        cbar = fig.colorbar(cax, ax=ax, pad=0.01, fraction=0.03)
        cbar.set_label(label_text, rotation=270, labelpad=20, fontsize=11, fontweight='bold')
        ax.set_title(title, fontsize=16, fontweight='bold', pad=15)
        
        ax.set_xticks(np.arange(-.5, len(matrix.columns), 1), minor=True)
        ax.set_yticks(np.arange(-.5, len(matrix.index), 1), minor=True)
        ax.grid(which="minor", color="gray", linestyle='-', linewidth=0.5, alpha=0.5)
        ax.tick_params(which="minor", bottom=False, left=False)

    draw_heatmap(axes[0], hetero_matrix, 
                 "LOO Log2FC of Hetero Rate (Normal -> 0.0)", 
                 "Log2(Observed / LOO Median)")
                 
    draw_heatmap(axes[1], imbal_matrix, 
                 "LOO Log2FC of Imbalance Rate (Red > 0.5 implies CNV)", 
                 "Log2(Observed / LOO Median)")
    
    axes[1].set_xlabel("Samples (Evaluated 100% Independently with Leave-One-Out)", fontsize=14, labelpad=10)
    fig.text(0.04, 0.5, 'Chromosomes', va='center', rotation='vertical', fontsize=14, fontweight='bold')
    
    plt.tight_layout(rect=[0.05, 0, 1, 1]) 
    
    output_heatmap_map = os.path.join(path, "LOO_Log2FC_Heatmaps.png")
    plt.savefig(output_heatmap_map, dpi=300)
    plt.close()
    
    print(f"[+] Leave-One-Out Heatmap Chart saved to: {output_heatmap_map}")