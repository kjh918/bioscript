import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from utils import log

def plot_gc_correction(gc_stats, out_path):
    if gc_stats is None: return
    x, y, fit = gc_stats
    plt.figure(figsize=(8, 5))
    plt.scatter(x, y, s=1, color='gray', alpha=0.2, label="Raw")
    idx = x.argsort()
    plt.plot(x[idx], fit[idx], color='red', lw=2, label="Lowess Fit")
    plt.title("GC Bias Correction (LOWESS)"); plt.legend(); plt.savefig(out_path); plt.close()

def plot_genome_profile(df, col, title, out_path):
    plt.figure(figsize=(18, 6))
    colors = ['#4E79A7', '#A0CBE8']
    ticks, labels = [], []
    df_plot = df.dropna(subset=[col]).reset_index(drop=True)
    
    for i, (chrom, g) in enumerate(df_plot.groupby("chrom", sort=False)):
        plt.scatter(g.index, g[col], s=2, color=colors[i%2], alpha=0.6)
        ticks.append(g.index.mean()); labels.append(chrom.replace("chr", ""))
        plt.axvline(g.index.max(), color='black', alpha=0.1, lw=0.5)
        
    plt.axhline(0, color='red', lw=1, linestyle="--")
    plt.xticks(ticks, labels, fontsize=8); plt.ylim(-2, 2)
    plt.title(title); plt.tight_layout(); plt.savefig(out_path); plt.close()

def plot_chrom_summary(chr_summary, out_path):
    """Aneuploidy 확인용 Z-score 막대 그래프"""
    if chr_summary is None or chr_summary.empty:
        print("[WARN] chr_summary is empty. Skipping plot.")
        return

    plt.figure(figsize=(12, 6))
    sns.barplot(data=chr_summary, x="chrom", y="z_score", palette="coolwarm")
    plt.axhline(3, color='red', linestyle='--', label="Trisomy (Z>3)")
    plt.axhline(-3, color='blue', linestyle='--', label="Monosomy (Z<-3)")
    plt.title("Chromosome-level Z-score (Aneuploidy Check)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()
def plot_genome_wide_cnv(bins_df, segments_df, output_path, run_id):
    # [CHECK] 필수 데이터 확인 (No default)
    if bins_df is None:
        raise ValueError("bins_df is mandatory for plotting.")
    if segments_df is None:
        raise ValueError("segments_df is mandatory for plotting.")
    if output_path is None:
        raise ValueError("output_path is mandatory.")
    if run_id is None:
        raise ValueError("run_id is mandatory.")

    plt.figure(figsize=(20, 6))
    
    current_pos = 0
    chrom_ticks = []
    chrom_names = []
    
    # [MODIFIED] 데이터의 염색체 순서 보장 및 상대 좌표 계산
    for chrom, group in bins_df.groupby("chrom", sort=False):
        n_bins = len(group)
        
        # [MODIFIED] group.index 대신 0부터 시작하는 상대 좌표(rel_idx) 생성
        # 이렇게 해야 current_pos와 합쳐졌을 때 연속적인 X축이 보장됨
        rel_idx = np.arange(n_bins)
        abs_idx = rel_idx + current_pos
        
        # Bin 산점도 (Log2 ratio)
        plt.scatter(abs_idx, group["log2_chrom_norm"], 
                    s=1, c='lightgrey', alpha=0.5)
        
        # 해당 염색체의 세그먼트 오버레이
        chrom_segs = segments_df[segments_df["chrom"] == chrom]
        
        # [MODIFIED] 세그먼트 매핑 시에도 상대 인덱스 기반으로 검색하여 정렬 일치
        bin_starts = group["start"].values
        bin_ends = group["end"].values
        
        for _, seg in chrom_segs.iterrows():
            # [MODIFIED] np.searchsorted를 사용하여 실제 bin 위치와 정확히 매칭 (O(log n))
            start_rel_pos = np.searchsorted(bin_starts, seg["start"])
            end_rel_pos = np.searchsorted(bin_ends, seg["end"], side='right') - 1
            
            # 인덱스 범위 초과 방지
            start_rel_pos = max(0, min(start_rel_pos, n_bins - 1))
            end_rel_pos = max(0, min(end_rel_pos, n_bins - 1))

            seg_start_x = start_rel_pos + current_pos
            seg_end_x = end_rel_pos + current_pos
            
            color = 'black'
            if seg["copy_number"] < 2: color = 'blue'
            elif seg["copy_number"] > 2: color = 'red'
            
            plt.hlines(y=seg["seg_mean"], xmin=seg_start_x, xmax=seg_end_x, 
                       colors=color, linewidth=3)

        # Tick 위치 계산 (해당 염색체 중앙)
        chrom_ticks.append(current_pos + n_bins // 2)
        chrom_names.append(chrom.replace("chr", ""))
        
        # [MODIFIED] 구분선 및 위치 업데이트
        current_pos += n_bins
        plt.axvline(x=current_pos, color='black', linestyle='--', linewidth=0.5)

    plt.xticks(chrom_ticks, chrom_names)
    plt.ylim(-3, 3) 
    plt.title(f"Genome-wide CNV Profile: {run_id}")
    plt.ylabel("Log2 Ratio (Normalized)")
    plt.xlabel("Chromosome")
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    log(f"Genome-wide plot saved to: {output_path}")