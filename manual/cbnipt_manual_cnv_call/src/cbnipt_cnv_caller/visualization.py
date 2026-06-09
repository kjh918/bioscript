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

    
def plot_comprehensive_baf_logr_qc(bins_df, segments_df, output_path, run_id):
    """
    [3-Panel Comprehensive View]
    Top: Copy Number (Log2 Ratio)
    Middle: B-Allele Frequency (BAF)
    Bottom: Allelic Imbalance & Hetero Rates (QC)
    모든 패널은 실제 물리적 유전체 좌표(Physical position)를 공유합니다.
    """
    if bins_df is None or segments_df is None or output_path is None or run_id is None:
        raise ValueError("All inputs are mandatory.")

    # [DATA] 인간 표준 24개 염색체 리스트 및 표준 최대 크기 (bp)
    ALL_CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ['chrX', 'chrY']
    CHROM_MAX_SIZES = {
        'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555,
        'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636,
        'chr9': 138394717, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
        'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345,
        'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
        'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415
    }

    # 3개의 행(Panel)을 가진 Figure 생성 (높이 비율 2:2:1)
    fig, axes = plt.subplots(3, 1, figsize=(28, 14), sharex=True, gridspec_kw={'height_ratios': [2, 2, 1]})
    ax_logr, ax_baf, ax_qc = axes

    current_genome_offset = 0
    chrom_ticks = []
    chrom_names = []

    for chrom in ALL_CHROMOSOMES:
        chrom_max_end = CHROM_MAX_SIZES.get(chrom, 100000000) 
        group = bins_df[bins_df["chrom"] == chrom].sort_values("start")
        
        if not group.empty:
            bin_mids = (group["start"].values + group["end"].values) / 2
            abs_x_positions = current_genome_offset + bin_mids
            
            # [Panel 1] Log2 Ratio (Depth)
            ax_logr.scatter(abs_x_positions, group["log2_chrom_norm"], 
                            s=1.5, c='lightgrey', alpha=0.6, edgecolors='none', zorder=1)
            
            # [Panel 2] BAF (Variant)
            # BAF는 0.5를 기준으로 대칭이 되도록 가이드라인을 그립니다.
            ax_baf.scatter(abs_x_positions, group["bin_BAF"], 
                           s=2.0, c='purple', alpha=0.5, edgecolors='none')
            
            # [Panel 3] QC & Imbalance Rates (0~1 scale)
            # 신뢰할 수 있는 구간(total_sites >= 5)만 라인으로 연결
            qc_valid = group[group['total_sites'] >= 5]
            if not qc_valid.empty:
                abs_x_qc = current_genome_offset + (qc_valid["start"].values + qc_valid["end"].values) / 2
                ax_qc.plot(abs_x_qc, qc_valid["imbalance_rate"], color='red', alpha=0.7, linewidth=1.5, label='Imbalance Rate' if chrom=='chr1' else "")
                ax_qc.plot(abs_x_qc, qc_valid["hetero_like_rate"], color='green', alpha=0.7, linewidth=1.5, label='Hetero Rate' if chrom=='chr1' else "")

            chrom_max_end = max(chrom_max_end, group["end"].max())

        # [Panel 1] 세그먼트 오버레이
        chrom_segs = segments_df[segments_df["chrom"] == chrom]
        for _, seg in chrom_segs.iterrows():
            seg_start_x = current_genome_offset + seg["start"]
            seg_end_x = current_genome_offset + seg["end"]
            
            color = 'black'
            if seg["copy_number"] < 2: color = 'blue'
            elif seg["copy_number"] > 2: color = 'red'
            
            ax_logr.hlines(y=seg["seg_mean"], xmin=seg_start_x, xmax=seg_end_x, colors=color, linewidth=4, zorder=4)

        # Tick & 가이드라인 설정
        chrom_ticks.append(current_genome_offset + (chrom_max_end / 2))
        chrom_names.append(chrom.replace("chr", ""))
        current_genome_offset += chrom_max_end
        
        # 각 패널에 염색체 경계선 추가
        for ax in axes:
            ax.axvline(x=current_genome_offset, color='black', linestyle='--', linewidth=0.5, alpha=0.6)

    # --- [Panel 1] 데코레이션 (Log2 Ratio) ---
    ax_logr.set_ylim(-2.5, 3.0)
    ax_logr.axhline(0, color='grey', linestyle='-', linewidth=1.0, alpha=0.5)
    ax_logr.axhline(0.58, color='red', linestyle=':', linewidth=1.0, alpha=0.5) # Duplication expected log2
    ax_logr.axhline(-1.0, color='blue', linestyle=':', linewidth=1.0, alpha=0.5) # Deletion expected log2
    ax_logr.set_ylabel("Log2 Ratio\n(Normalized Depth)", fontsize=12)
    ax_logr.set_title(f"Comprehensive Multi-omic CNV Profile: {run_id}", fontsize=16, fontweight='bold')
    
    # --- [Panel 2] 데코레이션 (BAF) ---
    ax_baf.set_ylim(-0.05, 1.05)
    ax_baf.axhline(0.5, color='gray', linestyle='-', linewidth=1.0, alpha=0.5) # Normal Hetero
    ax_baf.axhline(0.33, color='orange', linestyle=':', linewidth=1.0, alpha=0.5) # Trisomy (3n) BAF
    ax_baf.axhline(0.66, color='orange', linestyle=':', linewidth=1.0, alpha=0.5) # Trisomy (3n) BAF
    ax_baf.set_ylabel("B-Allele Frequency\n(Variant Space)", fontsize=12)
    
    # --- [Panel 3] 데코레이션 (Rates & QC) ---
    ax_qc.set_ylim(-0.05, 1.05)
    ax_qc.set_ylabel("AI / LOH Rates\n(0 to 1)", fontsize=12)
    ax_qc.set_xlabel("Chromosomes (Physical Position)", fontsize=14)
    ax_qc.legend(loc='upper right')

    # 공통 X축 데코레이션
    ax_qc.set_xticks(chrom_ticks)
    ax_qc.set_xticklabels(chrom_names, fontsize=11, fontweight='bold')
    ax_qc.set_xlim(0, current_genome_offset)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    log(f"Comprehensive 3-panel CNV plot saved to: {output_path}")