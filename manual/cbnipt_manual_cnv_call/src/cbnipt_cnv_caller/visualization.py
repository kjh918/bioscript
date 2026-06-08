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

def plot_genome_wide_cnv_with_gaps(bins_df, segments_df, output_path, run_id):
    """
    모든 염색체를 가로로 연속해서 그리되, 인덱스 개수가 아닌 '실제 염색체 크기(bp)'를 
    X축으로 사용하여 데이터가 없는 빈 공간(Centromere 등)이 차트에 공백으로 그대로 노출됩니다.
    """
    if bins_df is None or segments_df is None or output_path is None or run_id is None:
        raise ValueError("All inputs (bins_df, segments_df, output_path, run_id) are mandatory.")

    plt.figure(figsize=(25, 6)) # 빈 공간 표현을 위해 가로를 조금 더 확장
    
    current_genome_offset = 0
    chrom_ticks = []
    chrom_names = []
    
    # 1. 자연스러운 염색체 순서 정렬 보장 (chr1~22, X, Y)
    def chrom_key(c):
        c_str = str(c).replace('chr', '')
        if c_str == 'X': return 23
        if c_str == 'Y': return 24
        try: return int(c_str)
        except: return 99

    unique_chroms = sorted(bins_df['chrom'].unique(), key=chrom_key)

    # 2. 염색체별로 루프 돌며 절대적인 유전체 좌표계(Bp scale)로 매핑
    for chrom in unique_chroms:
        group = bins_df[bins_df["chrom"] == chrom].sort_values("start")
        if group.empty: continue
            
        # [CRITICAL OPTIMIZATION] 해당 염색체의 실제 시작과 끝 물리 좌표 획득
        # 이를 통해 데이터가 듬성듬성해도 실제 유전체 크기만큼의 공간을 확보함
        chrom_max_end = group["end"].max()
        
        # Bin 산점도 플로팅 (X축: 누적 오프셋 + 실제 Bin의 중간 좌표)
        bin_mids = (group["start"].values + group["end"].values) / 2
        abs_x_positions = current_genome_offset + bin_mids
        
        plt.scatter(abs_x_positions, group["log2_chrom_norm"], 
                    s=1.5, c='lightgrey', alpha=0.5, edgecolors='none')
        
        # 해당 염색체의 세그먼트 오버레이
        chrom_segs = segments_df[segments_df["chrom"] == chrom]
        
        for _, seg in chrom_segs.iterrows():
            # 세그먼트 역시 인덱스가 아닌 실제 bp 좌표를 누적 오프셋에 더해 절대 좌표로 변환
            seg_start_x = current_genome_offset + seg["start"]
            seg_end_x = current_genome_offset + seg["end"]
            
            # 카피 넘버별 컬러셋
            color = 'black'
            if seg["copy_number"] < 2: color = 'blue'
            elif seg["copy_number"] > 2: color = 'red'
            
            plt.hlines(y=seg["seg_mean"], xmin=seg_start_x, xmax=seg_end_x, 
                       colors=color, linewidth=3.5, zorder=3)

        # Tick 위치 계산 (해당 염색체의 물리적 중심점)
        chrom_ticks.append(current_genome_offset + (chrom_max_end / 2))
        chrom_names.append(chrom.replace("chr", ""))
        
        # [KEY CHANGES] 다음 염색체로 넘어가기 전, 실제 염색체의 Max bp 크기만큼 오프셋을 전진시킴
        # 이 처리를 통해 매핑 리드가 전혀 없는 Centromere 등은 자동으로 '빈 공간' 공백 레이아웃이 됨
        current_genome_offset += chrom_max_end
        
        # 염색체 경계 구분선
        plt.axvline(x=current_genome_offset, color='black', linestyle='--', linewidth=0.5, alpha=0.7)

    # 차트 스타일링 가이드 세팅
    plt.xticks(chrom_ticks, chrom_names, fontsize=10)
    plt.xlim(0, current_genome_offset)
    plt.ylim(-3, 3) 
    plt.title(f"Genome-wide CNV Profile (Physical Coordinate with Gaps): {run_id}", fontsize=14, fontweight='bold')
    plt.ylabel("Log2 Ratio (Normalized)", fontsize=12)
    plt.xlabel("Chromosome (Physical Position)", fontsize=12)
    plt.grid(True, axis='y', linestyle=':', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    log(f"True physical coordinate genome-wide plot saved to: {output_path}")

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from utils import log

def plot_genome_wide_cnv_with_annotations(bins_df, segments_df, output_path, run_id):
    """
    데이터 존재 여부와 상관없이 chr1~22, X, Y 전체 라벨과 공간을 무조건 고정하여 생성합니다.
    시퀀싱 리드가 없는 염색체(예: 여성의 chrY)는 깨끗한 공백(Gaps)으로 그래프에 표현됩니다.
    """
    if bins_df is None or segments_df is None or output_path is None or run_id is None:
        raise ValueError("All inputs are mandatory.")

    # [DATA] hg38 기준 대표적 염색체별 Centromere(p/q 경계) 대략적 좌표 (bp)
    CENTROMERES = {
        'chr1': 123400000, 'chr2': 93900000, 'chr3': 90900000, 'chr4': 50000000,
        'chr5': 48800000,  'chr6': 59100000, 'chr7': 60100000, 'chr8': 45200000,
        'chr9': 43000000,  'chr10': 39800000, 'chr11': 53400000, 'chr12': 35800000,
        'chr13': 17700000, 'chr14': 17200000, 'chr15': 19000000, 'chr16': 36800000,
        'chr17': 25100000, 'chr18': 18500000, 'chr19': 26200000, 'chr20': 28100000,
        'chr21': 12000000, 'chr22': 15000000, 'chrX': 61000000,  'chrY': 10400000
    }

    # [DATA] hg38 각 염색체별 표준 최대 크기 (bp) -> 데이터가 없어도 이 크기만큼 공간을 확보함
    CHROM_MAX_SIZES = {
        'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555,
        'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636,
        'chr9': 138394717, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
        'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345,
        'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
        'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415
    }

    # [DATA] 주요 염색체 이상 질환 Annotation 정보
    DISEASE_REGIONS = {
        'chr13': {'start': 0, 'end': 114364328, 'name': 'Patau Synd. (13q)', 'color': '#ff7f0e'},
        'chr18': {'start': 0, 'end': 80373285, 'name': 'Edwards Synd. (18q)', 'color': '#9467bd'},
        'chr21': {'start': 0, 'end': 46709983, 'name': 'Down Synd. (21q)', 'color': '#d62728'}
    }

    fig, ax = plt.subplots(figsize=(26, 7))
    
    current_genome_offset = 0
    chrom_ticks = []
    chrom_names = []
    
    # [FIXED] 데이터프레임의 unique 대신, 인간 표준 표준 24개 염색체 리스트를 순서대로 강제 지정
    ALL_CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ['chrX', 'chrY']

    # 게놈 전역 표준 루프 시작
    for chrom in ALL_CHROMOSOMES:
        # [FIXED] 표준 크기 사전에서 해당 염색체의 길이를 가져옴 (데이터가 없어도 최소 공간 유지)
        chrom_max_end = CHROM_MAX_SIZES.get(chrom, 100000000) 
        
        # 실제 데이터프레임에서 해당 염색체 추출
        group = bins_df[bins_df["chrom"] == chrom].sort_values("start")
        
        # 데이터가 있는 경우에만 산점도를 플로팅 (없으면 통과하여 빈 공간으로 남음)
        if not group.empty:
            bin_mids = (group["start"].values + group["end"].values) / 2
            abs_x_positions = current_genome_offset + bin_mids
            
            # 1. 기저 데이터 산점도 (Log2 ratio)
            ax.scatter(abs_x_positions, group["log2_chrom_norm"], 
                       s=1.2, c='lightgrey', alpha=0.4, edgecolors='none', zorder=1)
            
            # 실제 데이터가 있다면 데이터의 max 값으로 정밀 보정
            chrom_max_end = max(chrom_max_end, group["end"].max())

        # 2. p-arm / q-arm 경계 (Centromere) 시각화
        if chrom in CENTROMERES:
            centro_bp = CENTROMERES[chrom]
            if centro_bp < chrom_max_end:
                abs_centro_x = current_genome_offset + centro_bp
                ax.axvline(x=abs_centro_x, color='blue', linestyle=':', linewidth=0.6, alpha=0.3)
                
                if chrom in ['chr1', 'chr13', 'chr18', 'chr21']:
                    ax.text(abs_centro_x - (centro_bp/2), 2.7, 'p', 
                            fontsize=9, color='blue', alpha=0.4, ha='center', va='center')
                    ax.text(abs_centro_x + ((chrom_max_end - centro_bp)/2), 2.7, 'q', 
                            fontsize=9, color='blue', alpha=0.4, ha='center', va='center')

        # 3. 질병 연관 영역 (Disease ROI) 주석
        if chrom in DISEASE_REGIONS:
            disease = DISEASE_REGIONS[chrom]
            dis_start = current_genome_offset + disease['start']
            dis_end = current_genome_offset + disease['end']
            ax.axvspan(dis_start, dis_end, color=disease['color'], alpha=0.08, zorder=0)
            
            text_x = (dis_start + dis_end) / 2
            ax.text(text_x, 2.3, disease['name'], fontsize=10, color=disease['color'],
                    fontweight='bold', ha='center', bbox=dict(facecolor='white', alpha=0.8, edgecolor=disease['color'], boxstyle='round,pad=0.3'))

        # 4. 세그먼트 데이터 오버레이 (데이터가 있을 때만 루프 수행)
        chrom_segs = segments_df[segments_df["chrom"] == chrom]
        for _, seg in chrom_segs.iterrows():
            seg_start_x = current_genome_offset + seg["start"]
            seg_end_x = current_genome_offset + seg["end"]
            
            color = 'black'
            if seg["copy_number"] < 2: color = 'blue'
            elif seg["copy_number"] > 2: color = 'red'
            
            ax.hlines(y=seg["seg_mean"], xmin=seg_start_x, xmax=seg_end_x, 
                      colors=color, linewidth=3.5, zorder=4)

        # Tick & 가이드라인 업데이트 (중앙 정렬 좌표 바인딩)
        chrom_ticks.append(current_genome_offset + (chrom_max_end / 2))
        chrom_names.append(chrom.replace("chr", ""))
        
        # [CRITICAL] 데이터 유무에 상관없이 무조건 표준 염색체 크기만큼 오프셋 축을 전진시킴
        current_genome_offset += chrom_max_end
        
        # 염색체 대경계선 (Solid Black)
        ax.axvline(x=current_genome_offset, color='black', linestyle='--', linewidth=0.5, alpha=0.6)

    # 5. 최종 레이아웃 데코레이션
    ax.set_xticks(chrom_ticks)
    ax.set_xticklabels(chrom_names, fontsize=11, fontweight='bold')
    ax.set_xlim(0, current_genome_offset)
    ax.set_ylim(-2.5, 3.0)
    ax.axhline(0, color='grey', linestyle='-', linewidth=0.8, alpha=0.5)
    
    ax.set_title(f"Annotated Genome-wide CNV Profile: {run_id} (Fixed Standard Chromosome Axes)", fontsize=15, fontweight='bold', pad=15)
    ax.set_ylabel("Log2 Ratio (Normalized)", fontsize=12)
    ax.set_xlabel("Chromosomes (p/q Arms Split by Centromere Dotted Lines)", fontsize=12)
    ax.grid(True, axis='y', linestyle=':', alpha=0.4)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    log(f"Fixed standard coordinate genome-wide plot saved to: {output_path}")