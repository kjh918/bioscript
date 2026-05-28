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
    'cbNIPT_24_04_01_DS',
    'cbNIPT_24_04_02_DS',
    #'cbNIPT_24_04_02_DS_Ratio_1-1-Insilico_cbNIPT_24_04_04',
    #'cbNIPT_24_04_02_DS_Ratio_1-1-Insilico_cbNIPT_24_04_05',
    #'cbNIPT_24_04_02_DS_Ratio_1-1-Insilico_cbNIPT_24_04_06',
    #'cbNIPT_24_04_02_DS_Ratio_1-2-Insilico_cbNIPT_24_04_04',
    #'cbNIPT_24_04_02_DS_Ratio_1-2-Insilico_cbNIPT_24_04_05',
    #'cbNIPT_24_04_02_DS_Ratio_1-2-Insilico_cbNIPT_24_04_06',
    #'cbNIPT_24_04_02_DS_Ratio_2-1-Insilico_cbNIPT_24_04_04',
    #'cbNIPT_24_04_02_DS_Ratio_2-1-Insilico_cbNIPT_24_04_05',
    #'cbNIPT_24_04_02_DS_Ratio_2-1-Insilico_cbNIPT_24_04_06',
    'cbNIPT_24_04_03_DS',
    'cbNIPT_24_04_04',
    'cbNIPT_24_04_05',
    'cbNIPT_24_04_06',
    'cbNIPT_24_04_07',
    'cbNIPT_24_04_08',
    'cbNIPT_24_04_09',
]
sample_summary_accumulator = []

# [CHECK] sample_set 변수가 기존 환경에 이미 정의되어 있다고 가정합니다.
for sample in sample_set:
    path_list = Path(path).glob(f'**/baf/{sample}.position_baf_detail.tsv')
    
    for trans_fragment_path in path_list:
        DataId = str(Path(trans_fragment_path).name).split('.')[0]
        summary_path = str(trans_fragment_path).replace(".position_baf_detail.tsv", ".summary.tsv")
        print(summary_path)
        if not os.path.exists(summary_path):
            summary_path = str(trans_fragment_path).replace(".position_baf_detail.tsv", ".position_baf_detail.summary.tsv")
            #continue
            
        # 1. 샘플 데이터 로드
        summary_df = pd.read_csv(summary_path, sep='\t')
        
        # 2. 다중 유전 지표 일괄 Log2 정규화 수행 (크기 바이어스 및 뎁스 스케일링 보정)
        try:
            norm_df = normalize_all_metrics_with_sex_log2(
                df=summary_df, 
                depth_col="raw_total_depth", 
                sex_threshold=0.0001
            )
            
            # 3. 샘플별 비교를 위해 염색체(Chromosome) 레벨로 1차 평균 집계(Aggregating)
            # 이렇게 해야 빈(Bin) 단위 데이터가 염색체 단위로 요약되어 샘플 간 비교가 수월해집니다.
            chrom_mean = norm_df.groupby('chrom')[['trans_log2_norm', 'hetero_log2_norm', 'homo_log2_norm', 'raw_total_depth']].mean().reset_index()
            
            # 비교 식별용 샘플 고유 ID 주입
            chrom_mean['Sample_ID'] = DataId
            sample_summary_accumulator.append(chrom_mean)
            
        except Exception as e:
            print(f"[-] Error processing {DataId}: {e}")
            continue

# -----------------------------------------------------------------
# 4. [SAMPLE COMPARISON MATRIX] 모든 샘플을 가로/세로로 엮는 피벗 매트릭스 생성
# -----------------------------------------------------------------
if sample_summary_accumulator:
    # 전체 샘플 데이터 통합
    master_df = pd.concat(sample_summary_accumulator, ignore_index=True)

    print("\n" + "="*70)
    print("[★] Cross-Sample Comparison Matrices & Plots Generation...")
    print("="*70)
    
    # 자연스러운 염색체 순서 정렬 (chr1~22, X, Y)
    def chrom_key(c):
        c_str = str(c).replace('chr', '')
        if c_str == 'X': return 23
        if c_str == 'Y': return 24
        try: return int(c_str)
        except: return 99        
    # [Matrix 생성]
    raw_chroms = master_df['chrom'].unique()
    sorted_chroms = sorted(raw_chroms, key=chrom_key)
    
    comparison_trans_matrix = master_df.pivot(index='chrom', columns='Sample_ID', values='trans_log2_norm').reindex(sorted_chroms)
    comparison_hetero_matrix = master_df.pivot(index='chrom', columns='Sample_ID', values='hetero_log2_norm').reindex(sorted_chroms)
    comparison_homo_matrix = master_df.pivot(index='chrom', columns='Sample_ID', values='homo_log2_norm').reindex(sorted_chroms)
    comparison_depth_matrix = master_df.pivot(index='chrom', columns='Sample_ID', values='raw_total_depth').reindex(sorted_chroms)

    # 콘솔 출력 (Trans 확인용)
    print("\n[Matrix 1] Log2 Normalized Trans Fragments Comparison (Preview)")
    print("-" * 70)
    print(comparison_trans_matrix.round(3).head())
    
    # -----------------------------------------------------------------
    # [NEW] 샘플 간 교차 비교 시각화 (Multi-panel Line Plot)
    # -----------------------------------------------------------------
    print("\n[★] Generating Cross-Sample Profile Plot...")
    
    fig, axes = plt.subplots(3, 1, figsize=(18, 14), sharex=True)
    
    # 그릴 매트릭스와 타이틀 매핑
    plot_targets = [
        (comparison_trans_matrix, "Log2 Normalized Trans Fragments (Phasing Intensity)", axes[0]),
        (comparison_hetero_matrix, "Log2 Normalized Hetero Sites Count (Marker Survival)", axes[1]),
        (comparison_homo_matrix, "Log2 Normalized Homo Sites Count (Background Baseline)", axes[2])
    ]
    
    # 색상 팔레트 자동 지정 (샘플 구분을 위해)
    colors = plt.cm.tab10(np.linspace(0, 1, len(comparison_trans_matrix.columns)))
    
    for matrix, title, ax in plot_targets:
        # 각 샘플(열)별로 라인 플롯 생성
        for idx, sample_id in enumerate(matrix.columns):
            ax.plot(matrix.index, matrix[sample_id], marker='o', markersize=5, 
                    linewidth=1.5, alpha=0.8, color=colors[idx], label=sample_id)
            
        ax.set_title(title, fontsize=14, fontweight='bold', loc='left')
        ax.set_ylabel("Log2 Ratio", fontsize=12)
        
        # ---------------------------------------------------------
        # [수정된 부분] 데이터 범위에 맞춰 Y축을 자동으로 늘려주는 로직
        # ---------------------------------------------------------
        # 현재 매트릭스에 있는 모든 데이터의 진짜 최소/최대값 추출
        data_min = matrix.min().min()
        data_max = matrix.max().max()
        
        # 임상 판독선(-1.0 ~ +0.6)은 무조건 화면에 나오도록 최소 바운더리 설정
        plot_min = min(data_min - 0.5, -1.5) # 최소 -1.5 보장, 그보다 낮으면 더 내림
        plot_max = max(data_max + 0.5, 1.5)  # 최대 +1.5 보장, 그보다 높으면 더 올림
        
        ax.set_ylim(plot_min, plot_max)
        # ---------------------------------------------------------
        
        # 임상 판독용 진단 가이드라인 투사
        ax.axhline(0.0, color='black', linestyle='-', linewidth=1.2, alpha=0.6) # 2n Baseline
        ax.axhline(0.3, color='orange', linestyle='--', linewidth=1, alpha=0.7) # Mosaicism Threshold
        ax.axhline(0.6, color='red', linestyle='--', linewidth=1, alpha=0.7) # Full 3n Threshold
        ax.axhline(-1.0, color='blue', linestyle='--', linewidth=1, alpha=0.7) # 1n Threshold
        
        ax.grid(True, axis='y', linestyle=':', alpha=0.5)
        
        # 범례는 첫 번째 패널 우측 바깥에만 표시
        if ax == axes[0]:
            ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1), fontsize=10, title="Samples", title_fontsize=11)
            
    # X축 꾸미기 (염색체 이름)
    axes[2].set_xticks(range(len(sorted_chroms)))
    axes[2].set_xticklabels([c.replace('chr', '') for c in sorted_chroms], fontsize=11, fontweight='bold')
    axes[2].set_xlabel("Chromosomes", fontsize=14)
    
    plt.suptitle("Cross-Sample Genetic Metric Comparison Profiles", fontsize=18, fontweight='bold', y=0.96)
    plt.tight_layout()
    #plt.subplots_adjust(right=2) # 범례가 짤리지 않게 우측 여백 확보
    
    # 결과 차트 저장
    output_plot_path = os.path.join(path, "Cross_Sample_Comparison_Profiles.png")
    plt.savefig(output_plot_path, dpi=300)
    plt.close()
    
    print("\n[★] Calculating Absolute Delta Scores for Aneuploidy Detection...")
    
    # 1. 종합 신호 강도(Composite Signal) 산출 (Trans 70% + Hetero 30% 가중치)
    master_df['composite_signal'] = (0.7 * master_df['trans_log2_norm']) + (0.3 * master_df['hetero_log2_norm'])
    
    diagnostic_scores = []
    
    # 2. 샘플별 독립 Delta 연산 (단일 샘플 내 기저선과의 절대 차이)
    for sample_id, group in master_df.groupby('Sample_ID'):
        auto_mask = ~group['chrom'].isin(['chrX', 'chrY'])
        auto_group = group[auto_mask]
        
        if auto_group.empty: continue
            
        sample_baseline = auto_group['composite_signal'].median()
        
        group_delta = group.copy()
        group_delta['delta_score'] = group_delta['composite_signal'] - sample_baseline
        
        # [핵심 변환] 1. 샘플 내부의 전반적인 노이즈 경향성(MAD) 측정
        # 전체 상염색체들이 기저선으로부터 평균적으로 얼마나 요동치는지 수치화합니다.
        noise_mad = group_delta[auto_mask]['delta_score'].abs().median()
        
        # [핵심 변환] 2. 유의미한 도약(Trend) 임계값 설정
        # "일반적인 노이즈 요동폭의 3배 이상"이면서 "최소 0.1 이상"의 도약만 진짜 경향성으로 인정합니다.
        trend_threshold = max(noise_mad * 3.0, 0.1)
        
        # 3. 최대 도약 염색체 추출
        peak_idx = group_delta[auto_mask]['delta_score'].idxmax()
        peak_row = group_delta.loc[peak_idx]
        raw_max_delta = peak_row['delta_score']
        
        # 4. 경향성 검증 (필터 기각)
        # 단일 세포 벌크 믹스 노이즈로 인해 우연히 튄 최댓값이라면 임계값에 걸려 0% 처리됩니다.
        if raw_max_delta >= trend_threshold:
            peak_chrom = peak_row['chrom']
            valid_delta = raw_max_delta
            emf_percent = max(0.0, (valid_delta / 0.6) * 100)
        else:
            peak_chrom = "None (Noise Cancelled)"
            valid_delta = 0.0
            emf_percent = 0.0
            
        display_percent = min(emf_percent, 110.0)
        
        diagnostic_scores.append({
            'Sample_ID': sample_id,
            'Peak_Chrom': peak_chrom,
            'Baseline': sample_baseline,
            'Noise_MAD': noise_mad,
            'Max_Delta': valid_delta,
            'EMF_Percent': emf_percent,
            'Display_Percent': display_percent
        })
        
    score_df = pd.DataFrame(diagnostic_scores).sort_values(by='Display_Percent', ascending=False)
    
    print("\n[Diagnostic Delta-Score Table]")
    print(score_df[['Sample_ID', 'Peak_Chrom', 'Max_Delta', 'EMF_Percent']].round(3))
    
    # 4. 순수 도약(Delta) 기반 임상 판독 Bar Chart 시각화
    plt.figure(figsize=(14, 8))
    
    bars = plt.bar(score_df['Sample_ID'], score_df['Display_Percent'], width=0.6, alpha=0.85, edgecolor='black')
    
    # 임상 진단 기준선 투사 (0.6을 100%로 잡았을 때의 백분율)
    plt.axhline(20, color='gray', linestyle='--', linewidth=1.5, label='Noise Cutoff (20% Mosaic)')
    plt.axhline(80, color='red', linestyle='--', linewidth=1.5, label='Full Trisomy Cutoff (80% Mosaic)')
    
    for bar, _, row in zip(bars, score_df.index, score_df.itertuples()):
        val = row.Display_Percent
        chrom = row.Peak_Chrom
        
        # 도약 퍼센트에 따른 직관적 컬러링
        if val >= 80.0:
            bar.set_color('#d62728')  # Red
            diag_text = f"Full 3n\n({chrom})"
        elif val >= 20.0:
            bar.set_color('#ff7f0e')  # Orange
            diag_text = f"Mosaic\n({chrom})"
        else:
            bar.set_color('#2ca02c')  # Green
            diag_text = f"Normal"
            
        plt.text(bar.get_x() + bar.get_width()/2, val + 2, 
                 f"{row.EMF_Percent:.1f}%\n(Δ{row.Max_Delta:.2f})\n{diag_text}", 
                 ha='center', va='bottom', fontsize=10, fontweight='bold')
                 
    plt.title("Sample-Level Aneuploidy Detection via Baseline Delta Shift", fontsize=16, fontweight='bold', pad=20)
    plt.ylabel("Estimated Mosaic Fraction (%)", fontsize=14)
    plt.ylim(0, 130)
    plt.xticks(rotation=30, ha='right', fontsize=11)
    
    # 구역 색칠
    plt.axhspan(20, 80, color='orange', alpha=0.1)
    plt.axhspan(80, 130, color='red', alpha=0.1)
    
    plt.legend(loc='upper right', fontsize=12)
    plt.grid(axis='y', linestyle=':', alpha=0.6)
    plt.tight_layout()
    
    # 차트 저장
    score_plot_path = os.path.join(path, "Delta_Diagnostic_Profiles.png")
    plt.savefig(score_plot_path, dpi=300)
    plt.close()
    
    score_df.to_csv(os.path.join(path, "Delta_Diagnostic_Scores.tsv"), sep='\t', index=False)
    print(f"[+] Output Diagnostic Delta Chart saved to: {score_plot_path}")
    print("="*70)


else:
    print("\n[-] Failed to generate sample comparison matrix. No data parsed.")