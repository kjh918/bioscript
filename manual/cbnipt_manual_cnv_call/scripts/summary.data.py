import os
from pathlib import Path
import numpy as np
import pandas as pd

# 커스텀 정규화 모듈 임포트
from normalization import normalize_all_metrics_with_sex_log2

path = '/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Results/manual'

# 모든 샘플의 정규화 결과를 하나로 모으기 위한 리스트

sample_set = [
    'cbNIPT_24_04_01_DS',
    'cbNIPT_24_04_02_DS',
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
        summary_path = str(trans_fragment_path).replace(".position_baf_detail.tsv", ".position_baf_detail.summary.tsv")
        
        if not os.path.exists(summary_path):
            continue
            
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
    print("[★] Cross-Sample Comparison Matrices Generation Complete!")
    print("="*70)
    
    # [Matrix A] 샘플별 X 염색체별 Trans Fragments 비교 매트릭스 (행: 염색체, 열: 샘플)
    comparison_trans_matrix = master_df.pivot(index='chrom', columns='Sample_ID', values='trans_log2_norm')
    
    # 자연스러운 염색체 순서 정렬 (chr1~22, X, Y)
    def chrom_key(c):
        c_str = str(c).replace('chr', '')
        if c_str == 'X': return 23
        if c_str == 'Y': return 24
        try: return int(c_str)
        except: return 99
    
    sorted_chroms = sorted(comparison_trans_matrix.index, key=chrom_key)
    comparison_trans_matrix = comparison_trans_matrix.reindex(sorted_chroms)
    
    # [Matrix B] 샘플별 X 염색체별 Hetero Sites (ADO 우회 지표) 비교 매트릭스
    comparison_hetero_matrix = master_df.pivot(index='chrom', columns='Sample_ID', values='hetero_log2_norm').reindex(sorted_chroms)

    # [Matrix C] 샘플별 시퀀싱 총 뎁스(Yield) 비교 매트릭스
    comparison_depth_matrix = master_df.pivot(index='chrom', columns='Sample_ID', values='raw_total_depth').reindex(sorted_chroms)

    # 콘솔 출력 확인용 (Trans Fragment 비교 창)
    print("\n[Matrix 1] Log2 Normalized Trans Fragments Comparison (Row: Chrom, Col: Sample)")
    print("-" * 70)
    print(comparison_trans_matrix.round(4))
    
    print("\n[Matrix 2] Log2 Normalized Hetero Sites Count Comparison (Row: Chrom, Col: Sample)")
    print("-" * 70)
    print(comparison_hetero_matrix.round(4))

    # 후속 분석이나 액셀 검토를 위해 데이터프레임을 csv 파일로 바로 추출 보관할 수 있습니다.
    # comparison_trans_matrix.to_csv(os.path.join(path, "cross_sample_trans_matrix.tsv"), sep="\t")
    # comparison_hetero_matrix.to_csv(os.path.join(path, "cross_sample_hetero_matrix.tsv"), sep="\t")
else:
    print("\n[-] Failed to generate sample comparison matrix. No data parsed.")