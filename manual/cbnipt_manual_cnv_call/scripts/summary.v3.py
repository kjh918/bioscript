

import sys, os, argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
from pathlib import Path
from glob import glob

sys.path.append(str(Path(__file__).parent.parent))  # 프로젝트 루트 디렉토리 추가

from src.cbnipt_cnv_caller.visualization import plot_chromosome_overview, plot_bin_distribution, plot_score_heatmap, plot_final_call
from src.cbnipt_cnv_caller.utils import sort_chroms
from src.cbnipt_cnv_caller.normalization import normalize_and_estimate_sex
from src.cbnipt_cnv_caller.classification import apply_qc_filter, compute_chrom_summary, analyze_all_chromosomes



def process_and_run(raw_tsv_path: str):
    tsv_path = Path(raw_tsv_path)
    sample_id = tsv_path.name.split(".")[0]
    out_dir = tsv_path.parent / "Standardized_Results"
    out_dir.mkdir(parents=True, exist_ok=True)
    
    print("\n" + "=" * 75)
    print(f"[★] Unified Master Pipeline v3.0: {sample_id}")
    print("=" * 75)
    
    raw_df = pd.read_csv(tsv_path, sep="\t")
    
    # 1. 정규화 및 성별 추출
    print(" [STEP 1] Applying Normalization (LOO & Sex correction)...")
    norm_df, sex_tag, x_ratio, y_ratio = normalize_and_estimate_sex(raw_df, value_col="raw_count")
    print(sex_tag, x_ratio, y_ratio )
    #exit()
    # 2. 정규화 데이터 저장
    standardized_path = out_dir / f"{sample_id}.standardized.tsv"
    norm_df.to_csv(standardized_path, sep="\t", index=False)
    print(f" [STEP 2] Standardized TSV saved: {standardized_path.name}")
    
    # 3. 진단 스코어링
    print(" [STEP 3] Running Anomaly Detection Pipeline...")
    #df_f = apply_qc_filter(norm_df)
    df_f = norm_df

    summary = compute_chrom_summary(df_f)
    call_df = analyze_all_chromosomes(summary, sex_tag, x_ratio, y_ratio)
    print(call_df)
    
    
    print(f"\n{'Chrom':<8} {'Norm_Log2':>9} {'RobustZ':>8} {'Final':>7} Call")
    print("-" * 70)
    for _, row in call_df.set_index("chrom").reindex(sort_chroms(call_df["chrom"].tolist())).reset_index().iterrows():
        flag  = "⚠ " if row["call"] == "ABNORMAL" else ("△ " if row["call"] == "SUSPICIOUS" else "  ")
        log2  = f"{row['log2fc']:+.3f}" if not np.isnan(row["log2fc"]) else "   N/A"
        rz    = f"{row['robust_z']:+.2f}" if not np.isnan(row["robust_z"]) else "  N/A"
        print(f"{row['chrom']:<8} {log2:>9} {rz:>8} {row['final_score']:>7.3f}  {flag}{row['call']:<12} {row.get('detail','')}")

    # 4. 4종 전체 시각화
    print(" [STEP 4] Generating All Figures (Overview, Violin, Heatmap, Bar)...")
    plot_chromosome_overview(df_f, summary, call_df, sex_tag, sample_id, str(out_dir))
    plot_bin_distribution(df_f, call_df, sample_id, str(out_dir))
    plot_score_heatmap(call_df, sample_id, sex_tag, str(out_dir))
    plot_final_call(call_df, sex_tag, sample_id, str(out_dir))

    print(f"\n[Done] All process successfully completed. Output: {out_dir}")

if __name__ == "__main__":
    file_list = glob("/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Results/cbNIPT_24_04_02_DS_Ratio_1-1-Insilico_cbNIPT_24_04_04/data/cbNIPT_24_04_02_DS_Ratio_1-1-Insilico_cbNIPT_24_04_04.normalized.tsv")
    if not file_list: print("No files found to process.")
    else:
        for TSV_PATH in file_list[:1]: 
            print(TSV_PATH)
            process_and_run(TSV_PATH)