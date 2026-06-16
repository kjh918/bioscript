import os
import pysam
import numpy as np
import pandas as pd
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

from utils import log, ensure_dir, chrom_key
# [CRITICAL FIX] apply_qc_and_impute 임포트 제거
from processing import gc_correct_lowess, _parallel_chrom_worker
from normalization import normalize_by_chrom_with_sex, apply_qc_bin
from visualization import plot_comprehensive_baf_logr_qc
from segmentation import assign_cn_state, segment_one_cell


class CnvPipeline:
    def __init__(self, args):
        self.args = args
        self.dirs = {k: os.path.join(args.OutDir, k) for k in ["meta", "qc", "data", "plots", "segments", "baf"]}
        for d in self.dirs.values(): ensure_dir(d)

    def _load_or_run(self, target_file, run_func, *args, **kwargs):
        if os.path.exists(target_file):
            log(f"[Skip] Found existing file: {os.path.basename(target_file)}. Loading...")
            return pd.read_csv(target_file, sep="\t")
        log(f"[Run] Generating: {os.path.basename(target_file)}...")
        df = run_func(*args, **kwargs)
        df.to_csv(target_file, sep="\t", index=False)
        return df

    def load_target_bins(self):
        bin_file = self.args.AnnotatedBins
        if not os.path.exists(bin_file):
            raise FileNotFoundError(f"Cannot find bin file: {bin_file}")
        log(f"Loading pre-computed BED.GZ bins (with GC/Mappability & Tags): {bin_file}")
        df_bins = pd.read_csv(bin_file, sep="\t", compression='gzip', comment=None)
        if '#chrom' in df_bins.columns:
            df_bins.rename(columns={'#chrom': 'chrom'}, inplace=True)
        return df_bins

    # ---------------------------------------------------------
    # [모듈 1] One-Pass 병렬 추출: Coverage, BAF, 돌연변이 시그니처
    # ---------------------------------------------------------
    def extract_genetic_evidence(self, df_bins):
        target_file = os.path.join(self.dirs["data"], f"{self.args.SeqID}.bins_with_evidence.tsv")
        
        def _run():
            tmp_dir = os.path.join(self.dirs["data"], "tmp_workers")
            ensure_dir(tmp_dir)
            
            chrom_groups = [group.copy() for chrom, group in df_bins.groupby("chrom")]
            max_workers = getattr(self.args, 'Threads', 4)
            log(f"Parallel processing with {max_workers} threads (One-Pass BAM & VCF Extraction)...")
            
            valid_pos_files, valid_sum_files, valid_mut_files = [], [], []
            
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                thresholds = {
                    'homo_max': getattr(self.args, 'HomoMax', 0.1),
                    'hetero_min': getattr(self.args, 'HeteroMin', 0.4),
                    'hetero_max': getattr(self.args, 'HeteroMax', 0.6),
                    'homo_min': getattr(self.args, 'HomoMin', 0.9)
                }
                futures = {
                    executor.submit(_parallel_chrom_worker, group, self.args.BamPath, self.args.VcfFile, self.args.MinMapQ, thresholds, tmp_dir): group['chrom'].iloc[0]
                    for group in chrom_groups
                }
                for future in as_completed(futures):
                    pos_path, sum_path, mut_path = future.result()
                    if pos_path: valid_pos_files.append(pos_path)
                    if sum_path: valid_sum_files.append(sum_path)
                    if mut_path: valid_mut_files.append(mut_path)

            log("Merging worker results...")
            
            if valid_pos_files:
                full_pos_report = pd.concat([pd.read_parquet(f) for f in valid_pos_files], ignore_index=True)    
                baf_out_dir = os.path.join(self.dirs["baf"], f"{self.args.SeqID}_position_baf")
                log(f"Saving partitioned position details to {baf_out_dir}")
                full_pos_report.to_parquet(baf_out_dir, engine='pyarrow', partition_cols=['chrom'], index=False)
                
            if valid_mut_files:
                full_mut_report = pd.concat([pd.read_parquet(f) for f in valid_mut_files], ignore_index=True)
                mut_out_path = os.path.join(self.dirs["data"], f"{self.args.SeqID}.mutation_signatures.tsv")
                full_mut_report.to_csv(mut_out_path, sep="\t", index=False)
                log(f"Saved mutation signature distributions.")

            evidence_df = pd.concat([pd.read_parquet(f) for f in valid_sum_files], ignore_index=True)
            
            mean_density_informative = evidence_df['total_sites'].mean()
            evidence_df['Informative_OE_Ratio'] = evidence_df['total_sites'] / (mean_density_informative + 1e-9)
            evidence_df['Informative_Norm_Log2'] = np.log2(evidence_df['Informative_OE_Ratio'] + 1e-9)
            
            for f in valid_pos_files + valid_sum_files + valid_mut_files: 
                if os.path.exists(f): os.remove(f)
            if os.path.exists(tmp_dir): os.rmdir(tmp_dir)
            
            df_merged = df_bins.merge(evidence_df, on=["chrom", "start", "end"], how="left")
            
            df_merged['raw_count'] = df_merged['raw_count'].fillna(0)
            df_merged['breadth_ratio'] = df_merged['breadth_ratio'].fillna(0.0)
            
            return df_merged
            
        return self._load_or_run(target_file, _run)

    # ---------------------------------------------------------
    # [모듈 2] GC 보정 및 Wavelet Denoising (Imputation 제거됨)
    # ---------------------------------------------------------
    def correct_gc_and_normalize(self, df_merged):
        target_file = os.path.join(self.dirs["data"], f"{self.args.SeqID}.normalized.tsv")
        def _run():
            log("Applying Standard GC Correction and Normalization...")
            
            # 1. LOWESS 기반 GC 편향 보정 (있는 그대로의 데이터를 사용)
            df_gc, gc_stats = gc_correct_lowess(df_merged, frac=self.args.LowessFrac) 
            # 2. 보정값 정규화 (염색체 및 성별 기반 Baseline 세팅)
            norm_result = normalize_by_chrom_with_sex(df_gc, value_col="log2_chrom_norm")
            
            # Tuple 안전 추출
            if isinstance(norm_result, tuple):
                df_norm = norm_result[0]
            else:
                df_norm = norm_result
            
            # 복잡한 디노이징 없이 정직하게 정규화된 값만 반환
            return df_norm
        _run()
        return self._load_or_run(target_file, _run)

    # ---------------------------------------------------------
    # [모듈 3] 다차원 세그멘테이션 (Multi-Dimensional PELT)
    # ---------------------------------------------------------
    def perform_segmentation(self, df_global):
        target_file = os.path.join(self.dirs["data"], f"{self.args.SeqID}.segments.tsv")
        def _run():
            log("Running 1D PELT Segmentation on Copy Number (log2_chrom_norm)...")
            
            # 옛날 방식의 import로 원복 (segmentation.py 내부 함수명)
            from segmentation import segment_one_cell, assign_cn_state
            
            segments = segment_one_cell(
                meta=df_global, 
                df_signals=df_global, 
                signal_col="log2_chrom_norm", # 오직 Copy Number 하나만 봄!
                penalty=self.args.SegPenalty
            )
            return assign_cn_state(segments, baseline_ploidy=self.args.BaselinePloidy)
            
        _run()
        return self._load_or_run(target_file, _run)

    # ---------------------------------------------------------
    # 파이프라인 메인 실행 루프
    # ---------------------------------------------------------
    def run(self):
        log(f"--- Execution Start: {self.args.SeqID} ---")
        
        df_bins = self.load_target_bins()
        
        df_merged_evidence = self.extract_genetic_evidence(df_bins)
        
        df_filterd = apply_qc_bin(df_merged_evidence)
        print(df_filterd)

        df_global = self.correct_gc_and_normalize(df_filterd)
        
        segment_df = self.perform_segmentation(df_global)
        print(segment_df)

        plot_comprehensive_baf_logr_qc(
            df_bins, segment_df, os.path.join(self.dirs["plots"], f"{self.args.SeqID}.comprehensive_view.png"), self.args.SeqID
        )


