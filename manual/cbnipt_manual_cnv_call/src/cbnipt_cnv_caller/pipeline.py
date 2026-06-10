import os
import pysam
import numpy as np
import pandas as pd
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

from utils import log, ensure_dir
from processing import apply_low_quality_filter, gc_correct_lowess, _parallel_chrom_worker
from normalization import normalize_by_chrom_with_sex, normalize_all_metrics_with_sex_log2
from visualization import plot_gc_correction, plot_comprehensive_baf_logr_qc
from segmentation import segment_one_cell, assign_cn_state




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
        log(f"Loading pre-computed BED.GZ bins (with GC content): {bin_file}")
        df_bins = pd.read_csv(bin_file, sep="\t", compression='gzip', comment=None)
        if '#chrom' in df_bins.columns:
            df_bins.rename(columns={'#chrom': 'chrom'}, inplace=True)
        return df_bins

    # ---------------------------------------------------------
    # [모듈 1] One-Pass 병렬 추출: Coverage와 BAF 동시 추출
    # ---------------------------------------------------------
    def extract_genetic_evidence(self, df_bins):
        target_file = os.path.join(self.dirs["data"], f"{self.args.SeqID}.bins_with_evidence.tsv")
        
        def _run():
            tmp_dir = os.path.join(self.dirs["data"], "tmp_workers")
            ensure_dir(tmp_dir)
            
            chrom_groups = [group.copy() for chrom, group in df_bins.groupby("chrom")]
            max_workers = getattr(self.args, 'Threads', 4)
            log(f"Parallel processing with {max_workers} threads (One-Pass BAM & VCF Extraction)...")
            
            valid_pos_files, valid_sum_files = [], []
            
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
                    pos_path, sum_path = future.result()
                    if pos_path: valid_pos_files.append(pos_path)
                    if sum_path: valid_sum_files.append(sum_path)

            log("Merging worker results...")
            if valid_pos_files:
                full_pos_report = pd.concat([pd.read_parquet(f) for f in valid_pos_files], ignore_index=True)    
                # 저장될 폴더 경로 (예: baf/cbNIPT_24_04_01_DS_position_baf)
                baf_out_dir = os.path.join(self.dirs["baf"], f"{self.args.SeqID}_position_baf")
                log(f"Saving partitioned position details to {baf_out_dir}")
                
                # chrom 컬럼을 기준으로 자동으로 하위 폴더를 만들어 초고속 저장
                full_pos_report.to_parquet(
                    baf_out_dir, 
                    engine='pyarrow', 
                    partition_cols=['chrom'], 
                    index=False
                )
            evidence_df = pd.concat([pd.read_parquet(f) for f in valid_sum_files], ignore_index=True)
            
            mean_density_informative = evidence_df['total_sites'].mean()
            
            evidence_df['Informative_OE_Ratio'] = evidence_df['total_sites'] / (mean_density_informative + 1e-9)
            evidence_df['Informative_Norm_Log2'] = np.log2(evidence_df['Informative_OE_Ratio'] + 1e-9)
            
            # 임시 워커 파일 정리
            for f in valid_pos_files + valid_sum_files: os.remove(f)
            os.rmdir(tmp_dir)
            
            df_merged = df_bins.merge(evidence_df, on=["chrom", "start", "end"], how="left")
            
            # 결측치 처리: 새로 도입된 raw_count와 breadth_ratio 반영
            df_merged['raw_count'] = df_merged['raw_count'].fillna(0)
            df_merged['breadth_ratio'] = df_merged['breadth_ratio'].fillna(0.0)
            
            return df_merged
            
        return self._load_or_run(target_file, _run)

    # ---------------------------------------------------------
    # [모듈 2] Coverage 저품질 필터링 및 GC 보정 (LOWESS)
    # ---------------------------------------------------------
    def correct_gc_and_normalize(self, df_merged):
        target_file = os.path.join(self.dirs["data"], f"{self.args.SeqID}.normalized.tsv")
        
        def _run():
            # 1. 성염색체 보호 및 상염색체 고도화 필터링 (새로운 apply_low_quality_filter 사용)
            min_cov = getattr(self.args, 'MinCoverage', 0.5)
            df_clean = apply_low_quality_filter(df_merged, min_depth=self.args.MinDepth, min_coverage=min_cov)
            
            # 2. LOWESS 기반 GC 보정 (내부적으로 'raw_count'와 'gc' 컬럼 활용)
            df_gc, gc_stats = gc_correct_lowess(df_clean, frac=self.args.LowessFrac) 
            
            # 3. 보정값 정규화 (염색체 및 성별 기반)
            return normalize_by_chrom_with_sex(df_gc, value_col="log2_corrected")
            
        return self._load_or_run(target_file, _run)

    # ---------------------------------------------------------
    # [모듈 3] 세그멘테이션 및 통계
    # ---------------------------------------------------------
    def perform_segmentation(self, df_global):
        target_file = os.path.join(self.dirs["data"], f"{self.args.SeqID}.segments.tsv")
        def _run():
            segments = segment_one_cell(df_global, df_global["log2_chrom_norm"], self.args.SegPenalty)
            return assign_cn_state(segments, baseline_ploidy=self.args.BaselinePloidy)
        return self._load_or_run(target_file, _run)

    # ---------------------------------------------------------
    # 파이프라인 실행 로직 (Main)
    # ---------------------------------------------------------
    def run(self):
        log(f"--- Execution Start: {self.args.SeqID} ---")
        
        # 1. 뼈대 로드 (GC% 포함)
        df_bins = self.load_target_bins()
        
        # 2. 데이터 추출 (BAM 단일 패스 -> Coverage + BAF 동시 추출)
        df_merged_evidence = self.extract_genetic_evidence(df_bins)
        
        # 3. GC 편향 보정 (LOWESS) 및 정규화
        df_global = self.correct_gc_and_normalize(df_merged_evidence)
        
        # 4. 세그멘테이션 (CNV Call)
        segments = self.perform_segmentation(df_global)
        
        # 5. 최종 시각화
        plot_comprehensive_baf_logr_qc(
            df_global, segments, 
            os.path.join(self.dirs["plots"], f"{self.args.SeqID}.comprehensive_view.png"), 
            self.args.SeqID
        )
        
        log(f"--- [COMPLETE] Results saved to {self.args.OutDir} ---")