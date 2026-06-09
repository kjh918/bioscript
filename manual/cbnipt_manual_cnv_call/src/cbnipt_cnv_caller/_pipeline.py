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
from processing import process_bam_to_coverage, apply_low_quality_filter, gc_correct_lowess
from normalization import normalize_by_chrom_with_sex, normalize_all_metrics_with_sex_log2
from visualization import plot_gc_correction, plot_genome_wide_cnv_with_gaps
from segmentation import segment_one_cell, assign_cn_state

# [중요] 최적화된 유전적 증거 추출 워커 임포트
from utils import _parallel_chrom_worker

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
        log(f"Loading pre-computed BED.GZ bins: {bin_file}")
        df_bins = pd.read_csv(bin_file, sep="\t", compression='gzip', comment=None)
        if '#chrom' in df_bins.columns:
            df_bins.rename(columns={'#chrom': 'chrom'}, inplace=True)
        return df_bins

    def process_coverage_and_filter(self, filtered_bins):
        target_file = os.path.join(self.dirs["data"], f"{self.args.SeqID}.bins_with_reads_filtered.tsv")
        def _run():
            df_with_reads, filter_stats = process_bam_to_coverage(
                self.args.BamPath, filtered_bins, min_mapq=self.args.MinMapQ
            )
            df_clean, filter_stats_df = apply_low_quality_filter(
                df_with_reads, filter_stats, min_depth=self.args.MinDepth, min_coverage=self.args.MinCoverage
            )
            return df_clean
        return self._load_or_run(target_file, _run)

    # ---------------------------------------------------------
    # [모듈 3] 고성능 병렬 유전 증거 추출 (Parquet Disk-Spill)
    # ---------------------------------------------------------
    def extract_genetic_evidence(self, df_clean):
        target_file = os.path.join(self.dirs["data"], f"{self.args.SeqID}.bins_with_evidence.tsv")
        
        def _run():
            tmp_dir = os.path.join(self.dirs["data"], "tmp_workers")
            ensure_dir(tmp_dir)
            
            chrom_groups = [group.copy() for chrom, group in df_clean.groupby("chrom")]
            max_workers = getattr(self.args, 'Threads', 4)
            log(f"Parallel processing with {max_workers} threads...")
            
            valid_pos_files, valid_sum_files = [], []
            
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                futures = {
                    executor.submit(_parallel_chrom_worker, group, self.args.BamPath, self.args.VcfFile, self.args.MinMapQ, tmp_dir): group['chrom'].iloc[0]
                    for group in chrom_groups
                }
                for future in as_completed(futures):
                    pos_path, sum_path = future.result()
                    if pos_path: valid_pos_files.append(pos_path)
                    if sum_path: valid_sum_files.append(sum_path)

            # 병합
            log("Merging worker results...")
            full_pos_report = pd.concat([pd.read_parquet(f) for f in valid_pos_files], ignore_index=True)
            full_pos_report.to_csv(os.path.join(self.dirs["baf"], f"{self.args.SeqID}.position_baf_detail.tsv"), sep="\t", index=False)
            
            evidence_df = pd.concat([pd.read_parquet(f) for f in valid_sum_files], ignore_index=True)
            
            # 앵커 기반 정규화 및 Delta 연산
            ANCHOR_CHROMS = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5'] 
            anchor_df = evidence_df[evidence_df['chrom'].isin(ANCHOR_CHROMS)]
            anchor_hetero, anchor_homo = max(anchor_df['hetero_sites_count'].sum(), 1), max(anchor_df['homo_sites_count'].sum(), 1)
            
            evidence_df['Ratio_Hetero'] = evidence_df['hetero_sites_count'] / anchor_hetero
            evidence_df['Ratio_Homo'] = evidence_df['homo_sites_count'] / anchor_homo
            
            expected = evidence_df[~evidence_df['chrom'].isin(['chrX', 'chrY'])].groupby('chrom')[['Ratio_Hetero', 'Ratio_Homo']].median()
            evidence_df = evidence_df.merge(expected, on='chrom', suffixes=('', '_expected'))
            
            evidence_df['Hetero_Delta_Log2'] = np.log2((evidence_df['Ratio_Hetero'] + 1e-9) / (evidence_df['Ratio_Hetero_expected'] + 1e-9))
            evidence_df['Homo_Delta_Log2'] = np.log2((evidence_df['Ratio_Homo'] + 1e-9) / (evidence_df['Ratio_Homo_expected'] + 1e-9))
            
            # 임시 파일 정리
            for f in valid_pos_files + valid_sum_files: os.remove(f)
            os.rmdir(tmp_dir)
            
            return df_clean.merge(evidence_df, on=["chrom", "start", "end"], how="left")
            
        return self._load_or_run(target_file, _run)

    def correct_gc_and_normalize(self, df_clean):
        target_file = os.path.join(self.dirs["data"], f"{self.args.SeqID}.normalized.tsv")
        def _run():
            df_gc, gc_stats = gc_correct_lowess(df_clean, frac=self.args.LowessFrac)
            return normalize_by_chrom_with_sex(df_gc, value_col="log2_corrected")
        return self._load_or_run(target_file, _run)

    def perform_segmentation(self, df_global):
        target_file = os.path.join(self.dirs["data"], f"{self.args.SeqID}.segments.tsv")
        def _run():
            segments = segment_one_cell(df_global, df_global["log2_chrom_norm"], self.args.SegPenalty)
            return assign_cn_state(segments, baseline_ploidy=self.args.BaselinePloidy)
        return self._load_or_run(target_file, _run)

    def run(self):
        log(f"--- Execution Start: {self.args.SeqID} ---")
        filtered_bins = self.load_target_bins()
        df_clean = self.process_coverage_and_filter(filtered_bins)
        df_clean_with_evidence = self.extract_genetic_evidence(df_clean)
        df_global = self.correct_gc_and_normalize(df_clean_with_evidence)
        segments = self.perform_segmentation(df_global)
        
        plot_genome_wide_cnv_with_gaps(
            df_global, segments, 
            os.path.join(self.dirs["plots"], f"{self.args.SeqID}.genome_wide.png"), self.args.SeqID
        )
        log(f"--- [COMPLETE] Results saved to {self.args.OutDir} ---")