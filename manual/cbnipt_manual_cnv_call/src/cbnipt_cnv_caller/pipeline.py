import os
import sys
import pysam
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim

from concurrent.futures import ProcessPoolExecutor, as_completed

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

from utils import log, ensure_dir, chrom_key, sort_chroms
from processing import gc_correct_lowess_robust, _parallel_chrom_worker
from normalization import normalize_by_chrom_with_sex, apply_qc_bin
from visualization import (
    plot_comprehensive_baf_logr_qc,
    plot_chromosome_overview,
    plot_score_heatmap,
    plot_final_call,
    plot_gc_correction,
    plot_cohort_overview,
    load_cohort_call_dfs,
)
from segmentation import assign_cn_state, segment_one_cell, rolling_micro_cnv_segmentation
from classification import compute_chrom_summary, analyze_all_chromosomes, diagnose_clinical_markers

# ═════════════════════════════════════════════════════════════
# CNV Pipeline
# ═════════════════════════════════════════════════════════════
class CnvPipeline:
    def __init__(self, args):
        self.args = args

        self.dirs = {
            k: os.path.join(args.OutDir, k)
            for k in ["meta", "qc", "data", "plots", "segments", "baf"]
        }

        for d in self.dirs.values():
            ensure_dir(d)

        self._gender_tag = "Unknown"
        self._x_ratio = 0.0
        self._y_ratio = 0.0
        self._norm_qc = {}

    # ─────────────────────────────────────────────────────────────
    # [MAIN] 실행 루프 (Execution Flow)
    # ─────────────────────────────────────────────────────────────
    def run(self):
        log("=" * 60)
        log(f"  Full Pipeline: {self.args.SeqID}")
        log("=" * 60)

        # 1. 대상 Bin 파일 로드
        df_bins = self.load_target_bins()

        # 2. BAM/VCF 파일로부터 증거(카운트, BAF 등) 병렬 추출
        df_merged_evidence = self.extract_genetic_evidence(df_bins)

        # 3. 불량 Bin QC 필터링 (GC 보정 전에 수행하는 것이 정석)
        df_filtered = apply_qc_bin(df_merged_evidence, min_depth=1, min_coverage=0.3)

        # 4. GC 보정(LOWESS) 및 LOO 정규화 (기준점: 2.0)
        df_global = self.correct_gc_and_normalize(df_filtered)

        # 5. 성별 판별 및 X/Y Ratio 계산
        self._determine_fetal_sex(df_global)

        # 6. 염색체 수준 분류 (Classification)
        summary = self.classify_chromosomes(df_global)
        call_df = analyze_all_chromosomes(summary, self._gender_tag, self._x_ratio, self._y_ratio)

        # call_df 저장: 코호트(다중 샘플) 시각화가 나중에 이 파일을 읽어서 종합 뷰를 만든다.
        call_df_path = os.path.join(self.dirs["data"], f"{self.args.SeqID}.chrom_calls.tsv")
        call_df.to_csv(call_df_path, sep="\t", index=False)

        # 7. 전체 염색체 단위 분절화 (Macro CNV Segmentation)
        segment_df = self.perform_segmentation(df_global)

        # 8. 미세결실 탐지 (Micro CNV Segmentation)
        microdeletion_df = self.perform_microdeletion(df_global)

        # 9. 탐지된 결과와 Bin 데이터를 모아서 최종 임상 진단 리포트 생성
        clinical_report = self.generate_clinical_report(df_global, segment_df, microdeletion_df)
        print(clinical_report)

        # 10. 시각화 (Plots)
        self.generate_plots(df_global, summary, call_df, segment_df)

        log("=" * 60)
        log(f"  Pipeline Complete: {self.args.SeqID}")
        log("=" * 60)

        return {
            "df_global": df_global,
            "summary": summary,
            "call_df": call_df,
            "segment_df": segment_df,
            "microdeletion_df": microdeletion_df,
            "clinical_report": clinical_report,
        }

    # ─────────────────────────────────────────────────────────────
    # [UTIL] 캐싱 및 실행 제어기
    # ─────────────────────────────────────────────────────────────
    def _load_or_run(self, target_file, run_func, *args, **kwargs):
        if os.path.exists(target_file):
            log(f"[Skip] Found existing file: {os.path.basename(target_file)}. Loading...")
            return pd.read_csv(target_file, sep="\t")

        log(f"[Run] Generating: {os.path.basename(target_file)}...")
        df = run_func(*args, **kwargs)

        if df is not None and not df.empty:
            df.to_csv(target_file, sep="\t", index=False)

        return df

    # ─────────────────────────────────────────────────────────────
    # [모듈 1] Bin 로드
    # ─────────────────────────────────────────────────────────────
    def load_target_bins(self):
        bin_file = self.args.AnnotatedBins

        if not os.path.exists(bin_file):
            raise FileNotFoundError(f"Cannot find bin file: {bin_file}")

        log(f"Loading pre-computed BED.GZ bins: {bin_file}")

        df_bins = pd.read_csv(
            bin_file,
            sep="\t",
            compression="gzip",
            comment=None,
        )

        if "#chrom" in df_bins.columns:
            df_bins.rename(columns={"#chrom": "chrom"}, inplace=True)

        return df_bins

    # ─────────────────────────────────────────────────────────────
    # [모듈 2] One-Pass 병렬 추출 BAM + VCF
    # ─────────────────────────────────────────────────────────────
    def extract_genetic_evidence(self, df_bins):
        target_file = os.path.join(
            self.dirs["data"],
            f"{self.args.SeqID}.bins_with_evidence.tsv",
        )

        def _run():
            tmp_dir = os.path.join(self.dirs["data"], "tmp_workers")
            ensure_dir(tmp_dir)

            chrom_groups = [
                group.copy()
                for _, group in df_bins.groupby("chrom")
            ]

            max_workers = getattr(self.args, "Threads", 4)
            log(f"Parallel processing with {max_workers} threads...")

            valid_pos_files = []
            valid_sum_files = []
            valid_mut_files = []

            thresholds = {
                "homo_max": getattr(self.args, "HomoMax", 0.1),
                "hetero_min": getattr(self.args, "HeteroMin", 0.4),
                "hetero_max": getattr(self.args, "HeteroMax", 0.6),
                "homo_min": getattr(self.args, "HomoMin", 0.9),
            }

            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                futures = {
                    executor.submit(
                        _parallel_chrom_worker,
                        group,
                        self.args.BamPath,
                        self.args.VcfFile,
                        self.args.MinMapQ,
                        thresholds,
                        tmp_dir,
                    ): group["chrom"].iloc[0]
                    for group in chrom_groups
                }

                for future in as_completed(futures):
                    pos_path, sum_path, mut_path = future.result()

                    if pos_path:
                        valid_pos_files.append(pos_path)
                    if sum_path:
                        valid_sum_files.append(sum_path)
                    if mut_path:
                        valid_mut_files.append(mut_path)

            log("Merging worker results...")

            if valid_pos_files:
                full_pos = pd.concat(
                    [pd.read_parquet(f) for f in valid_pos_files],
                    ignore_index=True,
                )
                baf_out = os.path.join(
                    self.dirs["baf"],
                    f"{self.args.SeqID}_position_baf",
                )
                full_pos.to_parquet(
                    baf_out,
                    engine="pyarrow",
                    partition_cols=["chrom"],
                    index=False,
                )

            if valid_mut_files:
                full_mut = pd.concat(
                    [pd.read_parquet(f) for f in valid_mut_files],
                    ignore_index=True,
                )
                full_mut.to_csv(
                    os.path.join(
                        self.dirs["data"],
                        f"{self.args.SeqID}.mutation_signatures.tsv",
                    ),
                    sep="\t",
                    index=False,
                )

            if not valid_sum_files:
                raise RuntimeError("No valid summary files generated from workers.")

            evidence_df = pd.concat(
                [pd.read_parquet(f) for f in valid_sum_files],
                ignore_index=True,
            )

            for f in valid_pos_files + valid_sum_files + valid_mut_files:
                if os.path.exists(f):
                    os.remove(f)
            try:
                os.rmdir(tmp_dir)
            except OSError:
                pass

            df_merged = df_bins.merge(
                evidence_df,
                on=["chrom", "start", "end"],
                how="left",
            )

            df_merged["raw_count"] = df_merged["raw_count"].fillna(0)
            df_merged["breadth_ratio"] = df_merged["breadth_ratio"].fillna(0.0)

            return df_merged

        return self._load_or_run(target_file, _run)

    # ─────────────────────────────────────────────────────────────
    # [모듈 3] GC 보정 + LOO 정규화
    # ─────────────────────────────────────────────────────────────
    def correct_gc_and_normalize(self, df_merged):
        target_file = os.path.join(
            self.dirs["data"],
            f"{self.args.SeqID}.normalized.new.tsv",
        )

        def _run():
            log("Applying GC Correction (LOWESS)...")
            df_gc, gc_stats = gc_correct_lowess_robust(
                df_merged,
                frac=self.args.LowessFrac,
            )

            if gc_stats is not None:
                plot_gc_correction(
                    gc_stats,
                    os.path.join(
                        self.dirs["qc"],
                        f"{self.args.SeqID}.gc_correction.png",
                    ),
                )
            else:
                log("GC correction plot skipped: gc_stats is None.")

            log("Normalizing (LOO, sex-aware)...")

            res = normalize_by_chrom_with_sex(df_gc, value_col="log2_corrected")
            df_norm = res[0] if isinstance(res, tuple) else res

            return df_norm

        return self._load_or_run(target_file, _run)

    # ─────────────────────────────────────────────────────────────
    # [모듈 4] 성별 판별 로직 및 Ratio 갱신
    # ─────────────────────────────────────────────────────────────
    def _determine_fetal_sex(self, df_global):
        log("Determining Sex based on normalized copy number signals...")

        auto = df_global[~df_global["chrom"].isin(["chrX", "chrY"])]["copy_number_signal"]
        auto_med = auto.median()

        if auto_med == 0 or np.isnan(auto_med):
            auto_med = 2.0

        x_mask = df_global["chrom"] == "chrX"
        y_mask = df_global["chrom"] == "chrY"

        x_median = df_global.loc[x_mask, "copy_number_signal"].median() if x_mask.sum() > 0 else 0.0
        y_median = df_global.loc[y_mask, "copy_number_signal"].median() if y_mask.sum() > 0 else 0.0

        self._x_ratio = float(x_median / auto_med)
        self._y_ratio = float(y_median / auto_med)

        if x_median > 1.5 and y_median < 0.2:
            self._gender_tag = 'XX'
        elif x_median < 1.5 and y_median > 0.2:
            self._gender_tag = "XY"
        else:
            self._gender_tag = "Unknown"

        log(f"Determined Sex: {self._gender_tag} (X Copy: {x_median:.2f}, Y Copy: {y_median:.2f}, "
            f"X Ratio: {self._x_ratio:.2f}, Y Ratio: {self._y_ratio:.2f})")

    # ─────────────────────────────────────────────────────────────
    # [모듈 5] Classification
    # ─────────────────────────────────────────────────────────────
    def classify_chromosomes(self, df_global):
        log("Computing chromosome-level summary statistics...")
        summary = compute_chrom_summary(df_global)
        return summary

    # ─────────────────────────────────────────────────────────────
    # [모듈 6] 전체 염색체 단위 분절화 (Macro)
    # ─────────────────────────────────────────────────────────────
    def perform_segmentation(self, df_global):
        target_file = os.path.join(
            self.dirs["segments"],
            f"{self.args.SeqID}.segments.new.tsv",
        )

        def _run():
            log("Running 1D PELT Segmentation on log2_chrom_norm...")

            segments = segment_one_cell(
                meta=df_global,
                df_signals=df_global,
                signal_col="log2_chrom_norm",
                penalty=self.args.SegPenalty,
            )

            if segments is None or segments.empty:
                return pd.DataFrame()

            return assign_cn_state(
                segments,
                baseline_ploidy=self.args.BaselinePloidy,
                sex_tag=self._gender_tag,
            )
        return self._load_or_run(target_file, _run)

    # ─────────────────────────────────────────────────────────────
    # [모듈 7] 미세결실 탐지 (Micro)
    # ─────────────────────────────────────────────────────────────
    def perform_microdeletion(self, df_global):
        target_file = os.path.join(
            self.dirs["data"],
            f"{self.args.SeqID}.microdeletions.new.tsv",
        )

        def _run():
            log("Running rolling_micro_cnv_segmentation on normalized_count...")
            df_global["normalized_count"] = df_global["copy_number_signal"] / 2.0
            segments = rolling_micro_cnv_segmentation(
                df_global,
                signal_col="normalized_count",
                window=3,
                gain_threshold=1.3,
                loss_threshold=-0.7,
                min_bins=4,
                min_segment_bp=400_000,
                max_seg_mad=0.20,
            )
            if segments is None or segments.empty:
                return pd.DataFrame()

            return assign_cn_state(
                segments,
                baseline_ploidy=self.args.BaselinePloidy,
                sex_tag=self._gender_tag,
            )
        return self._load_or_run(target_file, _run)

    # ─────────────────────────────────────────────────────────────
    # [모듈 8] 임상 진단 리포트 생성
    # ─────────────────────────────────────────────────────────────
    def generate_clinical_report(self, df_global, segment_df, micro_df):
        """
        NIPT 임상 진단 리포트 생성 (마커 DB 교차 검증)
        """
        target_file = os.path.join(
            self.dirs["data"],
            f"{self.args.SeqID}.clinical_report.tsv"
        )

        def _run():
            log("Running Clinical Diagnosis against NIPT Marker Database...")

            marker_db_path = '/storage/home/jhkim/scripts/bioscript/manual/cbnipt_manual_cnv_call/src/cbnipt_cnv_caller/resources/gcx_nipt_target.tsv'
            if not os.path.exists(marker_db_path):
                log(f"[Error] Marker DB not found at {marker_db_path}")
                return pd.DataFrame()

            marker_df = pd.read_csv(marker_db_path, sep="\t")

            cnv_list = []
            if segment_df is not None and not segment_df.empty:
                cnv_list.append(segment_df)
            if micro_df is not None and not micro_df.empty:
                cnv_list.append(micro_df)

            combined_cnv = pd.concat(cnv_list, ignore_index=True) if cnv_list else pd.DataFrame()

            report_df = diagnose_clinical_markers(
                marker_df=marker_df,
                bin_df=df_global,
                cnv_df=combined_cnv,
                sex_tag=self._gender_tag
            )
            return report_df

        return self._load_or_run(target_file, _run)

    # ─────────────────────────────────────────────────────────────
    # [모듈 9] 시각화 (Plots)
    # ─────────────────────────────────────────────────────────────
    def generate_plots(self, df_global, summary, call_df, segment_df):
        gender_tag = self._gender_tag
        sample_id = self.args.SeqID
        plots_dir = self.dirs["plots"]

        log("Generating plots...")
        # [변경] plot_chromosome_overview 가 이제 bin distribution zoom까지 포함해서
        # 하나의 01_chromosome_overview.png 로 저장한다. (구 plot_bin_distribution 제거)
        plot_chromosome_overview(df_global, summary, call_df, gender_tag, sample_id, plots_dir)
        plot_score_heatmap(call_df, sample_id, gender_tag, plots_dir)
        plot_final_call(call_df, gender_tag, sample_id, plots_dir)

        if segment_df is not None and not segment_df.empty:
            plot_comprehensive_baf_logr_qc(
                df_global,
                segment_df,
                os.path.join(plots_dir, f"{sample_id}.comprehensive_view.png"),
                sample_id,
            )


# ═════════════════════════════════════════════════════════════
# [NEW] 여러 샘플을 돌린 뒤 종합 코호트 뷰를 뽑아내는 엔트리 함수
# ═════════════════════════════════════════════════════════════
def run_cohort_plot(run_dirs, out_path, sex_dict=None, score_col="final_score"):
    """
    run_dirs : dict {sample_id: OutDir}  -- 각 샘플을 CnvPipeline(args).run() 으로 이미 돌린 뒤의 OutDir
    out_path : 코호트 종합 plot 저장 경로 (.png)
    sex_dict : dict {sample_id: 'XX'/'XY'/'Unknown'} (선택)

    사용 예:
        run_dirs = {
            "SAMPLE_A": "/path/to/out/SAMPLE_A",
            "SAMPLE_B": "/path/to/out/SAMPLE_B",
        }
        run_cohort_plot(run_dirs, "/path/to/out/cohort_overview.png")
    """
    call_df_dict = load_cohort_call_dfs(run_dirs)
    plot_cohort_overview(call_df_dict, out_path, score_col=score_col, sex_dict=sex_dict)
    return call_df_dict