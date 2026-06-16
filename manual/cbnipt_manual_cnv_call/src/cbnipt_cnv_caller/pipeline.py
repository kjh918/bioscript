import os
import pysam
import numpy as np
import pandas as pd
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

from utils import log, ensure_dir, chrom_key, sort_chroms
from processing import gc_correct_lowess, _parallel_chrom_worker
from normalization import normalize_by_chrom_with_sex, apply_qc_bin
from visualization import (
    plot_comprehensive_baf_logr_qc,
    plot_chromosome_overview,
    plot_bin_distribution,
    plot_score_heatmap,
    plot_final_call,
    plot_gc_correction,
)
from segmentation import assign_cn_state, segment_one_cell
from classification import compute_chrom_summary, analyze_all_chromosomes


class CnvPipeline:
    def __init__(self, args):
        self.args = args
        self.dirs = {k: os.path.join(args.OutDir, k)
                     for k in ["meta", "qc", "data", "plots", "segments", "baf"]}
        for d in self.dirs.values():
            ensure_dir(d)

        # normalize 단계에서 채워지는 성별/비율 정보
        self._gender_tag = "Unknown"
        self._x_ratio    = 0.0
        self._y_ratio    = 0.0

    # ─────────────────────────────────────────────────────────────
    # 내부 유틸
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
    # [모듈 0] Bin 로드
    # ─────────────────────────────────────────────────────────────
    def load_target_bins(self):
        bin_file = self.args.AnnotatedBins
        if not os.path.exists(bin_file):
            raise FileNotFoundError(f"Cannot find bin file: {bin_file}")
        log(f"Loading pre-computed BED.GZ bins: {bin_file}")
        df_bins = pd.read_csv(bin_file, sep="\t", compression="gzip", comment=None)
        if "#chrom" in df_bins.columns:
            df_bins.rename(columns={"#chrom": "chrom"}, inplace=True)
        return df_bins

    # ─────────────────────────────────────────────────────────────
    # [모듈 1] One-Pass 병렬 추출 (BAM + VCF)
    # ─────────────────────────────────────────────────────────────
    def extract_genetic_evidence(self, df_bins):
        target_file = os.path.join(self.dirs["data"],
                                   f"{self.args.SeqID}.bins_with_evidence.tsv")

        def _run():
            tmp_dir = os.path.join(self.dirs["data"], "tmp_workers")
            ensure_dir(tmp_dir)

            chrom_groups = [group.copy() for _, group in df_bins.groupby("chrom")]
            max_workers  = getattr(self.args, "Threads", 4)
            log(f"Parallel processing with {max_workers} threads...")

            valid_pos_files, valid_sum_files, valid_mut_files = [], [], []
            thresholds = {
                "homo_max":   getattr(self.args, "HomoMax",   0.1),
                "hetero_min": getattr(self.args, "HeteroMin", 0.4),
                "hetero_max": getattr(self.args, "HeteroMax", 0.6),
                "homo_min":   getattr(self.args, "HomoMin",   0.9),
            }

            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                futures = {
                    executor.submit(
                        _parallel_chrom_worker,
                        group, self.args.BamPath, self.args.VcfFile,
                        self.args.MinMapQ, thresholds, tmp_dir
                    ): group["chrom"].iloc[0]
                    for group in chrom_groups
                }
                for future in as_completed(futures):
                    pos_path, sum_path, mut_path = future.result()
                    if pos_path: valid_pos_files.append(pos_path)
                    if sum_path: valid_sum_files.append(sum_path)
                    if mut_path: valid_mut_files.append(mut_path)

            log("Merging worker results...")

            if valid_pos_files:
                full_pos = pd.concat(
                    [pd.read_parquet(f) for f in valid_pos_files], ignore_index=True)
                baf_out = os.path.join(self.dirs["baf"],
                                       f"{self.args.SeqID}_position_baf")
                full_pos.to_parquet(baf_out, engine="pyarrow",
                                    partition_cols=["chrom"], index=False)

            if valid_mut_files:
                full_mut = pd.concat(
                    [pd.read_parquet(f) for f in valid_mut_files], ignore_index=True)
                full_mut.to_csv(
                    os.path.join(self.dirs["data"],
                                 f"{self.args.SeqID}.mutation_signatures.tsv"),
                    sep="\t", index=False)

            evidence_df = pd.concat(
                [pd.read_parquet(f) for f in valid_sum_files], ignore_index=True)
            mean_density = evidence_df["total_sites"].mean()
            evidence_df["Informative_OE_Ratio"]  = (
                evidence_df["total_sites"] / (mean_density + 1e-9))
            evidence_df["Informative_Norm_Log2"] = np.log2(
                evidence_df["Informative_OE_Ratio"] + 1e-9)

            for f in valid_pos_files + valid_sum_files + valid_mut_files:
                if os.path.exists(f): os.remove(f)
            try:
                os.rmdir(tmp_dir)
            except OSError:
                pass

            df_merged = df_bins.merge(evidence_df, on=["chrom", "start", "end"], how="left")
            df_merged["raw_count"]     = df_merged["raw_count"].fillna(0)
            df_merged["breadth_ratio"] = df_merged["breadth_ratio"].fillna(0.0)
            return df_merged

        return self._load_or_run(target_file, _run)

    # ─────────────────────────────────────────────────────────────
    # [모듈 2] GC 보정 + 성별-인식 LOO 정규화
    #
    # normalize_by_chrom_with_sex() 가 수행하는 핵심 보정:
    #   log2_chrom_norm = log2(obs / ref_median) - expected_log2
    #
    # 따라서 정상 염색체(포함 male chrX/Y)는 모두 0.0 근방에 놓입니다.
    # classification 의 robust-Z 는 이 0.0 기준으로만 계산합니다.
    # ─────────────────────────────────────────────────────────────
    def correct_gc_and_normalize(self, df_merged):
        target_file = os.path.join(self.dirs["data"],
                                   f"{self.args.SeqID}.normalized.tsv")

        def _run():
            log("Applying GC Correction (LOWESS)...")
            df_gc, gc_stats = gc_correct_lowess(df_merged, frac=self.args.LowessFrac)
            plot_gc_correction(
                gc_stats,
                os.path.join(self.dirs["qc"], f"{self.args.SeqID}.gc_correction.png"))

            log("Normalizing (LOO, sex-aware, expected_log2 subtracted)...")
            # 반환: (df, gender_tag, y_ratio)
            df_norm, gender_tag, x_ratio, y_ratio = normalize_by_chrom_with_sex(
                df_gc, value_col="log2_corrected")

            self._gender_tag = gender_tag
            self._y_ratio    = y_ratio
            # x_ratio: raw_count 기준 (성별 컷오프 판단용)
            self._x_ratio    = x_ratio # self._calc_x_ratio(df_gc, value_col="log2_corrected")

            return df_norm

        return self._load_or_run(target_file, _run)

    # ─────────────────────────────────────────────────────────────
    # 공통 헬퍼: x_ratio 계산 (raw_count 기준)
    # ─────────────────────────────────────────────────────────────
    @staticmethod
    def _calc_x_ratio(df, value_col="raw_count"):
        """상염색체 median 대비 chrX median 비율 (성염색체 컷오프 판단용)."""
        auto = df[~df["chrom"].isin(["chrX", "chrY"])][value_col]
        auto_med = auto.median()
        if auto_med == 0:
            return 0.0
        x_med = df[df["chrom"] == "chrX"][value_col].median()
        return float(x_med / auto_med)

    # ─────────────────────────────────────────────────────────────
    # [모듈 3] Classification (염색체 수 이상 탐지)
    # ─────────────────────────────────────────────────────────────
    def classify_chromosomes(self, df_global):
        """
        summary  → analyze_all_chromosomes 실행 후 결과 저장.
        robust-Z 는 log2_chrom_norm(≈0 기준) 기반으로만 계산됩니다.
        """
        target_file = os.path.join(self.dirs["data"],
                                   f"{self.args.SeqID}.chrom_calls.tsv")

        log("Computing chromosome-level summary statistics...")
        summary = compute_chrom_summary(df_global)

        log(f"Classification — Sex={self._gender_tag}, "
            f"xRatio={self._x_ratio:.3f}, yRatio={self._y_ratio:.3f}")
        call_df = analyze_all_chromosomes(
            summary, self._gender_tag, self._x_ratio, self._y_ratio)

        if not os.path.exists(target_file):
            call_df.to_csv(target_file, sep="\t", index=False)

        return summary, call_df

    # ─────────────────────────────────────────────────────────────
    # [모듈 4] Segmentation
    # ─────────────────────────────────────────────────────────────
    def perform_segmentation(self, df_global):
        target_file = os.path.join(self.dirs["segments"],
                                   f"{self.args.SeqID}.segments.tsv")

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
    # [모듈 5] 시각화 (plots/ 에 저장)
    # ─────────────────────────────────────────────────────────────
    def generate_plots(self, df_global, summary, call_df, segment_df):
        gender_tag = self._gender_tag
        sample_id  = self.args.SeqID
        plots_dir  = self.dirs["plots"]

        log("Generating plots...")
        plot_chromosome_overview(df_global, summary, call_df,
                                 gender_tag, sample_id, plots_dir)
        plot_bin_distribution(df_global, call_df, sample_id, plots_dir)
        plot_score_heatmap(call_df, sample_id, gender_tag, plots_dir)
        plot_final_call(call_df, gender_tag, sample_id, plots_dir)

        if segment_df is not None and not segment_df.empty:
            plot_comprehensive_baf_logr_qc(
                df_global,
                segment_df,
                os.path.join(plots_dir,
                             f"{sample_id}.comprehensive_view.png"),
                sample_id,
            )

    # ─────────────────────────────────────────────────────────────
    # [모듈 6] 콘솔 출력
    # ─────────────────────────────────────────────────────────────
    def _print_call_table(self, call_df):
        ordered = (call_df
                   .set_index("chrom")
                   .reindex(sort_chroms(call_df["chrom"].tolist()))
                   .reset_index())
        print(f"\n{'Chrom':<8} {'Log2':>7} {'RobustZ':>8} {'Final':>7}  Call")
        print("-" * 70)
        for _, row in ordered.iterrows():
            flag = ("⚠  " if row["call"] == "ABNORMAL"
                    else "△  " if row["call"] == "SUSPICIOUS"
                    else "   ")
            log2 = f"{row['log2fc']:+.3f}" if not np.isnan(row["log2fc"]) else "   N/A"
            rz   = f"{row['robust_z']:+.2f}" if not np.isnan(row["robust_z"]) else "  N/A"
            print(f"{row['chrom']:<8} {log2:>7} {rz:>8} {row['final_score']:>7.3f}  "
                  f"{flag}{row['call']:<12} {row.get('detail', '')}")

    # ═════════════════════════════════════════════════════════════
    # [실행 루프 B] 기존 normalized TSV → Classification + 시각화
    #
    # 입력 TSV 는 raw_count 컬럼을 포함해야 합니다.
    # normalize_by_chrom_with_sex() 로 log2_chrom_norm 을 생성(0.0 기준)한 뒤
    # 나머지 흐름은 run() 과 완전히 동일합니다.
    # ═════════════════════════════════════════════════════════════
    def run_from_normalized(self, normalized_tsv: str):
        """
        Parameters
        ----------
        normalized_tsv : str
            raw_count 컬럼을 포함한 TSV 파일 경로. (보통 bins_with_evidence.tsv)
            여기서부터 QC 필터 -> GC 보정 -> 정규화 -> 분류 -> 시각화의 흐름을 다시 탑니다.
        """
        log("=" * 60)
        log(f"  From Normalized TSV: {self.args.SeqID}")
        log("=" * 60)

        if not os.path.exists(normalized_tsv):
            raise FileNotFoundError(f"Normalized TSV not found: {normalized_tsv}")

        log(f"Loading: {os.path.basename(normalized_tsv)}")
        raw_df = pd.read_csv(normalized_tsv, sep="\t")

        # ── 1. QC 필터 (메인 run() 과 동일하게 성염색체 보호 필터 적용) ──
        log("Applying QC Bin Filter...")
        df_f = apply_qc_bin(raw_df)

        # ── 2. GC 보정 (LOWESS) ──────────────────────────────────
        log("Applying GC Correction (LOWESS) before normalization...")
        # [CRITICAL FIX] raw_count를 바로 쓰지 않고 반드시 GC 보정을 거칩니다.
        df_gc, _ = gc_correct_lowess(df_f, frac=self.args.LowessFrac)

        # ── 3. 정규화 (LOO, sex-aware, expected_log2 subtracted) ──
        log("Normalizing (LOO, sex-aware, expected_log2 subtracted)...")
        # [CRITICAL FIX] raw_count가 아닌, 보정된 "log2_corrected"를 사용합니다!
        df_norm, gender_tag, y_ratio = normalize_by_chrom_with_sex(
            df_gc, value_col="log2_corrected")

        self._gender_tag = gender_tag
        self._y_ratio    = y_ratio
        # X 염색체 비율도 GC 보정된 값을 기준으로 계산해야 정확합니다.
        self._x_ratio    = self._calc_x_ratio(df_gc, value_col="log2_corrected")

        log(f"Sex={gender_tag}  xRatio={self._x_ratio:.3f}  yRatio={y_ratio:.3f}")

        # 정규화 결과 저장 (캐싱)
        std_path = os.path.join(self.dirs["data"],
                                f"{self.args.SeqID}.standardized.tsv")
        df_norm.to_csv(std_path, sep="\t", index=False)
        log(f"Standardized TSV saved: {os.path.basename(std_path)}")

        # ── 4. Classification / Segmentation / 시각화 (run() 과 완벽 공유) ──
        summary, call_df = self.classify_chromosomes(df_norm)
        self._print_call_table(call_df)

        segment_df = self.perform_segmentation(df_norm)
        self.generate_plots(df_norm, summary, call_df, segment_df)

        log("=" * 60)
        log(f"  run_from_normalized Complete: {self.args.SeqID}")
        log(f"  Output: {self.dirs['plots']}")
        log("=" * 60)

        return {"df_global": df_norm, "summary": summary,
                "call_df": call_df, "segment_df": segment_df}

                
    # ═════════════════════════════════════════════════════════════
    # [실행 루프 A] BAM → 전체 파이프라인
    # ═════════════════════════════════════════════════════════════
    def run(self):
        log("=" * 60)
        log(f"  Full Pipeline: {self.args.SeqID}")
        log("=" * 60)

        df_bins            = self.load_target_bins()
        df_merged_evidence = self.extract_genetic_evidence(df_bins)
        df_filtered        = apply_qc_bin(df_merged_evidence)

        # GC 보정 + LOO 정규화 → log2_chrom_norm (0.0 기준)
        # 동시에 self._gender_tag / _x_ratio / _y_ratio 세팅
        df_global = self.correct_gc_and_normalize(df_filtered)

        # Classification (robust-Z 기준)
        summary, call_df = self.classify_chromosomes(df_global)
        self._print_call_table(call_df)

        # Segmentation (동일 log2_chrom_norm / gender_tag 사용)
        segment_df = self.perform_segmentation(df_global)

        self.generate_plots(df_global, summary, call_df, segment_df)

        log("=" * 60)
        log(f"  Pipeline Complete: {self.args.SeqID}")
        log("=" * 60)

        return {"df_global": df_global, "summary": summary,
                "call_df": call_df, "segment_df": segment_df}
