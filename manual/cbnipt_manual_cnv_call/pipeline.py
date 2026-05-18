import os
import pandas as pd
import sys
sys.path.append(os.path.dirname(__file__))  # 현재 디렉토리를 모듈 검색 경로에 추가
from utils import log, ensure_dir
from binning import generate_bins, annotate_bin_metadata, get_chromosomes, apply_final_filters
from processing import  process_bam_to_coverage, apply_low_quality_filter, gc_correct_lowess


from normalization import normalize_by_chrom_with_sex
from visualization import plot_gc_correction, plot_genome_wide_cnv
from segmentation import segment_one_cell, assign_cn_state
from ado import process_fragments_for_baf, fetch_sites_from_vcf_gz, summarize_by_position
from utils import log

class CnvPipeline:
    def __init__(self, args):
        self.args = args
        self.dirs = {k: os.path.join(args.OutDir, k) for k in ["meta", "qc", "data", "plots", "segments"]}
        for d in self.dirs.values(): ensure_dir(d)
    
    def run(self):
        log(f"--- Execution Start: {self.args.RunID} ---")
        
        # [Step 1-1] Binning & Annotation
        chroms = get_chromosomes(fasta_path=self.args.ReferenceFasta, include_sex=self.args.IncludeSexChrom)
        bins = generate_bins(fasta_path=self.args.ReferenceFasta, bin_size=self.args.BinSize, chromosomes=chroms)
        bins_annotated = annotate_bin_metadata(bins, fasta_path=self.args.ReferenceFasta, mappability_bw=self.args.MappabilityBW)
        bins_annotated.to_csv(os.path.join(self.dirs["meta"], f"{self.args.RunID}.bins_annotated.tsv"), sep="\t", index=False)
        
        filtered_bins = apply_final_filters(bins_annotated, self.args)
        filtered_bins.to_csv(os.path.join(self.dirs["data"], f"{self.args.RunID}.bins_annotated_filtered.tsv"), sep="\t", index=False)
        
        genetic_bin_summaries = []
        all_position_reports = []

        # [CRITICAL] 염색체별로 그룹화하여 VCF 파일 접근 최적화
        for chrom, group in filtered_bins.groupby("chrom"):
            # 염색체별 VCF 경로 설정
            vcf_path = f"/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Resources/reference/KOVA_v7/kova_sites_vcf/KOVA_v7_SNPOnly_QCSites_af.{chrom}.vcf.gz"
            
            if not os.path.exists(vcf_path):
                log(f"Warning: VCF not found for {chrom}, skipping haplotype analysis.")
                continue
                
            log(f"Processing Genetic Evidence for {chrom}...")
            
            for r in group.itertuples():
                # 1. 해당 영역 SNP 추출
                site_df = fetch_sites_from_vcf_gz(vcf_path, r.chrom, r.start, r.end)
                if site_df.empty: continue

                # 2. 프래그먼트 분석 (TER/CER 계산의 기초 데이터)
                fragment_df, _ = process_fragments_for_baf(
                    self.args.BamPath, 
                    {"chrom": r.chrom, "start": r.start, "end": r.end}, 
                    site_df, 
                    min_mapq=30, 
                    min_baseq=20
                )
                print(fragment_df)
                
                if fragment_df.empty: continue

        # [Step 1-3] 유전적 증거 데이터 병합 및 저장
        if genetic_bin_summaries:
            genetic_evidence_df = pd.DataFrame(genetic_bin_summaries)
            # 메인 데이터프레임에 병합
            filtered_bins = filtered_bins.merge(genetic_evidence_df, on=["chrom", "start", "end"], how="left")
            # QC용 데이터 저장
            genetic_evidence_df.to_csv(os.path.join(self.dirs["data"], f"{self.args.RunID}.genetic_evidence.tsv"), sep="\t", index=False)
            
        if all_position_reports:
            full_pos_report = pd.concat(all_position_reports)
            full_pos_report.to_csv(os.path.join(self.dirs["baf"], f"{self.args.RunID}.position_baf_detail.tsv"), sep="\t", index=False)


        # [Step 1-2] BAM 처리 및 3중 필터링 적용
        df_with_reads, filter_stats = process_bam_to_coverage(
            bam_path=self.args.BamPath, bins_df=filtered_bins, min_mapq=self.args.MinMapQ
        )
        df_with_reads.to_csv(os.path.join(self.dirs["qc"], f"{self.args.RunID}.raw_data.tsv"), sep="\t", index=False)
        df_clean, filter_stats_df = apply_low_quality_filter(df_with_reads, filter_stats, min_depth=self.args.MinDepth, min_coverage=self.args.MinCoverage)
        df_clean.to_csv(os.path.join(self.dirs["data"], f"{self.args.RunID}.bins_with_reads_filtered.tsv"), sep="\t", index=False)
        filter_stats_df.to_csv(os.path.join(self.dirs["qc"], f"{self.args.RunID}.filter_stats.tsv"), sep="\t", index=False)

        # [Step 3] GC 보정 (processing.py)
        df_gc, gc_stats = gc_correct_lowess(df_clean, frac=self.args.LowessFrac)
        plot_gc_correction(gc_stats, os.path.join(self.dirs["plots"], "gc_correction.png"))
        df_gc.to_csv(os.path.join(self.dirs["data"], f"{self.args.RunID}.gc_corrected.tsv"), sep="\t", index=False)
        
        log("Calculating Global Normalization & Chromosome Summary...")
        df_global = normalize_by_chrom_with_sex(df_gc, value_col="log2_corrected")
        df_global.to_csv(os.path.join(self.dirs["data"], f"{self.args.RunID}.normalized.tsv"), sep="\t", index=False)
        print(df_global)

        ## [Step 6-7] 세그멘테이션 및 저장
        segments = segment_one_cell(
            df_global, 
            df_global["log2_chrom_norm"], 
            self.args.SegPenalty
        )
        
        # [Step 6] 복제수(Copy Number) 할당
        # [MODIFIED] BaselinePloidy를 기반으로 정수형 CN 계산
        segments = assign_cn_state(
            segments, 
            baseline_ploidy=self.args.BaselinePloidy
        )
        segments.to_csv(os.path.join(self.dirs["data"], f"{self.args.RunID}.segments.tsv"), sep="\t", index=False)

        # [Step 7] 게놈 와이드 시각화
        # [MODIFIED] Bin 데이터와 Segment 데이터를 결합하여 리포트용 그래프 생성
        plot_genome_wide_cnv(
            bins_df=df_global, 
            segments_df=segments, 
            output_path=os.path.join(self.dirs["plots"], f"{self.args.RunID}.genome_wide.png"),
            run_id=self.args.RunID
        )
        log(f"--- [COMPLETE] Results saved to {self.args.OutDir} ---")