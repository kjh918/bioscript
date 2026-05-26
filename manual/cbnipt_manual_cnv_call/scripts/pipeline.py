import os
import pysam
import numpy as np
import pandas as pd
import sys

sys.path.append(os.path.dirname(__file__))  # 현재 디렉토리를 모듈 검색 경로에 추가
from utils import log, ensure_dir
from binning import generate_bins, annotate_bin_metadata, get_chromosomes, apply_final_filters
from processing import process_bam_to_coverage, apply_low_quality_filter, gc_correct_lowess
from normalization import normalize_by_chrom_with_sex, normalize_all_metrics_with_sex_log2
from visualization import plot_gc_correction, plot_genome_wide_cnv_with_gaps
from segmentation import segment_one_cell, assign_cn_state

# ---------------------------------------------------------
# Pipeline 내부에서 고속 연산을 처리하기 위한 최적화 헬퍼 함수들
# ---------------------------------------------------------
def _fetch_sites_with_handle(vcf_handle, chrom, start, end):
    sites = []
    try:
        for record in vcf_handle.fetch(chrom, start, end):
            kova = record.info.get("KOVA_AF")
            kova_af = kova[0] if kova and isinstance(kova, (list, tuple)) else (kova if kova else 0.0)
            gnomad = record.info.get("GNOMAD_AF")
            gnomad_af = gnomad[0] if gnomad and isinstance(gnomad, (list, tuple)) else (gnomad if gnomad else 0.0)
            pop_af = kova_af if kova_af > 0 else gnomad_af
            if pop_af > 0:
                sites.append({"chrom": record.chrom, "pos": record.pos, "ref": record.ref, "alt": record.alts[0], "pop_af": pop_af})
    except: pass
    return pd.DataFrame(sites)

def _summarize_and_classify_bin(fragment_df, site_df, bin_id, hetero_range):
    if fragment_df.empty: return pd.DataFrame(), pd.DataFrame()
    pos_stats = {
        row.pos: {
            "bin_id": bin_id, "chrom": row.chrom, "pos": row.pos, "ref": row.ref, "alt": row.alt, "pop_af": row.pop_af,
            "ref_depth": 0, "alt_depth": 0, "cis_support": 0, "trans_support": 0, "total_fragments": 0
        } for row in site_df.itertuples()
    }
    for frag in fragment_df.to_dict('records'):
        for pos, data in frag['obs_dict'].items():
            if pos not in pos_stats: continue
            stat = pos_stats[pos]
            stat["total_fragments"] += 1
            if data['base'] == stat['alt'].upper(): stat["alt_depth"] += 1
            elif data['base'] == stat['ref'].upper(): stat["ref_depth"] += 1
            if frag['is_cis_alt']: stat["cis_support"] += 1
            if frag['is_trans']: stat["trans_support"] += 1

    report_df = pd.DataFrame(pos_stats.values())
    if report_df.empty: return pd.DataFrame(), pd.DataFrame()
    
    # [MODIFIED] 1. 필터링 전 원천(Raw) 데이터 총량 합산 보존 (정규화 자산)
    raw_ref_sum = int(report_df["ref_depth"].sum())
    raw_alt_sum = int(report_df["alt_depth"].sum())
    raw_total_depth = raw_ref_sum + raw_alt_sum
    raw_total_fragments_sum = int(report_df["total_fragments"].sum())

    # 포지션별 기본 지표 연산
    report_df["observed_baf"] = (report_df["alt_depth"] / (report_df["ref_depth"] + report_df["alt_depth"])).fillna(0)
    
    # [MODIFIED] 2. QC 조건 만족 데이터 격리 관리 오퍼레이션
    # 포지션당 읽힌 독립 파편 수가 2개 이상인 클린 지점만 추출
    qc_mask = report_df['total_fragments'] >= 2
    qc_group = report_df[qc_mask]
    
    # QC 패스 순수 프래그먼트 양 격리 카운트 
    qc_pass_fragments = int(qc_group["total_fragments"].sum()) if not qc_group.empty else 0

    h_min, h_max = hetero_range
    report_df['group'] = 'Other'
    report_df.loc[(report_df['pop_af'] >= h_min) & (report_df['pop_af'] <= h_max), 'group'] = 'Hetero'
    report_df.loc[(report_df['pop_af'] <= 0.1) | (report_df['pop_af'] >= 0.9), 'group'] = 'Homo'

    # QC 통과 그룹 내에서만 안정적인 BAF 평균/메디안/개수 집계
    summary_grouped = qc_group.groupby('group')['observed_baf'].agg(['mean', 'median', 'count']) if not qc_group.empty else pd.DataFrame()
    
    res = {"bin_id": bin_id}
    for g in ['Hetero', 'Homo']:
        if not summary_grouped.empty and g in summary_grouped.index:
            res[f"{g.lower()}_baf_mean"] = summary_grouped.loc[g, 'mean']
            res[f"{g.lower()}_baf_median"] = summary_grouped.loc[g, 'median']
            res[f"{g.lower()}_sites_count"] = int(summary_grouped.loc[g, 'count'])
        else:
            res[f"{g.lower()}_baf_mean"] = res[f"{g.lower()}_baf_median"] = res[f"{g.lower()}_sites_count"] = 0
            
    # [MODIFIED] 3. 마스터 요약 프레임에 자산 연동 및 분모 제어 스케일링 적용
    total_trans = int(fragment_df['is_trans'].sum())
    total_cis = int(fragment_df['is_cis_alt'].sum())
    
    res["total_trans_fragments"] = total_trans
    res["total_cis_fragments"] = total_cis
    res["raw_ref_depth_sum"] = raw_ref_sum
    res["raw_alt_depth_sum"] = raw_alt_sum
    res["raw_total_depth"] = raw_total_depth
    res["raw_total_fragments_sum"] = raw_total_fragments_sum
    res["qc_pass_fragments"] = qc_pass_fragments  # 정규화 모듈용 명시적 격리 지표
    
    # QC 패스된 알짜배기 분모 풀 기준의 스케일링 비율 연산 (Double Counting 오류 차단)
    denom = qc_pass_fragments if qc_pass_fragments > 0 else 1
    report_df['TER'] = (report_df['trans_support'] / denom).fillna(0)
    report_df['CER'] = (report_df['cis_support'] / denom).fillna(0)
    
    return report_df, pd.DataFrame([res])

def _process_bin_with_handles(bam_handle, bin_info, site_df, min_mapq, min_baseq, hetero_range):
    chrom, start, end = bin_info['chrom'], bin_info['start'], bin_info['end']
    bin_sites = site_df[(site_df['chrom'] == chrom) & (site_df['pos'] >  start) & (site_df['pos'] < end)]
    if bin_sites.empty: return pd.DataFrame(), pd.DataFrame()
    site_lookup = {row.pos: {"ref": row.ref, "alt": row.alt} for row in bin_sites.itertuples()}

    fragments = {}
    for read in bam_handle.fetch(chrom, start, end):
        if read.is_unmapped or read.is_duplicate or read.is_secondary or read.is_supplementary: continue
        if read.mapping_quality < min_mapq: continue
        qname = read.query_name
        obs = {}
        for q_pos, r_pos in read.get_aligned_pairs(matches_only=True):
            g_pos = r_pos + 1
            if g_pos in site_lookup:
                qual = read.query_qualities[q_pos]
                if qual < min_baseq: continue
                obs[g_pos] = {'base': read.query_sequence[q_pos].upper(), 'qual': qual}

        if not obs: continue
        if qname not in fragments: fragments[qname] = obs
        else:
            for g_pos, data in obs.items():
                if g_pos not in fragments[qname] or data['qual'] > fragments[qname][g_pos]['qual']:
                    fragments[qname][g_pos] = data

    frag_list = []
    for qname, obs in fragments.items():
        status_list = []
        for p, d in obs.items():
            site = site_lookup[p]
            if d['base'] == site['alt'].upper(): status_list.append("ALT")
            elif d['base'] == site['ref'].upper(): status_list.append("REF")
        if status_list:
            n_alt, n_total = status_list.count("ALT"), len(status_list)
            frag_list.append({
                "obs_dict": obs,
                "is_cis_alt": (n_total >= 2 and n_alt == n_total),
                "is_trans": (n_total >= 2 and 0 < n_alt < n_total)
            })
    return _summarize_and_classify_bin(pd.DataFrame(frag_list), bin_sites, f"{chrom}:{start}-{end}", hetero_range)


class CnvPipeline:
    def __init__(self, args):
        self.args = args
        self.dirs = {k: os.path.join(args.OutDir, k) for k in ["meta", "qc", "data", "plots", "segments", "baf"]}
        for d in self.dirs.values(): ensure_dir(d)
    
    def run(self):
        log(f"--- Execution Start: {self.args.SeqID} ---")
        
        # ---------------------------------------------------------
        # [Step 1-1] Binning & Annotation
        # ---------------------------------------------------------
        chroms = get_chromosomes(fasta_path=self.args.ReferenceFasta, include_sex=self.args.IncludeSexChrom)
        bins = generate_bins(fasta_path=self.args.ReferenceFasta, bin_size=self.args.BinSize, chromosomes=chroms)
        bins_annotated = annotate_bin_metadata(bins, fasta_path=self.args.ReferenceFasta, mappability_bw=self.args.MappabilityBW)
        bins_annotated.to_csv(os.path.join(self.dirs["meta"], f"{self.args.SeqID}.bins_annotated.tsv"), sep="\t", index=False)
        
        filtered_bins = apply_final_filters(bins_annotated, self.args)
        filtered_bins.to_csv(os.path.join(self.dirs["data"], f"{self.args.SeqID}.bins_annotated_filtered.tsv"), sep="\t", index=False)
        
        # ---------------------------------------------------------
        # [Step 1-2] BAM Coverage Processing & Quality Filters (순서 전진 배치)
        # [REASON] 노이즈 빈을 먼저 제거해야 후속 BAF/TER 루프 연산량이 급감하고 무결성이 보장됨
        # ---------------------------------------------------------
        log("Processing BAM Coverage & Low-Quality Bin Filtering first...")
        df_with_reads, filter_stats = process_bam_to_coverage(
            bam_path=self.args.BamPath, bins_df=filtered_bins, min_mapq=self.args.MinMapQ
        )
        df_with_reads.to_csv(os.path.join(self.dirs["qc"], f"{self.args.SeqID}.raw_data.tsv"), sep="\t", index=False)
        
        df_clean, filter_stats_df = apply_low_quality_filter(
            df_with_reads, filter_stats, min_depth=self.args.MinDepth, min_coverage=self.args.MinCoverage
        )
        df_clean.to_csv(os.path.join(self.dirs["data"], f"{self.args.SeqID}.bins_with_reads_filtered.tsv"), sep="\t", index=False)
        filter_stats_df.to_csv(os.path.join(self.dirs["qc"], f"{self.args.SeqID}.filter_stats.tsv"), sep="\t", index=False)

        # ---------------------------------------------------------
        # [Step 1-3] Genetic Evidence Analysis (정제된 df_clean 기반으로 수행)
        # ---------------------------------------------------------
        log("Starting High-Speed Genetic Evidence Loop on Clean Bins...")
        genetic_bin_summaries = []
        all_position_reports = []

        # BAM 핸들을 오픈하여 정제된 df_clean의 빈 영역만 고속 순회
        with pysam.AlignmentFile(self.args.BamPath, "rb", threads=4) as bam_handle:
            for chrom, group in df_clean.groupby("chrom"):
                vcf_path = f"/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Resources/reference/KOVA_v7/kova_sites_vcf/KOVA_v7_SNPOnly_QCSites_af.{chrom}.vcf.gz"
                
                if not os.path.exists(vcf_path):
                    log(f"Warning: VCF not found for {chrom}, skipping haplotype analysis.")
                    continue
                    
                log(f"Processing Genetic Evidence for clean {chrom}...")
                
                with pysam.VariantFile(vcf_path) as vcf_handle:
                    for r in group.itertuples():
                        # 1. 해당 영역 SNP 추출
                        site_df = _fetch_sites_with_handle(vcf_handle, r.chrom, r.start, r.end)
                        if site_df.empty: continue

                        # 2. 프래그먼트 분석 및 정밀 BAF/TER 추출
                        report_df, summary_df = _process_bin_with_handles(
                            bam_handle, 
                            {"chrom": r.chrom, "start": r.start, "end": r.end}, 
                            site_df, 
                            min_mapq=self.args.MinMapQ, 
                            min_baseq=20,
                            hetero_range=(0.4, 0.6)
                        )
                        if not report_df.empty and not summary_df.empty:
                            all_position_reports.append(report_df)
                            
                            summary_dict = summary_df.iloc[0].to_dict()
                            summary_dict["chrom"] = r.chrom
                            summary_dict["start"] = r.start
                            summary_dict["end"] = r.end
                            genetic_bin_summaries.append(summary_dict)

        # ---------------------------------------------------------
        # [Step 1-4] 유전적 증거 데이터 병합, 포지션 정규화 및 저장
        # ---------------------------------------------------------
        if genetic_bin_summaries:
            genetic_evidence_df = pd.DataFrame(genetic_bin_summaries)
            # [CRITICAL] 퀄리티 통과한 df_clean에 유전 증거를 레프트 머지합니다.
            df_clean = df_clean.merge(genetic_evidence_df, on=["chrom", "start", "end"], how="left")
            genetic_evidence_df.to_csv(os.path.join(self.dirs["data"], f"{self.args.SeqID}.genetic_evidence.tsv"), sep="\t", index=False)
            
        if all_position_reports:
            full_pos_report = pd.concat(all_position_reports, ignore_index=True)
            full_pos_report.to_csv(os.path.join(self.dirs["baf"], f"{self.args.SeqID}.position_baf_detail.tsv"), sep="\t", index=False)
            
            # 포지션 레벨 이원화 정규화 엔진 가동
            log("Running Position-Level Dual-Normalization for detailed BAF/TER plotting...")
            pos_norm_df = full_pos_report.copy()
            try:
                pos_norm_df["pos_raw_depth"] = (pos_norm_df["ref_depth"] + pos_norm_df["alt_depth"]).replace(0, 1)
                pos_norm_df["pos_qc_fragments"] = pos_norm_df["total_fragments"].apply(lambda x: x if x >= 2 else 1)
                
                autosomes = pos_norm_df[~pos_norm_df["chrom"].isin(["chrX", "chrY"])]
                if not autosomes.empty:
                    pos_norm_df["density_trans"] = pos_norm_df["trans_support"] / pos_norm_df["pos_qc_fragments"]
                    base_trans = max(autosomes["trans_support"].sum() / autosomes["total_fragments"].apply(lambda x: x if x >= 2 else 1).sum(), 1e-8)
                    pos_norm_df["trans_log2_norm"] = np.log2((pos_norm_df["density_trans"] + 1e-10) / base_trans)
                    
                    chrY_data = pos_norm_df[pos_norm_df["chrom"] == "chrY"]
                    is_male = len(chrY_data) > 5 if not chrY_data.empty else False
                    
                    for chrom, group in pos_norm_df.groupby("chrom"):
                        target_offset = 0.0
                        if chrom == "chrX" and is_male: target_offset = -1.0
                        elif chrom == "chrY" and is_male: target_offset = -1.0
                        elif chrom == "chrY" and not is_male: continue
                        pos_norm_df.loc[pos_norm_df["chrom"] == chrom, "trans_log2_norm"] += target_offset
                
                normalized_pos_path = os.path.join(self.dirs["baf"], f"{self.args.SeqID}.position_baf_normalized.tsv")
                pos_norm_df.to_csv(normalized_pos_path, sep="\t", index=False)
                log(f"[+] Normalized Position-Level BAF report saved to: {normalized_pos_path}")
            except Exception as pos_err:
                log(f"Warning: Position-level normalization failed due to: {pos_err}.")

        # ---------------------------------------------------------
        # [Step 3] GC 보정 (유전 증거 정보가 탑재된 df_clean 기반으로 연산 연속 진행)
        # ---------------------------------------------------------
        df_gc, gc_stats = gc_correct_lowess(df_clean, frac=self.args.LowessFrac)
        plot_gc_correction(gc_stats, os.path.join(self.dirs["plots"], "gc_correction.png"))
        df_gc.to_csv(os.path.join(self.dirs["data"], f"{self.args.SeqID}.gc_corrected.tsv"), sep="\t", index=False)
        
        log("Calculating Global Normalization & Chromosome Summary...")
        df_global = normalize_by_chrom_with_sex(df_gc, value_col="log2_corrected")
        df_global.to_csv(os.path.join(self.dirs["data"], f"{self.args.SeqID}.normalized.tsv"), sep="\t", index=False)

        # [Step 4] 세그멘테이션 실행
        segments = segment_one_cell(
            df_global, 
            df_global["log2_chrom_norm"], 
            self.args.SegPenalty
        )
        
        # [Step 5] 복제수(Copy Number) 할당
        segments = assign_cn_state(
            segments, 
            baseline_ploidy=self.args.BaselinePloidy
        )
        segments.to_csv(os.path.join(self.dirs["data"], f"{self.args.SeqID}.segments.tsv"), sep="\t", index=False)

        # [Step 6] 게놈 와이드 시각화
        plot_genome_wide_cnv_with_gaps(
            bins_df=df_global, 
            segments_df=segments, 
            output_path=os.path.join(self.dirs["plots"], f"{self.args.SeqID}.genome_wide.png"),
            run_id=self.args.SeqID
        )
        log(f"--- [COMPLETE] Results saved to {self.args.OutDir} ---")