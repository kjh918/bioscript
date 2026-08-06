import numpy as np
import pandas as pd
import pysam
import os
from utils import log
from statsmodels.nonparametric.smoothers_lowess import lowess

from utils import chrom_key
from rules import CFG

_BIN_EXTRACTION_CFG = CFG["BIN_EXTRACTION"]

# -------------------------------------------------------------------------
# [데이터 스키마 정의] Breadth와 Midpoint Read Count 추가
# -------------------------------------------------------------------------
BIN_SUMMARY_SCHEMA = {
    "bin_id": "object",
    "raw_count": "int32",               "breadth_ratio": "float32",

    "total_sites": "int32",
    "ref_sum": "int32",                 "alt_sum": "int32",         "other_sum": "int32",
    "total_depth": "int32",             "bin_BAF": "float32",

    "pop_hetero_count": "int32",        "pop_homo_count": "int32",

    "hetero_like_count": "int32",       "imbalance_count": "int32", "homo_like_count": "int32",
    "hetero_like_rate": "float32",      "imbalance_rate": "float32","homo_like_rate": "float32",
    "MAD_BAF": "float32",

    # 노이즈 척도
    "on_target_noise_rate": "float32",  "off_target_noise_rate": "float32",

    "total_fragments_sum": "int32",     "qc_pass_fragments": "int32",
    "total_trans_fragments": "int32",   "total_cis_fragments": "int32",

    "raw_bin_TER": "float32",           "raw_bin_CER": "float32",
    "adj_bin_TER": "float32",           "adj_bin_CER": "float32"
}


def process_bin_with_handles(bam_handle, bin_info, site_df, min_mapq, min_baseq, thresholds):
    chrom, start, end = bin_info['chrom'], bin_info['start'], bin_info['end']
    bin_len = end - start

    if site_df.empty or 'chrom' not in site_df.columns:
        bin_sites = pd.DataFrame(columns=["chrom", "pos", "ref", "alt", "pop_af"])
    else:
        bin_sites = site_df[(site_df['chrom'] == chrom) & (site_df['pos'] > start) & (site_df['pos'] < end)]

    site_lookup = {row.pos: {"ref": row.ref, "alt": row.alt} for row in bin_sites.itertuples()} if not bin_sites.empty else {}

    fragments = {}
    raw_count = 0
    coverage_mask = np.zeros(bin_len, dtype=bool)

    total_mapped_bases = 0
    total_nm_errors = 0

    for read in bam_handle.fetch(chrom, start, end):
        if read.is_unmapped or read.is_duplicate or read.is_secondary or read.is_supplementary: continue
        if read.mapping_quality < min_mapq: continue

        # 전체 매핑 염기 수와 NM 태그(에러 수) 수집
        total_mapped_bases += read.query_alignment_length
        if read.has_tag("NM"):
            total_nm_errors += read.get_tag("NM")

        ref_start = max(read.reference_start, start)
        ref_end = min(read.reference_end, end)
        if ref_start < ref_end:
            coverage_mask[ref_start - start : ref_end - start] = True

        mid = (read.reference_start + read.reference_end) // 2
        if start <= mid < end:
            raw_count += 1

        if not site_lookup: continue

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

    breadth_ratio = float(np.sum(coverage_mask) / bin_len)

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

    return summarize_and_classify_bin(pd.DataFrame(frag_list), bin_sites, f"{chrom}:{start}-{end}", raw_count, breadth_ratio, thresholds, total_mapped_bases, total_nm_errors)

def summarize_and_classify_bin(fragment_df, site_df, bin_id, raw_count, breadth_ratio, thresholds, total_mapped_bases, total_nm_errors):
    res = {k: (0 if 'int' in v else 0.0) for k, v in BIN_SUMMARY_SCHEMA.items()}
    res["bin_id"] = bin_id
    res["raw_count"] = raw_count
    res["breadth_ratio"] = breadth_ratio

    mut_dist_df = pd.DataFrame()

    if fragment_df.empty or site_df.empty:
        return pd.DataFrame(), pd.DataFrame([res]), mut_dist_df

    # 1. 포지션별 기본 통계 적재 (other_depth 추가)
    pos_stats = {
        row.pos: {
            "bin_id": bin_id, "chrom": row.chrom, "pos": row.pos, "ref": row.ref, "alt": row.alt, "pop_af": row.pop_af,
            "ref_depth": 0, "alt_depth": 0, "other_depth": 0,
            "cis_support": 0, "trans_support": 0, "total_fragments": 0
        } for row in site_df.itertuples()
    }

    for frag in fragment_df.to_dict('records'):
        for pos, data in frag['obs_dict'].items():
            if pos not in pos_stats: continue
            stat = pos_stats[pos]
            stat["total_fragments"] += 1
            if data['base'] == stat['alt'].upper(): stat["alt_depth"] += 1
            elif data['base'] == stat['ref'].upper(): stat["ref_depth"] += 1
            else: stat["other_depth"] += 1

            if frag['is_cis_alt']: stat["cis_support"] += 1
            if frag['is_trans']: stat["trans_support"] += 1

    report_df = pd.DataFrame(pos_stats.values())
    if report_df.empty: return pd.DataFrame(), pd.DataFrame([res]), mut_dist_df

    report_df["total_depth"] = report_df["ref_depth"] + report_df["alt_depth"] + report_df["other_depth"]
    report_df["BAF"] = (report_df["alt_depth"] / report_df["total_depth"]).fillna(0)

    cond_pop = [
        (report_df['pop_af'] >= thresholds['hetero_min']) & (report_df['pop_af'] <= thresholds['hetero_max']),
        (report_df['pop_af'] <= thresholds['homo_max']) | (report_df['pop_af'] >= thresholds['homo_min'])
    ]
    report_df['pop_class'] = np.select(cond_pop, ['pop_hetero_informative', 'pop_homo_like'], default='intermediate')

    cond_baf = [
        (report_df['BAF'] >= thresholds['hetero_min']) & (report_df['BAF'] <= thresholds['hetero_max']),
        (report_df['BAF'] <= thresholds['homo_max']) | (report_df['BAF'] >= thresholds['homo_min'])
    ]
    report_df['baf_class'] = np.select(cond_baf, ['hetero_like', 'homo_like'], default='imbalance')

    pop_counts = report_df['pop_class'].value_counts()
    res["pop_hetero_count"] = int(pop_counts.get('pop_hetero_informative', 0))
    res["pop_homo_count"] = int(pop_counts.get('pop_homo_like', 0))
    target_df = report_df[report_df['pop_class'] == 'pop_hetero_informative']
    total_sites = len(target_df)

    res["total_sites"] = total_sites
    if total_sites > 0:
        res["ref_sum"] = int(target_df['ref_depth'].sum())
        res["alt_sum"] = int(target_df['alt_depth'].sum())
        res["other_sum"] = int(target_df['other_depth'].sum())
        res["total_depth"] = int(target_df['total_depth'].sum())

        # Mutation Signature Calculation
        complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        mutation_types = []
        for r_row in target_df.itertuples():
            ref = r_row.ref.upper()
            alt = r_row.alt.upper()
            if ref in ['A', 'C', 'G', 'T'] and alt in ['A', 'C', 'G', 'T'] and ref != alt:
                if ref in ['A', 'G']:
                    std_ref, std_alt = complement_map[ref], complement_map[alt]
                else:
                    std_ref, std_alt = ref, alt
                mutation_types.append(f"{std_ref}>{std_alt}")

        if mutation_types:
            mut_counts = pd.Series(mutation_types).value_counts()
            mut_dist_df = pd.DataFrame(mut_counts).reset_index()
            mut_dist_df.columns = ['mutation_type', 'count']
            mut_dist_df['bin_id'] = bin_id
            mut_dist_df['percentage'] = mut_dist_df['count'] / len(mutation_types)
            for s_type in ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']:
                if s_type not in mut_dist_df['mutation_type'].values:
                    mut_dist_df = pd.concat([mut_dist_df, pd.DataFrame([{'mutation_type': s_type, 'count': 0, 'bin_id': bin_id, 'percentage': 0.0}])], ignore_index=True)
            mut_dist_df = mut_dist_df.sort_values('mutation_type').reset_index(drop=True)

        # Noise Rates
        res["on_target_noise_rate"] = float(res["other_sum"] / res["total_depth"]) if res["total_depth"] > 0 else 0.0
        off_target_errors = max(0, total_nm_errors - res["alt_sum"] - res["other_sum"])
        off_target_bases = max(1, total_mapped_bases - res["total_depth"])
        res["off_target_noise_rate"] = float(off_target_errors / off_target_bases)

        res["bin_BAF"] = float(res["alt_sum"] / res["total_depth"]) if res["total_depth"] > 0 else 0.0

        baf_counts = target_df['baf_class'].value_counts()
        res["hetero_like_count"] = int(baf_counts.get('hetero_like', 0))
        res["imbalance_count"] = int(baf_counts.get('imbalance', 0))
        res["homo_like_count"] = int(baf_counts.get('homo_like', 0))

        res["hetero_like_rate"] = float(res["hetero_like_count"] / total_sites)
        res["imbalance_rate"] = float(res["imbalance_count"] / total_sites)
        res["homo_like_rate"] = float(res["homo_like_count"] / total_sites)
        res["MAD_BAF"] = float(np.median(np.abs(target_df['BAF'] - 0.5)))

    ##############################################

    qc_mask = report_df['total_fragments'] >= 2
    qc_pass_fragments = int(report_df[qc_mask]["total_fragments"].sum()) if not report_df[qc_mask].empty else 0

    res.update({
        "total_fragments_sum": int(report_df["total_fragments"].sum()),
        "qc_pass_fragments": qc_pass_fragments,
        "total_trans_fragments": int(fragment_df['is_trans'].sum()),
        "total_cis_fragments": int(fragment_df['is_cis_alt'].sum())
    })

    denom = qc_pass_fragments if qc_pass_fragments > 0 else 1
    res["raw_bin_TER"] = float(res["total_trans_fragments"] / denom)
    res["raw_bin_CER"] = float(res["total_cis_fragments"] / denom)

    # Hetero/Homo 발생 확률 기반 Adjusted 에러율 보정
    res["adj_bin_TER"] = float(res["raw_bin_TER"] / (res["hetero_like_rate"] + 1e-4))
    res["adj_bin_CER"] = float(res["raw_bin_CER"] / (res["hetero_like_rate"] + res["homo_like_rate"] + 1e-4))

    summary_df = pd.DataFrame([res])
    for col, dtype in BIN_SUMMARY_SCHEMA.items():
        if col in summary_df.columns:
            summary_df[col] = summary_df[col].astype(dtype)

    return report_df, summary_df, mut_dist_df


def apply_low_quality_filter(df, min_depth=None, min_coverage=None):
    """
    Parquet 워커에서 병합되어 나온 최종 DataFrame(df)을 필터링합니다.
    성염색체는 보호하고 상염색체만 필터링합니다.
    min_depth/min_coverage 미지정 시 CFG["min_depth"]/CFG["min_coverage"] 사용.
    """
    if df is None or df.empty: return pd.DataFrame()

    if min_depth is None:
        min_depth = CFG["min_depth"]
    if min_coverage is None:
        min_coverage = CFG["min_coverage"]

    initial_count = len(df)
    sex_chrom_mask = df["chrom"].isin(["chrX", "chrY"])
    df_autosomes = df[~sex_chrom_mask].copy()
    df_sex_chroms = df[sex_chrom_mask].copy()

    # Depth 및 Breadth 필터링
    df_auto_filtered = df_autosomes[
        (df_autosomes["raw_count"] >= min_depth) &
        (df_autosomes["breadth_ratio"] >= min_coverage)
    ].copy()

    df_filtered = pd.concat([df_auto_filtered, df_sex_chroms], ignore_index=True)

    def _chrom_key(c):
        c_str = str(c).replace('chr', '')
        if c_str == 'X': return 23
        if c_str == 'Y': return 24
        try: return int(c_str)
        except: return 99

    df_filtered['chrom_sort_key'] = df_filtered['chrom'].apply(_chrom_key)
    df_filtered = df_filtered.sort_values(by=['chrom_sort_key', 'start']).drop(columns=['chrom_sort_key']).reset_index(drop=True)

    log(f"Quality Filter: Removed {initial_count - len(df_filtered)} bins.")
    return df_filtered

def fetch_sites_with_handle(vcf_handle, chrom, start, end):
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

    # 데이터가 없어도 무조건 컬럼(Header) 구조는 만들어서 반환하도록 강제
    df = pd.DataFrame(sites)
    if df.empty:
        return pd.DataFrame(columns=["chrom", "pos", "ref", "alt", "pop_af"])
    return df

def _parallel_chrom_worker(chrom_df, bam_path, vcf_file, min_mapq, thresholds, tmp_dir):
    chrom = chrom_df['chrom'].iloc[0]
    if not os.path.exists(vcf_file): return None, None, None

    tmp_pos_path = os.path.join(tmp_dir, f"tmp_pos_{chrom}.parquet")
    tmp_sum_path = os.path.join(tmp_dir, f"tmp_sum_{chrom}.parquet")
    tmp_mut_path = os.path.join(tmp_dir, f"tmp_mut_{chrom}.parquet")

    local_position_reports = []
    local_bin_summaries = []
    local_mut_distributions = []

    # bin extraction 관련 QC threshold: 명시적으로 안 넘어오면 config 기본값 사용
    min_baseq = thresholds.get("min_baseq", _BIN_EXTRACTION_CFG["min_baseq"])

    with pysam.AlignmentFile(bam_path, "rb", threads=1) as bam_handle, \
         pysam.VariantFile(vcf_file) as vcf_handle:

        for r in chrom_df.itertuples():
            site_df = fetch_sites_with_handle(vcf_handle, r.chrom, r.start, r.end)

            report_df, summary_df, mut_dist_df = process_bin_with_handles(
                bam_handle,
                {"chrom": r.chrom, "start": r.start, "end": r.end},
                site_df,
                min_mapq=min_mapq,
                min_baseq=min_baseq,
                thresholds=thresholds
            )

            if not report_df.empty: local_position_reports.append(report_df)
            if not mut_dist_df.empty: local_mut_distributions.append(mut_dist_df)

            if not summary_df.empty:
                summary_dict = summary_df.iloc[0].to_dict()
                summary_dict.update({"chrom": r.chrom, "start": r.start, "end": r.end})
                local_bin_summaries.append(summary_dict)

    if local_position_reports:
        pd.concat(local_position_reports, ignore_index=True).to_parquet(tmp_pos_path, engine='pyarrow', index=False)
    if local_bin_summaries:
        pd.DataFrame(local_bin_summaries).to_parquet(tmp_sum_path, engine='pyarrow', index=False)
    if local_mut_distributions:
        pd.concat(local_mut_distributions, ignore_index=True).to_parquet(tmp_mut_path, engine='pyarrow', index=False)

    del local_position_reports, local_bin_summaries, local_mut_distributions

    return (
        tmp_pos_path if os.path.exists(tmp_pos_path) else None,
        tmp_sum_path if os.path.exists(tmp_sum_path) else None,
        tmp_mut_path if os.path.exists(tmp_mut_path) else None
    )


def gc_correct_lowess_robust(
    df,
    count_col="raw_count",
    gc_col="gc",
    frac=0.5,
    pseudocount=0.25,
    min_count_for_fit=1,
    gc_range=(0.25, 0.75),
    mappability_col=None,
    blacklist_col=None,
):
    """
    Low-pass/WGA용 robust GC correction.

    output:
      log2_corrected_count : GC 보정된 log2 normalized count
      gc_fit               : GC별 fitted value
      gc_valid_for_fit     : LOWESS fit에 사용된 bin 여부
    """

    import numpy as np
    from statsmodels.nonparametric.smoothers_lowess import lowess

    df = df.copy()

    counts = df[count_col].astype(float).values
    gc = df[gc_col].astype(float).values

    # ----------------------------------
    # 1. library size normalization
    # ----------------------------------
    total = np.nansum(counts)

    if total <= 0:
        df["log2_corrected"] = np.nan
        df["gc_fit"] = np.nan
        df["gc_valid_for_fit"] = False
        return df, None

    # CPM scale
    norm_count = counts / total * 1_000_000

    y = np.log2(norm_count + pseudocount)

    # ----------------------------------
    # 2. valid bin for GC fit
    # ----------------------------------
    valid = (
        np.isfinite(gc) &
        np.isfinite(y) &
        (gc >= gc_range[0]) &
        (gc <= gc_range[1]) &
        (counts >= min_count_for_fit)
    )

    if mappability_col is not None and mappability_col in df.columns:
        valid &= df[mappability_col].fillna(0).values >= 0.75

    if blacklist_col is not None and blacklist_col in df.columns:
        valid &= ~df[blacklist_col].fillna(False).values

    # fit할 bin이 너무 적으면 correction 생략
    if valid.sum() < 100:
        df["log2_corrected"] = y
        df["gc_fit"] = np.nan
        df["gc_valid_for_fit"] = valid
        return df, (gc[valid], y[valid], None)

    # ----------------------------------
    # 3. LOWESS fit
    # ----------------------------------
    fit_valid = lowess(
        y[valid],
        gc[valid],
        frac=frac,
        return_sorted=False
    )

    # ----------------------------------
    # 4. 모든 bin에 대해 GC fit 예측 (LOWESS는 valid 위치만 반환하므로 interpolation)
    # ----------------------------------
    order = np.argsort(gc[valid])
    gc_sorted = gc[valid][order]
    fit_sorted = fit_valid[order]

    fit_all = np.interp(
        gc,
        gc_sorted,
        fit_sorted,
        left=fit_sorted[0],
        right=fit_sorted[-1]
    )

    # ----------------------------------
    # 5. GC corrected log2 normalized count
    # ----------------------------------
    baseline = np.nanmedian(fit_valid)

    df["log2_corrected"] = y - fit_all + baseline
    df["gc_fit"] = fit_all
    df["gc_valid_for_fit"] = valid
    correction = fit_all - baseline
    correction = np.clip(correction, -0.35, 0.35)

    df["log2_corrected"] = y - correction
    df["gc_correction"] = correction
    
    return df, (gc[valid], y[valid], fit_valid)
