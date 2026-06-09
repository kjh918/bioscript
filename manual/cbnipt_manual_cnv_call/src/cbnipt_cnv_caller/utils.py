import os
import pysam
import pandas as pd
import numpy as np

def log(msg): 
    print(f"[INFO] {msg}", flush=True)

def ensure_dir(path): 
    os.makedirs(path, exist_ok=True)


BIN_SUMMARY_SCHEMA = {
    "bin_id": "object",
    "raw_ref_depth_sum": "int32", "raw_alt_depth_sum": "int32", "raw_total_depth": "int32",
    "raw_total_fragments_sum": "int32", "qc_pass_fragments": "int32",
    
    "hetero_baf_mean": "float32", "hetero_baf_median": "float32", "hetero_sites_count": "int32",
    "homo_baf_mean": "float32","homo_baf_median": "float32","homo_sites_count": "int32",
    "other_low_baf_median": "float32","other_low_baf_mean": "float32",   "other_low_sites_count": "int32",
    "other_low_baf_median": "float32","other_high_baf_mean": "float32",  "other_high_sites_count": "int32",
    
    "total_trans_fragments": "int32","total_cis_fragments": "int32",
    "raw_bin_TER": "float32",
    "raw_bin_CER": "float32"
}
def summarize_and_classify_bin(fragment_df, site_df, bin_id, thresholds):
    """
    thresholds 예시: {'homo_max': 0.1, 'hetero_min': 0.4, 'hetero_max': 0.6, 'homo_min': 0.9}
    """
    if fragment_df.empty: 
        # 빈 결과일 경우 스키마에 맞춰 빈 DataFrame 반환
        empty_res = {k: (0 if 'int' in v else 0.0) for k, v in BIN_SUMMARY_SCHEMA.items()}
        empty_res["bin_id"] = bin_id
        return pd.DataFrame(), pd.DataFrame([empty_res])
    
    # 1. 포지션별 통계 산출
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
    report_df["observed_baf"] = (report_df["alt_depth"] / (report_df["ref_depth"] + report_df["alt_depth"])).fillna(0)
    
    # 2. 그룹화 (Logic: Others 분리)
    conditions = [
        (report_df['pop_af'] <= thresholds['homo_max']) | (report_df['pop_af'] >= thresholds['homo_min']),
        (report_df['pop_af'] >= thresholds['hetero_min']) & (report_df['pop_af'] <= thresholds['hetero_max']),
        (report_df['pop_af'] > thresholds['homo_max']) & (report_df['pop_af'] < thresholds['hetero_min']),
        (report_df['pop_af'] > thresholds['hetero_max']) & (report_df['pop_af'] < thresholds['homo_min'])
    ]
    choices = ['Homo', 'Hetero', 'Other_Low', 'Other_High']
    report_df['group'] = np.select(conditions, choices, default='Other')
    
    # 메모리 최적화
    report_df['group'] = report_df['group'].astype('category')
    
    # 3. QC 필터링 및 요약
    qc_group = report_df[report_df['total_fragments'] >= 2]
    qc_pass_fragments = int(qc_group["total_fragments"].sum()) if not qc_group.empty else 0
    
    # 그룹별 Aggregation
    summary_grouped = qc_group.groupby('group', observed=False)['observed_baf'].agg(['mean', 'median', 'count'])
    
    # 4. 결과 매핑
    res = {k: (0 if 'int' in v else 0.0) for k, v in BIN_SUMMARY_SCHEMA.items()}
    res["bin_id"] = bin_id
    
    # 그룹별 매핑
    group_map = {
        'Hetero': ('hetero_baf_mean', 'hetero_baf_median', 'hetero_sites_count'),
        'Homo': ('homo_baf_mean', 'homo_baf_median', 'homo_sites_count'),
        'Other_Low': ('other_low_baf_mean', 'other_low_baf_median', 'other_low_sites_count'),
        'Other_High': ('other_high_baf_mean', 'other_high_baf_median', 'other_high_sites_count')
    }
    
    for grp_name, (m, med, cnt) in group_map.items():
        if grp_name in summary_grouped.index:
            res[m] = float(summary_grouped.loc[grp_name, 'mean'])
            res[med] = float(summary_grouped.loc[grp_name, 'median'])
            res[cnt] = int(summary_grouped.loc[grp_name, 'count'])
            
    # 나머지 기본 지표 업데이트
    res.update({
        "raw_ref_depth_sum": int(report_df["ref_depth"].sum()),
        "raw_alt_depth_sum": int(report_df["alt_depth"].sum()),
        "raw_total_depth": int(report_df["ref_depth"].sum() + report_df["alt_depth"].sum()),
        "raw_total_fragments_sum": int(report_df["total_fragments"].sum()),
        "qc_pass_fragments": qc_pass_fragments,
        "total_trans_fragments": int(fragment_df['is_trans'].sum()),
        "total_cis_fragments": int(fragment_df['is_cis_alt'].sum())
    })
    
    denom = qc_pass_fragments if qc_pass_fragments > 0 else 1
    res["raw_bin_TER"] = float(res["total_trans_fragments"] / denom)
    res["raw_bin_CER"] = float(res["total_cis_fragments"] / denom)
    
    summary_df = pd.DataFrame([res])
    
    # 스키마 타입 적용 (데이터 타입 강제)
    for col, dtype in BIN_SUMMARY_SCHEMA.items():
        if col in summary_df.columns:
            summary_df[col] = summary_df[col].astype(dtype)
            
    return report_df, summary_df



def process_bin_with_handles(bam_handle, bin_info, site_df, min_mapq, min_baseq, hetero_range):
    chrom, start, end = bin_info['chrom'], bin_info['start'], bin_info['end']
    bin_sites = site_df[(site_df['chrom'] == chrom) & (site_df['pos'] > start) & (site_df['pos'] < end)]
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
    return summarize_and_classify_bin(pd.DataFrame(frag_list), bin_sites, f"{chrom}:{start}-{end}", hetero_range)


# -------------------------------------------------------------------------
# [최적화 2] 독립 병렬 워커: TSV 대신 Parquet 파일로 고속 Disk-Spill 저장
# -------------------------------------------------------------------------