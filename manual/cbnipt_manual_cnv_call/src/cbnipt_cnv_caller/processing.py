import numpy as np
import pandas as pd
import pysam
import os 
from utils import log
from statsmodels.nonparametric.smoothers_lowess import lowess

# -------------------------------------------------------------------------
# [데이터 스키마 정의] Breadth와 Midpoint Read Count 추가
# -------------------------------------------------------------------------
BIN_SUMMARY_SCHEMA = {
    "bin_id": "object",
    "raw_count": "int32",               "breadth_ratio": "float32",

    "total_sites": "int32",             # Population 기준 Informative SNP 총 개수
    "ref_sum": "int32",                 "alt_sum": "int32",
    "total_depth": "int32",             "bin_BAF": "float32",
    
    "pop_hetero_count": "int32",        # AF 0.4~0.6 인 마커 수 (total_sites와 동일)
    "pop_homo_count": "int32",          # AF <=0.1 or >=0.9 인 마커 수 (Reference 유사/고정 영역)
    
    # [NEW] 관측된 절대 개수 (신뢰도 가중치용)
    "hetero_like_count": "int32",       "imbalance_count": "int32",
    "homo_like_count": "int32",
    
    # 관측된 비율
    "hetero_like_rate": "float32",      "imbalance_rate": "float32",
    "homo_like_rate": "float32",        "MAD_BAF": "float32",
    
    "total_fragments_sum": "int32",     "qc_pass_fragments": "int32",
    "total_trans_fragments": "int32",   "total_cis_fragments": "int32",
    "raw_bin_TER": "float32",           "raw_bin_CER": "float32"
}

def process_bin_with_handles(bam_handle, bin_info, site_df, min_mapq, min_baseq, thresholds):
    """
    One-Pass BAM 처리
    Midpoint 카운팅과 Breadth of Coverage 계산, VCF 파싱을 단 하나의 for 루프에서 동시에 해결합니다.
    """
    chrom, start, end = bin_info['chrom'], bin_info['start'], bin_info['end']
    bin_len = end - start
    
    # [CRITICAL FIX] bin_sites가 무조건 존재하도록 강제 할당
    if site_df.empty or 'chrom' not in site_df.columns:
        # 변이가 없거나 컬럼이 깨져서 들어와도 빈 DataFrame으로 안전하게 초기화
        bin_sites = pd.DataFrame(columns=["chrom", "pos", "ref", "alt", "pop_af"])
    else:
        # 정상적인 경우 필터링 진행
        bin_sites = site_df[(site_df['chrom'] == chrom) & (site_df['pos'] > start) & (site_df['pos'] < end)]
        
    site_lookup = {row.pos: {"ref": row.ref, "alt": row.alt} for row in bin_sites.itertuples()} if not bin_sites.empty else {}

    fragments = {}
    raw_count = 0  
    coverage_mask = np.zeros(bin_len, dtype=bool)

    for read in bam_handle.fetch(chrom, start, end):
        if read.is_unmapped or read.is_duplicate or read.is_secondary or read.is_supplementary: continue
        if read.mapping_quality < min_mapq: continue
        
        # 1. Breadth 계산
        ref_start = max(read.reference_start, start)
        ref_end = min(read.reference_end, end)
        if ref_start < ref_end:
            coverage_mask[ref_start - start : ref_end - start] = True
            
        # 2. Midpoint 카운팅
        mid = (read.reference_start + read.reference_end) // 2
        if start <= mid < end:
            raw_count += 1

        # 3. VCF Phasing 및 BAF 추출 (변이가 없으면 이 아래는 자동 스킵)
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
            
    return summarize_and_classify_bin(pd.DataFrame(frag_list), bin_sites, f"{chrom}:{start}-{end}", raw_count, breadth_ratio, thresholds)
 

def summarize_and_classify_bin(fragment_df, site_df, bin_id, raw_count, breadth_ratio, thresholds):
    res = {k: (0 if 'int' in v else 0.0) for k, v in BIN_SUMMARY_SCHEMA.items()}
    res["bin_id"] = bin_id
    res["raw_count"] = raw_count
    res["breadth_ratio"] = breadth_ratio
    
    if fragment_df.empty or site_df.empty: 
        return pd.DataFrame(), pd.DataFrame([res])
    
    # 1. 포지션별 기본 통계 적재
    pos_stats = {
        row.pos: {
            "bin_id": bin_id, "chrom": row.chrom, "pos": row.pos, "ref": row.ref, "alt": row.alt, "pop_af": row.pop_af,
            "ref_depth": 0, "alt_depth": 0, "cis_support": 0, "trans_support": 0, "total_fragments": 0
        } for row in site_df.itertuples()
    }
    
    # Fragment 정보 파싱
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
    if report_df.empty: return pd.DataFrame(), pd.DataFrame([res])

    # Site-level 계산
    report_df["total_depth"] = report_df["ref_depth"] + report_df["alt_depth"]
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
    
    # Bin-level Summary (pop_hetero_informative 기준 추출)
    pop_counts = report_df['pop_class'].value_counts()
    res["pop_hetero_count"] = int(pop_counts.get('pop_hetero_informative', 0))
    res["pop_homo_count"] = int(pop_counts.get('pop_homo_like', 0))
    target_df = report_df[report_df['pop_class'] == 'pop_hetero_informative']
    total_sites = len(target_df)
    
    res["total_sites"] = total_sites
    if total_sites > 0:
        res["ref_sum"] = int(target_df['ref_depth'].sum())
        res["alt_sum"] = int(target_df['alt_depth'].sum())
        res["total_depth"] = int(target_df['total_depth'].sum())
        res["bin_BAF"] = float(res["alt_sum"] / res["total_depth"]) if res["total_depth"] > 0 else 0.0
        
        # [NEW] 절대 카운트(Count) 먼저 추출
        baf_counts = target_df['baf_class'].value_counts()
        res["hetero_like_count"] = int(baf_counts.get('hetero_like', 0))
        res["imbalance_count"] = int(baf_counts.get('imbalance', 0))
        res["homo_like_count"] = int(baf_counts.get('homo_like', 0))
        
        # [NEW] 추출된 카운트를 바탕으로 비율(Rate) 계산
        res["hetero_like_rate"] = float(res["hetero_like_count"] / total_sites)
        res["imbalance_rate"] = float(res["imbalance_count"] / total_sites)
        res["homo_like_rate"] = float(res["homo_like_count"] / total_sites)
        
        res["MAD_BAF"] = float(np.median(np.abs(target_df['BAF'] - 0.5)))
        
    # Fragment Phasing Summary
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

    summary_df = pd.DataFrame([res])
    for col, dtype in BIN_SUMMARY_SCHEMA.items():
        if col in summary_df.columns:
            summary_df[col] = summary_df[col].astype(dtype)
            
    return report_df, summary_df


# [수정] 인자에 thresholds 추가
# -------------------------------------------------------------------------
# 기존 필터링 및 GC 보정 함수 (그대로 사용 가능합니다)
# -------------------------------------------------------------------------
def apply_low_quality_filter(df, min_depth=1, min_coverage=0.5):
    """
    Parquet 워커에서 병합되어 나온 최종 DataFrame(df)을 필터링합니다.
    성염색체는 보호하고 상염색체만 필터링합니다.
    """
    if df is None or df.empty: return pd.DataFrame()
        
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
    
    def chrom_key(c):
        c_str = str(c).replace('chr', '')
        if c_str == 'X': return 23
        if c_str == 'Y': return 24
        try: return int(c_str)
        except: return 99

    df_filtered['chrom_sort_key'] = df_filtered['chrom'].apply(chrom_key)
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
    
    # [수정] 데이터가 없어도 무조건 컬럼(Header) 구조는 만들어서 반환하도록 강제
    df = pd.DataFrame(sites)
    if df.empty:
        return pd.DataFrame(columns=["chrom", "pos", "ref", "alt", "pop_af"])
    return df

def _parallel_chrom_worker(chrom_df, bam_path, vcf_file, min_mapq, thresholds, tmp_dir):
    """
    [수정 사항]
    1. 파라미터로 thresholds 딕셔너리를 받습니다.
    2. site_df가 비어있어도 BAM Coverage를 구하기 위해 continue 하지 않고 무조건 함수를 호출합니다.
    """
    chrom = chrom_df['chrom'].iloc[0]
    
    if not os.path.exists(vcf_file):
        return None, None
        
    tmp_pos_path = os.path.join(tmp_dir, f"tmp_pos_{chrom}.parquet")
    tmp_sum_path = os.path.join(tmp_dir, f"tmp_sum_{chrom}.parquet")
    
    local_position_reports = []
    local_bin_summaries = []
    
    with pysam.AlignmentFile(bam_path, "rb", threads=1) as bam_handle, \
         pysam.VariantFile(vcf_file) as vcf_handle:
         
        for r in chrom_df.itertuples():
            site_df = fetch_sites_with_handle(vcf_handle, r.chrom, r.start, r.end)

            report_df, summary_df = process_bin_with_handles(
                bam_handle, 
                {"chrom": r.chrom, "start": r.start, "end": r.end}, 
                site_df, 
                min_mapq=min_mapq, 
                min_baseq=20, 
                thresholds=thresholds  # <- hetero_range 대신 thresholds 딕셔너리 전달
            )
            
            # 변이가 있는 경우에만 Position Report 적재
            if not report_df.empty:
                local_position_reports.append(report_df)
                
            # Summary는 Coverage 정보를 항상 담고 있으므로 무조건 적재
            if not summary_df.empty:
                summary_dict = summary_df.iloc[0].to_dict()
                summary_dict.update({"chrom": r.chrom, "start": r.start, "end": r.end})
                local_bin_summaries.append(summary_dict)

    # 파이썬 리스트 메모리 해제 및 Parquet 엔진으로 직렬화 저장
    if local_position_reports:
        pd.concat(local_position_reports, ignore_index=True).to_parquet(tmp_pos_path, engine='pyarrow', index=False)
        
    if local_bin_summaries:
        pd.DataFrame(local_bin_summaries).to_parquet(tmp_sum_path, engine='pyarrow', index=False)

    del local_position_reports
    del local_bin_summaries

    return (tmp_pos_path if os.path.exists(tmp_pos_path) else None, 
            tmp_sum_path if os.path.exists(tmp_sum_path) else None)

def gc_correct_lowess(df, frac=0.2):
    x, y = df["gc"].values, np.log2(df["raw_count"].values + 1)
    valid = np.isfinite(x) & (y > 0)
    fit = lowess(y[valid], x[valid], frac=frac, return_sorted=False)
    df["log2_corrected"] = np.nan
    df.loc[valid, "log2_corrected"] = y[valid] - fit + np.nanmedian(fit)
    return df, (x[valid], y[valid], fit)