import pysam
import numpy as np
import pandas as pd
import os
from utils import log

# [기존 요약 및 분석 함수들은 변동 없음]
def summarize_and_classify(fragment_df, site_df, bin_id, hetero_range=(0.4, 0.6)):
    if fragment_df.empty: return pd.DataFrame(), pd.DataFrame()
    pos_stats = {
        row.pos: {
            "bin_id": bin_id, "chrom": row.chrom, "pos": row.pos,
            "ref": row.ref, "alt": row.alt, "pop_af": row.pop_af,
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

    report_df["observed_baf"] = report_df["alt_depth"] / (report_df["ref_depth"] + report_df["alt_depth"])
    report_df["observed_baf"] = report_df["observed_baf"].fillna(0)
    report_df['TER'] = (report_df['trans_support'] / report_df['total_fragments']).fillna(0)
    report_df['CER'] = (report_df['cis_support'] / report_df['total_fragments']).fillna(0)

    h_min, h_max = hetero_range
    report_df['group'] = 'Other'
    report_df.loc[(report_df['pop_af'] >= h_min) & (report_df['pop_af'] <= h_max), 'group'] = 'Hetero'
    report_df.loc[(report_df['pop_af'] <= 0.1) | (report_df['pop_af'] >= 0.9), 'group'] = 'Homo'

    summary_grouped = report_df[report_df['total_fragments'] >= 2].groupby('group')['observed_baf'].agg(['mean', 'median', 'count'])
    
    res = {"bin_id": bin_id}
    for g in ['Hetero', 'Homo']:
        if g in summary_grouped.index:
            res[f"{g.lower()}_baf_mean"] = summary_grouped.loc[g, 'mean']
            res[f"{g.lower()}_baf_median"] = summary_grouped.loc[g, 'median']
            res[f"{g.lower()}_sites_count"] = summary_grouped.loc[g, 'count']
        else:
            res[f"{g.lower()}_baf_mean"] = res[f"{g.lower()}_baf_median"] = 0
            res[f"{g.lower()}_sites_count"] = 0
    
    res["total_trans_fragments"] = fragment_df['is_trans'].sum()
    res["total_cis_fragments"] = fragment_df['is_cis_alt'].sum()
    summary_df = pd.DataFrame([res])

    return report_df, summary_df

def process_bin_efficiently(bam_handle, bin_info, site_df, min_mapq, min_baseq, hetero_range=(0.4, 0.6)):
    chrom, start, end = bin_info['chrom'], bin_info['start'], bin_info['end']
    
    bin_sites = site_df[(site_df['chrom'] == chrom) & (site_df['pos'] >= start) & (site_df['pos'] < end)]
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

        if qname not in fragments:
            fragments[qname] = obs
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

    return summarize_and_classify(pd.DataFrame(frag_list), bin_sites, f"{chrom}:{start}-{end}", hetero_range)


# ---------------------------------------------------------
# [NEW] 마커 정보 프리페칭 (메모리 캐싱 기법)
# ---------------------------------------------------------
def preload_marker_info(vcf_path, bin_dataframe):
    """
    모든 Bin 영역에 대한 VCF 마커 정보를 사전에 로드하여 딕셔너리 캐시로 반환합니다.
    BAM 루프를 돌기 전 단 한 번만 실행됩니다.
    """
    log(f"Pre-loading VCF marker information from: {os.path.basename(vcf_path)}")
    marker_cache = {}
    
    with pysam.VariantFile(vcf_path) as vcf:
        for _, row in bin_dataframe.iterrows():
            chrom, start, end = row['chrom'], row['start'], row['end']
            bin_key = (chrom, start, end)
            sites = []
            
            try:
                for record in vcf.fetch(chrom, start, end):
                    kova = record.info.get("KOVA_AF")
                    kova_af = kova[0] if kova and isinstance(kova, (list, tuple)) else (kova if kova else 0.0)
                    
                    gnomad = record.info.get("GNOMAD_AF")
                    gnomad_af = gnomad[0] if gnomad and isinstance(gnomad, (list, tuple)) else (gnomad if gnomad else 0.0)

                    pop_af = kova_af if kova_af > 0 else gnomad_af
                    if pop_af > 0:
                        sites.append({
                            "chrom": record.chrom, "pos": record.pos,
                            "ref": record.ref, "alt": record.alts[0],
                            "pop_af": pop_af
                        })
            except Exception as e:
                pass
                
            marker_cache[bin_key] = pd.DataFrame(sites)
            
    log(f"Pre-loading complete. Cached {len(marker_cache)} bins.")
    return marker_cache


# ---------------------------------------------------------
# [MODIFIED] VCF 경로 대신 marker_cache를 주입받는 run 함수
# ---------------------------------------------------------
def run(bam_path, marker_cache, bin_dataframe, min_mapq=30, min_baseq=20, hetero_range=(0.4, 0.6)):
    """
    전체 파이프라인 실행 함수 (VCF I/O 제거 버전)
    """
    rawdata_results = []
    summary_results = []
    
    with pysam.AlignmentFile(bam_path, "rb", threads=4) as bam_handle:
        for _, row in bin_dataframe.iterrows():
            chrom, start, end = row['chrom'], row['start'], row['end']
            bin_key = (chrom, start, end)
            
            # [OPTIMIZED] VCF 파일 탐색 없이 메모리 캐시에서 즉시 추출 (O(1))
            site_df = marker_cache.get(bin_key, pd.DataFrame())
            if site_df.empty: continue
            
            report_df, summary_df = process_bin_efficiently(
                bam_handle, 
                {"chrom": chrom, "start": start, "end": end}, 
                site_df, min_mapq, min_baseq, hetero_range
            )
            if not report_df.empty and not summary_df.empty:
                rawdata_results.append(report_df)
                summary_results.append(summary_df)
    
    if not rawdata_results or not summary_results:
        return pd.DataFrame(), pd.DataFrame()
        
    return pd.concat(rawdata_results, ignore_index=True), pd.concat(summary_results, ignore_index=True)


if __name__ == "__main__":
    bam_path_lits = [
        "/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Resources/BAM/cbNIPT_24_04_13.analysisReady.bam", 
        "/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Resources/BAM/cbNIPT_24_04_14.analysisReady.bam", 
        "/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Resources/BAM/cbNIPT_24_04_15.analysisReady.bam", 
    ]

    chrom_list = ['chr13'] 

    for chrom in chrom_list:
        vcf_path = '/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Resources/reference/KOVA_v7/kova_sites_vcf/KOVA_v7_SNPOnly_QCSites_af.chr13.vcf.gz'
        
        # 1. 윈도우/빈 설정 데이터프레임 생성
        bin_dataframe = pd.DataFrame([{
            "chrom": chrom, "start": 28470000, "end": 28570000
        }])
        
        # ---------------------------------------------------------
        # [CRITICAL OPTIMIZATION] 마커 정보 프리페치 실행 (BAM 루프 바깥에서 딱 한 번만!)
        # ---------------------------------------------------------
        marker_cache = preload_marker_info(vcf_path, bin_dataframe)
        
        # 2. 캐시를 재사용하여 여러 BAM 파일 고속 처리
        for bam_path in bam_path_lits:
            log(f"Processing BAM: {os.path.basename(bam_path)}")
            
            # 주입되는 인자가 vcf_path에서 marker_cache로 변경됨
            results_df, summary_df = run(bam_path, marker_cache, bin_dataframe, min_mapq=30, min_baseq=20, hetero_range=(0.4, 0.6))
            print(results_df)
            print(f"\n=== Result for {os.path.basename(bam_path)} ===")
            print(summary_df.columns)

            exit()