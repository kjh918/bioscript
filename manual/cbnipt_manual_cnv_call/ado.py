import pysam
import numpy as np
import pandas as pd
import os
from utils import log

def summarize_and_classify(fragment_df, site_df, bin_id, hetero_range=(0.4, 0.6)):
    """
    포지션별 집계와 Hetero/Homo 분류 통계를 단일 패스로 처리합니다.
    """
    if fragment_df.empty:
        return pd.DataFrame(), pd.DataFrame()

    # 1. 사이트 데이터 준비 (O(1) 조회를 위한 사전)
    pos_stats = {
        row.pos: {
            "bin_id": bin_id, "chrom": row.chrom, "pos": row.pos,
            "ref": row.ref, "alt": row.alt, "pop_af": row.pop_af,
            "ref_depth": 0, "alt_depth": 0, "cis_support": 0, "trans_support": 0, "total_fragments": 0
        } for row in site_df.itertuples()
    }

    # 2. 프래그먼트 데이터 집계
    # itertuples() 대신 리스트화된 dict 접근이 속도면에서 유리할 수 있음
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

    # 3. 지표 계산 (Vectorized)
    report_df["observed_baf"] = report_df["alt_depth"] / (report_df["ref_depth"] + report_df["alt_depth"])
    report_df["observed_baf"] = report_df["observed_baf"].fillna(0)
    report_df['TER'] = (report_df['trans_support'] / report_df['total_fragments']).fillna(0)
    report_df['CER'] = (report_df['cis_support'] / report_df['total_fragments']).fillna(0)

    # 4. Hetero/Homo 그룹화 통계
    h_min, h_max = hetero_range
    report_df['group'] = 'Other'
    report_df.loc[(report_df['pop_af'] >= h_min) & (report_df['pop_af'] <= h_max), 'group'] = 'Hetero'
    report_df.loc[(report_df['pop_af'] <= 0.1) | (report_df['pop_af'] >= 0.9), 'group'] = 'Homo'

    # Groupby를 이용한 요약 통계 계산
    summary_grouped = report_df[report_df['total_fragments'] >= 2].groupby('group')['observed_baf'].agg(['mean', 'median', 'count'])
    
    # 결과 요약 생성
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
    """
    열려있는 BAM 핸들을 사용하여 특정 Bin의 Haplotype을 분석합니다.
    """
    chrom, start, end = bin_info['chrom'], bin_info['start'], bin_info['end']
    
    # 1. 사이트 룩업 생성
    bin_sites = site_df[(site_df['chrom'] == chrom) & (site_df['pos'] >= start) & (site_df['pos'] < end)]
    if bin_sites.empty: return pd.DataFrame(), pd.DataFrame()
    site_lookup = {row.pos: {"ref": row.ref, "alt": row.alt} for row in bin_sites.itertuples()}

    # 2. 프래그먼트 수집 (메모리 절약형: Read 객체 전체를 저장하지 않음)
    fragments = {}
    for read in bam_handle.fetch(chrom, start, end):
        if read.is_unmapped or read.is_duplicate or read.is_secondary or read.is_supplementary: continue
        if read.mapping_quality < min_mapq: continue
            
        qname = read.query_name
        # 필요한 정보만 미리 파싱하여 저장
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
            # R1/R2 정보 통합 (더 높은 퀄리티 우선)
            for g_pos, data in obs.items():
                if g_pos not in fragments[qname] or data['qual'] > fragments[qname][g_pos]['qual']:
                    fragments[qname][g_pos] = data

    # 3. 프래그먼트 데이터 구조화
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

def fetch_sites_from_vcf_gz(vcf_path, chrom, start, end):
    """
    vcf.gz에서 INFO 필드의 AF를 포함하여 사이트를 추출합니다.
    """
    vcf = pysam.VariantFile(vcf_path)
    sites = []
    try:
        for record in vcf.fetch(chrom, start, end):
            # KOVA_AF가 튜플 형태일 경우 첫 번째 값 추출
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
    except: pass
    vcf.close()
    return pd.DataFrame(sites)

def run(bam_path, vcf_path, bin_dataframe, min_mapq=30, min_baseq=20, hetero_range=(0.4, 0.6)):
    """
    전체 파이프라인 실행 함수
    """
    rawdata_results = []
    summary_results = []
    with pysam.AlignmentFile(bam_path, "rb", threads=4) as bam_handle:
        for _, row in bin_dataframe.iterrows():
            chrom, start, end = row['chrom'], row['start'], row['end']
            site_df = fetch_sites_from_vcf_gz(vcf_path, chrom, start, end)
            if site_df.empty: continue
            
            report_df, summary_df = process_bin_efficiently(bam_handle, {"chrom": chrom, "start": start, "end": end}, site_df, min_mapq, min_baseq, hetero_range)
            if not report_df.empty and not summary_df.empty:
                rawdata_results.append(report_df)
                summary_results.append(summary_df)
    
    return pd.concat(rawdata_results), pd.concat(summary_results)




if __name__ == "__main__":
    # 테스트용 실행 코드 (실제 사용 시 cli.py에서 호출)
    bam_path_lits = [
    #    "/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Resources/BAM/cbNIPT_24_04_01_DS.analysisReady.bam",# WGA
    #    "/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Resources/BAM/cbNIPT_24_04_02_DS.analysisReady.bam", # WGA
    #    "/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Resources/BAM/cbNIPT_24_04_03_DS.analysisReady.bam", # WGA
        #"/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Resources/BAM/cbNIPT_24_04_04.analysisReady.bam", # WGA
        "/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Resources/BAM/cbNIPT_24_04_13.analysisReady.bam", # WGS
        "/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Resources/BAM/cbNIPT_24_04_14.analysisReady.bam", # WGS
        "/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Resources/BAM/cbNIPT_24_04_15.analysisReady.bam", # WGS

        #"/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Resources/BAM/cbNIPT_26_02_01_HTR8-Svneo.analysisReady.bam", # WGS
        #'/storage/home/jhkim/Projects/cbNIPT/260330-GCX-cbNIPT-InsilicoMosaicism/Results/insilico/cbNIPT_24_04_03_DS_Ratio_1-1-Insilico_cbNIPT_24_04_06.bam',
        #"/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Resources/BAM/cbNIPT_24_04_13.analysisReady.bam"
    ]

    chrom_list = ['chr13'] 

    for chrom in chrom_list:
        
        #vcf_path = f"/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Resources/reference/KOVA_v7/kova_sites_vcf/KOVA_v7_SNPOnly_QCSites_af.{chrom}.vcf.gz"
        vcf_path = '/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Resources/reference/KOVA_v7/kova_sites_vcf/KOVA_v7_SNPOnly_QCSites_af.chr13.vcf.gz'
        bin_info = {"chrom": chrom, "start": 28470000, "end": 28570000}
        site_df = fetch_sites_from_vcf_gz(vcf_path, bin_info["chrom"], bin_info["start"], bin_info["end"])

        for bam_path in bam_path_lits:
            log(f"Processing BAM: {os.path.basename(bam_path)}")
            
            with pysam.AlignmentFile(bam_path, "rb", threads=4) as bam_handle:
                for chrom in chrom_list:

                    site_df = fetch_sites_from_vcf_gz(vcf_path, chrom, 28470000, 28570000)
                    
                    bin_info = {"chrom": chrom, "start": 28470000, "end": 28570000}
                    pos_report, summary = process_bin_efficiently(bam_handle, bin_info, site_df, 30, 20)
                    
                    pos_report.to_csv('temp.csv',sep='\t')
                    if not summary.empty:
                        print(summary)
                        