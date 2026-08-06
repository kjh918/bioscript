import numpy as np
import pandas as pd
import pysam
import os
from utils import log
from statsmodels.nonparametric.smoothers_lowess import lowess

from utils import chrom_key
from rules import CFG

_BIN_EXTRACTION_CFG = CFG["BIN_EXTRACTION"]

# NIPT fragment size threshold (bp)
SHORT_FRAG_MAX = 150  # <= 150bp : short (fetal-enriched)
LONG_FRAG_MIN  = 151  # >= 151bp : long  (maternal-enriched)

# -------------------------------------------------------------------------
# [데이터 스키마 정의] short / long fragment 분리 지표
# -------------------------------------------------------------------------
def _make_frag_schema(prefix: str) -> dict:
    """short_ / long_ prefix 지표 블록을 생성합니다."""
    return {
        f"{prefix}raw_count":            "int32",
        f"{prefix}breadth_ratio":        "float32",

        f"{prefix}total_sites":          "int32",
        f"{prefix}ref_sum":              "int32",
        f"{prefix}alt_sum":              "int32",
        f"{prefix}other_sum":            "int32",
        f"{prefix}total_depth":          "int32",
        f"{prefix}bin_BAF":              "float32",

        f"{prefix}pop_hetero_count":     "int32",
        f"{prefix}pop_homo_count":       "int32",

        f"{prefix}hetero_like_count":    "int32",
        f"{prefix}imbalance_count":      "int32",
        f"{prefix}homo_like_count":      "int32",
        f"{prefix}hetero_like_rate":     "float32",
        f"{prefix}imbalance_rate":       "float32",
        f"{prefix}homo_like_rate":       "float32",
        f"{prefix}MAD_BAF":              "float32",

        f"{prefix}on_target_noise_rate": "float32",
        f"{prefix}off_target_noise_rate":"float32",

        f"{prefix}total_fragments_sum":  "int32",
        f"{prefix}qc_pass_fragments":    "int32",
        f"{prefix}total_trans_fragments":"int32",
        f"{prefix}total_cis_fragments":  "int32",

        f"{prefix}raw_bin_TER":          "float32",
        f"{prefix}raw_bin_CER":          "float32",
        f"{prefix}adj_bin_TER":          "float32",
        f"{prefix}adj_bin_CER":          "float32",
    }


BIN_SUMMARY_SCHEMA = {
    "bin_id": "object",
    **_make_frag_schema("short_"),
    **_make_frag_schema("long_"),
}


# -------------------------------------------------------------------------
# 내부 헬퍼: fragment_list → 단일 prefix 블록 summary dict
# -------------------------------------------------------------------------
def _summarize_fragment_group(
    fragment_df: pd.DataFrame,
    site_df: pd.DataFrame,
    bin_id: str,
    raw_count: int,
    breadth_ratio: float,
    thresholds: dict,
    total_mapped_bases: int,
    total_nm_errors: int,
    prefix: str,
) -> dict:
    """
    fragment_df / site_df 를 받아 {prefix}* 컬럼 dict를 반환합니다.
    summarize_and_classify_bin의 핵심 로직을 그대로 유지하되
    반환값을 prefix dict로 래핑합니다.
    """
    p = prefix  # 타이핑 단축

    res = {f"{p}{k}": (0 if "int" in v else 0.0)
           for k, v in _make_frag_schema(p).items()}
    # _make_frag_schema 키에는 이미 prefix가 붙어 있으므로 그대로 사용
    res[f"{p}raw_count"]    = raw_count
    res[f"{p}breadth_ratio"] = breadth_ratio

    if fragment_df.empty or site_df.empty:
        return res

    # 1. 포지션별 기본 통계
    pos_stats = {
        row.pos: {
            "bin_id": bin_id, "chrom": row.chrom, "pos": row.pos,
            "ref": row.ref, "alt": row.alt, "pop_af": row.pop_af,
            "ref_depth": 0, "alt_depth": 0, "other_depth": 0,
            "cis_support": 0, "trans_support": 0, "total_fragments": 0,
        }
        for row in site_df.itertuples()
    }

    for frag in fragment_df.to_dict("records"):
        for pos, data in frag["obs_dict"].items():
            if pos not in pos_stats:
                continue
            stat = pos_stats[pos]
            stat["total_fragments"] += 1
            if data["base"] == stat["alt"].upper():
                stat["alt_depth"] += 1
            elif data["base"] == stat["ref"].upper():
                stat["ref_depth"] += 1
            else:
                stat["other_depth"] += 1
            if frag["is_cis_alt"]:
                stat["cis_support"] += 1
            if frag["is_trans"]:
                stat["trans_support"] += 1

    report_df = pd.DataFrame(pos_stats.values())
    if report_df.empty:
        return res

    report_df["total_depth"] = (
        report_df["ref_depth"] + report_df["alt_depth"] + report_df["other_depth"]
    )
    report_df["BAF"] = (report_df["alt_depth"] / report_df["total_depth"]).fillna(0)

    cond_pop = [
        (report_df["pop_af"] >= thresholds["hetero_min"]) & (report_df["pop_af"] <= thresholds["hetero_max"]),
        (report_df["pop_af"] <= thresholds["homo_max"])   | (report_df["pop_af"] >= thresholds["homo_min"]),
    ]
    report_df["pop_class"] = np.select(
        cond_pop, ["pop_hetero_informative", "pop_homo_like"], default="intermediate"
    )

    cond_baf = [
        (report_df["BAF"] >= thresholds["hetero_min"]) & (report_df["BAF"] <= thresholds["hetero_max"]),
        (report_df["BAF"] <= thresholds["homo_max"])   | (report_df["BAF"] >= thresholds["homo_min"]),
    ]
    report_df["baf_class"] = np.select(
        cond_baf, ["hetero_like", "homo_like"], default="imbalance"
    )

    pop_counts  = report_df["pop_class"].value_counts()
    res[f"{p}pop_hetero_count"] = int(pop_counts.get("pop_hetero_informative", 0))
    res[f"{p}pop_homo_count"]   = int(pop_counts.get("pop_homo_like", 0))

    target_df   = report_df[report_df["pop_class"] == "pop_hetero_informative"]
    total_sites = len(target_df)
    res[f"{p}total_sites"] = total_sites

    if total_sites > 0:
        res[f"{p}ref_sum"]     = int(target_df["ref_depth"].sum())
        res[f"{p}alt_sum"]     = int(target_df["alt_depth"].sum())
        res[f"{p}other_sum"]   = int(target_df["other_depth"].sum())
        res[f"{p}total_depth"] = int(target_df["total_depth"].sum())

        td = res[f"{p}total_depth"]
        res[f"{p}on_target_noise_rate"] = float(res[f"{p}other_sum"] / td) if td > 0 else 0.0

        off_errors = max(0, total_nm_errors - res[f"{p}alt_sum"] - res[f"{p}other_sum"])
        off_bases  = max(1, total_mapped_bases - td)
        res[f"{p}off_target_noise_rate"] = float(off_errors / off_bases)

        res[f"{p}bin_BAF"] = float(res[f"{p}alt_sum"] / td) if td > 0 else 0.0

        baf_counts = target_df["baf_class"].value_counts()
        res[f"{p}hetero_like_count"] = int(baf_counts.get("hetero_like", 0))
        res[f"{p}imbalance_count"]   = int(baf_counts.get("imbalance",   0))
        res[f"{p}homo_like_count"]   = int(baf_counts.get("homo_like",   0))

        res[f"{p}hetero_like_rate"] = float(res[f"{p}hetero_like_count"] / total_sites)
        res[f"{p}imbalance_rate"]   = float(res[f"{p}imbalance_count"]   / total_sites)
        res[f"{p}homo_like_rate"]   = float(res[f"{p}homo_like_count"]   / total_sites)
        res[f"{p}MAD_BAF"]          = float(np.median(np.abs(target_df["BAF"] - 0.5)))

    qc_mask = report_df["total_fragments"] >= 2
    qc_pass = int(report_df[qc_mask]["total_fragments"].sum()) if not report_df[qc_mask].empty else 0

    res[f"{p}total_fragments_sum"]   = int(report_df["total_fragments"].sum())
    res[f"{p}qc_pass_fragments"]     = qc_pass
    res[f"{p}total_trans_fragments"] = int(fragment_df["is_trans"].sum())
    res[f"{p}total_cis_fragments"]   = int(fragment_df["is_cis_alt"].sum())

    denom = qc_pass if qc_pass > 0 else 1
    res[f"{p}raw_bin_TER"] = float(res[f"{p}total_trans_fragments"] / denom)
    res[f"{p}raw_bin_CER"] = float(res[f"{p}total_cis_fragments"]   / denom)

    hl = res[f"{p}hetero_like_rate"]
    ho = res[f"{p}homo_like_rate"]
    res[f"{p}adj_bin_TER"] = float(res[f"{p}raw_bin_TER"] / (hl + 1e-4))
    res[f"{p}adj_bin_CER"] = float(res[f"{p}raw_bin_CER"] / (hl + ho + 1e-4))

    return res


# -------------------------------------------------------------------------
# 메인 bin 처리 함수 (BAM fetch 1회, short/long 동시 분기)
# -------------------------------------------------------------------------
def process_bin_with_handles(bam_handle, bin_info, site_df, min_mapq, min_baseq, thresholds):
    chrom, start, end = bin_info["chrom"], bin_info["start"], bin_info["end"]
    bin_len = end - start

    if site_df.empty or "chrom" not in site_df.columns:
        bin_sites = pd.DataFrame(columns=["chrom", "pos", "ref", "alt", "pop_af"])
    else:
        bin_sites = site_df[
            (site_df["chrom"] == chrom) &
            (site_df["pos"] > start) &
            (site_df["pos"] < end)
        ]

    site_lookup = (
        {row.pos: {"ref": row.ref, "alt": row.alt} for row in bin_sites.itertuples()}
        if not bin_sites.empty else {}
    )

    # short / long 분리 카운터
    short_raw_count = 0
    long_raw_count  = 0

    short_coverage_mask = np.zeros(bin_len, dtype=bool)
    long_coverage_mask  = np.zeros(bin_len, dtype=bool)

    short_mapped_bases = 0;  short_nm_errors = 0
    long_mapped_bases  = 0;  long_nm_errors  = 0

    short_fragments: dict = {}
    long_fragments:  dict = {}

    # ── BAM fetch: 단 1회 ──────────────────────────────────────────────
    for read in bam_handle.fetch(chrom, start, end):
        if (read.is_unmapped or read.is_duplicate or
                read.is_secondary or read.is_supplementary):
            continue
        if read.mapping_quality < min_mapq:
            continue

        # fragment 크기: paired-end tlen, single-end은 query_length 사용
        if read.is_paired and read.template_length != 0:
            frag_size = abs(read.template_length)
        else:
            frag_size = read.query_length or 0

        is_short = (frag_size <= SHORT_FRAG_MAX)

        # coverage mask 업데이트
        ref_start = max(read.reference_start, start)
        ref_end   = min(read.reference_end,   end)
        if ref_start < ref_end:
            sl = slice(ref_start - start, ref_end - start)
            if is_short:
                short_coverage_mask[sl] = True
            else:
                long_coverage_mask[sl]  = True

        # midpoint-based raw count
        mid = (read.reference_start + read.reference_end) // 2
        if start <= mid < end:
            if is_short:
                short_raw_count += 1
            else:
                long_raw_count  += 1

        # NM / mapped bases
        aln_len = read.query_alignment_length
        nm      = read.get_tag("NM") if read.has_tag("NM") else 0
        if is_short:
            short_mapped_bases += aln_len;  short_nm_errors += nm
        else:
            long_mapped_bases  += aln_len;  long_nm_errors  += nm

        # ── SNP site 관찰 ────────────────────────────────────────────
        if not site_lookup:
            continue

        qname = read.query_name
        obs   = {}
        for q_pos, r_pos in read.get_aligned_pairs(matches_only=True):
            g_pos = r_pos + 1
            if g_pos in site_lookup:
                qual = read.query_qualities[q_pos]
                if qual < min_baseq:
                    continue
                obs[g_pos] = {"base": read.query_sequence[q_pos].upper(), "qual": qual}

        if not obs:
            continue

        target_frags = short_fragments if is_short else long_fragments
        if qname not in target_frags:
            target_frags[qname] = obs
        else:
            for g_pos, data in obs.items():
                if g_pos not in target_frags[qname] or data["qual"] > target_frags[qname][g_pos]["qual"]:
                    target_frags[qname][g_pos] = data

    # ── breadth ratio ────────────────────────────────────────────────
    short_breadth = float(np.sum(short_coverage_mask) / bin_len)
    long_breadth  = float(np.sum(long_coverage_mask)  / bin_len)

    # ── fragment DataFrame 변환 (공통 헬퍼) ─────────────────────────
    def _build_frag_df(frag_dict: dict) -> pd.DataFrame:
        frag_list = []
        for qname, obs in frag_dict.items():
            status_list = []
            for p, d in obs.items():
                site = site_lookup[p]
                if   d["base"] == site["alt"].upper(): status_list.append("ALT")
                elif d["base"] == site["ref"].upper(): status_list.append("REF")
            if status_list:
                n_alt   = status_list.count("ALT")
                n_total = len(status_list)
                frag_list.append({
                    "obs_dict":   obs,
                    "is_cis_alt": (n_total >= 2 and n_alt == n_total),
                    "is_trans":   (n_total >= 2 and 0 < n_alt < n_total),
                })
        return pd.DataFrame(frag_list)

    short_frag_df = _build_frag_df(short_fragments)
    long_frag_df  = _build_frag_df(long_fragments)

    bin_id = f"{chrom}:{start}-{end}"

    # ── 두 그룹 각각 요약 ────────────────────────────────────────────
    short_res = _summarize_fragment_group(
        short_frag_df, bin_sites, bin_id,
        short_raw_count, short_breadth, thresholds,
        short_mapped_bases, short_nm_errors,
        prefix="short_",
    )
    long_res = _summarize_fragment_group(
        long_frag_df, bin_sites, bin_id,
        long_raw_count, long_breadth, thresholds,
        long_mapped_bases, long_nm_errors,
        prefix="long_",
    )

    # ── 단일 summary_df 조합 ─────────────────────────────────────────
    row = {"bin_id": bin_id, **short_res, **long_res}
    summary_df = pd.DataFrame([row])

    for col, dtype in BIN_SUMMARY_SCHEMA.items():
        if col in summary_df.columns:
            summary_df[col] = summary_df[col].astype(dtype)

    # position report는 short+long 합산 (BAF 분석용), 없으면 빈 DataFrame
    # 필요 시 short_frag_df / long_frag_df 기반 report_df도 분리 가능
    return pd.DataFrame(), summary_df


# -------------------------------------------------------------------------
# Quality filter (기존과 동일, short_raw_count 기준 사용)
# -------------------------------------------------------------------------
def apply_low_quality_filter(df, min_depth=None, min_coverage=None):
    if df is None or df.empty:
        return pd.DataFrame()

    if min_depth    is None: min_depth    = CFG["min_depth"]
    if min_coverage is None: min_coverage = CFG["min_coverage"]

    initial_count  = len(df)
    sex_chrom_mask = df["chrom"].isin(["chrX", "chrY"])
    df_auto        = df[~sex_chrom_mask].copy()
    df_sex         = df[ sex_chrom_mask].copy()

    # short_raw_count + long_raw_count 합산으로 필터 (총 depth 기준)
    total_count = df_auto["short_raw_count"] + df_auto["long_raw_count"]
    breadth     = df_auto[["short_breadth_ratio", "long_breadth_ratio"]].max(axis=1)

    df_auto_filtered = df_auto[
        (total_count >= min_depth) &
        (breadth     >= min_coverage)
    ].copy()

    df_filtered = pd.concat([df_auto_filtered, df_sex], ignore_index=True)

    def _chrom_key(c):
        c_str = str(c).replace("chr", "")
        if c_str == "X": return 23
        if c_str == "Y": return 24
        try: return int(c_str)
        except: return 99

    df_filtered["chrom_sort_key"] = df_filtered["chrom"].apply(_chrom_key)
    df_filtered = (
        df_filtered
        .sort_values(["chrom_sort_key", "start"])
        .drop(columns=["chrom_sort_key"])
        .reset_index(drop=True)
    )

    log(f"Quality Filter: Removed {initial_count - len(df_filtered)} bins.")
    return df_filtered


# -------------------------------------------------------------------------
# fetch_sites_with_handle (기존 동일)
# -------------------------------------------------------------------------
def fetch_sites_with_handle(vcf_handle, chrom, start, end):
    sites = []
    try:
        for record in vcf_handle.fetch(chrom, start, end):
            kova      = record.info.get("KOVA_AF")
            kova_af   = kova[0] if kova and isinstance(kova, (list, tuple)) else (kova if kova else 0.0)
            gnomad    = record.info.get("GNOMAD_AF")
            gnomad_af = gnomad[0] if gnomad and isinstance(gnomad, (list, tuple)) else (gnomad if gnomad else 0.0)
            pop_af    = kova_af if kova_af > 0 else gnomad_af
            if pop_af > 0:
                sites.append({
                    "chrom": record.chrom, "pos": record.pos,
                    "ref": record.ref, "alt": record.alts[0], "pop_af": pop_af,
                })
    except:
        pass

    df = pd.DataFrame(sites)
    if df.empty:
        return pd.DataFrame(columns=["chrom", "pos", "ref", "alt", "pop_af"])
    return df


# -------------------------------------------------------------------------
# 병렬 워커 (기존 구조 유지, parquet 단일 파일)
# -------------------------------------------------------------------------
def _parallel_chrom_worker(chrom_df, bam_path, vcf_file, min_mapq, thresholds, tmp_dir):
    chrom = chrom_df["chrom"].iloc[0]
    if not os.path.exists(vcf_file):
        return None, None

    tmp_sum_path = os.path.join(tmp_dir, f"tmp_sum_{chrom}.parquet")

    local_bin_summaries = []
    min_baseq = thresholds.get("min_baseq", _BIN_EXTRACTION_CFG["min_baseq"])

    with pysam.AlignmentFile(bam_path, "rb", threads=1) as bam_handle, \
         pysam.VariantFile(vcf_file) as vcf_handle:

        for r in chrom_df.itertuples():
            site_df = fetch_sites_with_handle(vcf_handle, r.chrom, r.start, r.end)

            _, summary_df = process_bin_with_handles(
                bam_handle,
                {"chrom": r.chrom, "start": r.start, "end": r.end},
                site_df,
                min_mapq=min_mapq,
                min_baseq=min_baseq,
                thresholds=thresholds,
            )

            if not summary_df.empty:
                summary_dict = summary_df.iloc[0].to_dict()
                summary_dict.update({"chrom": r.chrom, "start": r.start, "end": r.end})
                local_bin_summaries.append(summary_dict)

    if local_bin_summaries:
        pd.DataFrame(local_bin_summaries).to_parquet(tmp_sum_path, engine="pyarrow", index=False)

    return (
        tmp_sum_path if os.path.exists(tmp_sum_path) else None,
        None,  # position report 미사용 (필요 시 복구)
    )


# -------------------------------------------------------------------------
# GC correction (기존과 동일, count_col 파라미터로 short/long 선택)
# -------------------------------------------------------------------------
def gc_correct_lowess_robust(
    df,
    count_col="short_raw_count",   # "short_raw_count" | "long_raw_count"
    gc_col="gc",
    frac=0.5,
    pseudocount=0.25,
    min_count_for_fit=1,
    gc_range=(0.25, 0.75),
    mappability_col=None,
    blacklist_col=None,
):
    """
    short / long 각각 호출해서 GC 보정.
    count_col 파라미터로 short_raw_count 또는 long_raw_count 지정.

    output columns:
      log2_corrected  : GC 보정된 log2 normalized count
      gc_fit          : GC별 fitted value
      gc_valid_for_fit: LOWESS fit에 사용된 bin 여부
      gc_correction   : correction term (clipped)
    """
    df = df.copy()

    counts = df[count_col].astype(float).values
    gc     = df[gc_col].astype(float).values

    total = np.nansum(counts)
    if total <= 0:
        df["log2_corrected"]   = np.nan
        df["gc_fit"]           = np.nan
        df["gc_valid_for_fit"] = False
        return df, None

    norm_count = counts / total * 1_000_000
    y = np.log2(norm_count + pseudocount)

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

    if valid.sum() < 100:
        df["log2_corrected"]   = y
        df["gc_fit"]           = np.nan
        df["gc_valid_for_fit"] = valid
        return df, (gc[valid], y[valid], None)

    fit_valid = lowess(y[valid], gc[valid], frac=frac, return_sorted=False)

    order      = np.argsort(gc[valid])
    gc_sorted  = gc[valid][order]
    fit_sorted = fit_valid[order]

    fit_all  = np.interp(gc, gc_sorted, fit_sorted,
                         left=fit_sorted[0], right=fit_sorted[-1])
    baseline = np.nanmedian(fit_valid)

    correction = fit_all - baseline
    correction = np.clip(correction, -0.35, 0.35)

    df["log2_corrected"]   = y - correction
    df["gc_fit"]           = fit_all
    df["gc_valid_for_fit"] = valid
    df["gc_correction"]    = correction

    return df, (gc[valid], y[valid], fit_valid)