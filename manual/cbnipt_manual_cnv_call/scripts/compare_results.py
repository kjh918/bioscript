import os
import glob
import pandas as pd
import numpy as np


def load_segment_file(path, sample_id=None):
    """
    CNV segment TSV file load.
    Required columns:
    chrom, start, end, copy_number, copy_number_signal # [MODIFIED] Added copy_number_signal
    """

    df = pd.read_csv(path, sep="\t")

    # [MODIFIED] Added copy_number_signal to required columns
    required_cols = ["chrom", "start", "end", "copy_number", "copy_number_signal"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"{path} missing columns: {missing}")

    df = df.copy()
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["copy_number"] = df["copy_number"].astype(int)
    df["copy_number_signal"] = df["copy_number_signal"].astype(float) # [MODIFIED] Cast signal to float

    if sample_id is None:
        sample_id = os.path.basename(path).split(".")[0]

    df["sample_id"] = sample_id

    return df


def find_best_overlap(ref_row, query_df):
    """
    기준 영역 ref_row와 query_df 내 segment 중 overlap 되는 것들을 찾고,
    overlap 길이가 가장 긴 segment를 반환.
    """

    chrom = ref_row["chrom"]
    ref_start = int(ref_row["start"])
    ref_end = int(ref_row["end"])

    sub = query_df[query_df["chrom"] == chrom].copy()

    if sub.empty:
        return None

    # overlap 조건
    sub = sub[
        (sub["start"] < ref_end) &
        (sub["end"] > ref_start)
    ].copy()

    if sub.empty:
        return None

    sub["overlap_start"] = np.maximum(sub["start"], ref_start)
    sub["overlap_end"] = np.minimum(sub["end"], ref_end)
    sub["overlap_bp"] = sub["overlap_end"] - sub["overlap_start"]

    sub = sub[sub["overlap_bp"] > 0].copy()

    if sub.empty:
        return None

    sub["ref_len_bp"] = ref_end - ref_start
    sub["overlap_fraction"] = sub["overlap_bp"] / sub["ref_len_bp"]

    # 가장 많이 겹치는 segment 선택
    best = sub.sort_values(
        ["overlap_bp", "overlap_fraction"],
        ascending=False
    ).iloc[0]

    return best


def compare_to_reference(
    ref_df,
    query_df,
    query_sample_id,
    min_overlap_fraction=0.5
):
    """
    기준 샘플 segment 기준으로 query sample copy_number 및 signal 비교.
    """

    records = []

    for _, ref_row in ref_df.iterrows():
        best = find_best_overlap(ref_row, query_df)

        ref_cn = int(ref_row["copy_number"])
        ref_signal = float(ref_row["copy_number_signal"]) # [MODIFIED] Extract reference signal

        if best is None:
            query_cn = np.nan
            query_signal = np.nan # [MODIFIED] Handle missing signal
            query_call = "NO_OVERLAP"
            cn_diff = np.nan
            is_match = False
            overlap_bp = 0
            overlap_fraction = 0
            query_start = np.nan
            query_end = np.nan
        else:
            query_cn = int(best["copy_number"])
            query_signal = float(best["copy_number_signal"]) # [MODIFIED] Extract query signal
            overlap_bp = int(best["overlap_bp"])
            overlap_fraction = float(best["overlap_fraction"])
            query_start = int(best["start"])
            query_end = int(best["end"])

            if overlap_fraction < min_overlap_fraction:
                query_call = "LOW_OVERLAP"
                is_match = False
            else:
                query_call = "COMPARED"
                is_match = ref_cn == query_cn

            cn_diff = query_cn - ref_cn

        records.append({
            "query_sample_id": query_sample_id,

            "chrom": ref_row["chrom"],
            "ref_start": int(ref_row["start"]),
            "ref_end": int(ref_row["end"]),
            "ref_len_bp": int(ref_row["end"] - ref_row["start"]),
            "ref_copy_number": ref_cn,
            "ref_copy_number_signal": ref_signal, # [MODIFIED] Added to records

            "query_start": query_start,
            "query_end": query_end,
            "query_copy_number": query_cn,
            "query_copy_number_signal": query_signal, # [MODIFIED] Added to records

            "overlap_bp": overlap_bp,
            "overlap_fraction": overlap_fraction,

            "cn_diff": cn_diff,
            "abs_cn_diff": abs(cn_diff) if not pd.isna(cn_diff) else np.nan,

            "is_match": is_match,
            "compare_status": query_call
        })

    return pd.DataFrame(records)


def summarize_comparison(detail_df):
    """
    샘플별 일치율 / 불일치 정도 / signal correlation 요약.
    """

    rows = []

    for sample_id, sub in detail_df.groupby("query_sample_id"):
        comparable = sub[sub["compare_status"] == "COMPARED"].copy()

        n_ref_regions = len(sub)
        n_compared = len(comparable)
        n_no_overlap = (sub["compare_status"] == "NO_OVERLAP").sum()
        n_low_overlap = (sub["compare_status"] == "LOW_OVERLAP").sum()

        if n_compared == 0:
            match_rate = np.nan
            mismatch_rate = np.nan
            mean_abs_cn_diff = np.nan
            median_abs_cn_diff = np.nan
            max_abs_cn_diff = np.nan
            signal_correlation = np.nan # [MODIFIED] Default nan for correlation
            n_match = 0
            n_mismatch = 0
        else:
            n_match = comparable["is_match"].sum()
            n_mismatch = n_compared - n_match

            match_rate = n_match / n_compared
            mismatch_rate = n_mismatch / n_compared
            mean_abs_cn_diff = comparable["abs_cn_diff"].mean()
            median_abs_cn_diff = comparable["abs_cn_diff"].median()
            max_abs_cn_diff = comparable["abs_cn_diff"].max()
            
            # [MODIFIED] Calculate Pearson correlation for copy_number_signal in overlapping regions
            if n_compared > 1:
                signal_correlation = comparable["ref_copy_number_signal"].corr(comparable["query_copy_number_signal"])
            else:
                signal_correlation = np.nan # Correlation requires at least 2 points

        rows.append({
            "query_sample_id": sample_id,
            "n_ref_regions": n_ref_regions,
            "n_compared_regions": n_compared,
            "n_match": n_match,
            "n_mismatch": n_mismatch,
            "n_no_overlap": n_no_overlap,
            "n_low_overlap": n_low_overlap,

            "match_rate": match_rate,
            "mismatch_rate": mismatch_rate,

            "mean_abs_cn_diff": mean_abs_cn_diff,
            "median_abs_cn_diff": median_abs_cn_diff,
            "max_abs_cn_diff": max_abs_cn_diff,
            
            "signal_correlation": signal_correlation # [MODIFIED] Added correlation to summary
        })

    summary_df = pd.DataFrame(rows)

    if not summary_df.empty:
        summary_df = summary_df.sort_values(
            ["match_rate", "mean_abs_cn_diff"],
            ascending=[False, True]
        )

    return summary_df


def compare_all_samples_to_reference(
    segment_files,
    reference_sample_id,
    output_prefix="cnv_reference_comparison",
    min_overlap_fraction=0.5
):
    """
    segment_files:
      list of segment tsv paths

    reference_sample_id:
      기준 샘플 ID.
      파일명에서 첫 번째 '.' 앞부분과 매칭함.

    output:
      1) {output_prefix}.detail.tsv
      2) {output_prefix}.summary.tsv
    """

    sample_to_path = {
        os.path.basename(p).split(".")[0]: p
        for p in segment_files
    }

    if reference_sample_id not in sample_to_path:
        raise ValueError(
            f"Reference sample not found: {reference_sample_id}\n"
            f"Available samples: {list(sample_to_path.keys())}"
        )

    ref_path = sample_to_path[reference_sample_id]
    ref_df = load_segment_file(ref_path, sample_id=reference_sample_id)

    all_detail = []

    for sample_id, path in sample_to_path.items():
        if sample_id == reference_sample_id:
            continue

        print(f"[COMPARE] ref={reference_sample_id} vs query={sample_id}")

        query_df = load_segment_file(path, sample_id=sample_id)

        detail_df = compare_to_reference(
            ref_df=ref_df,
            query_df=query_df,
            query_sample_id=sample_id,
            min_overlap_fraction=min_overlap_fraction
        )

        all_detail.append(detail_df)

    if not all_detail:
        raise ValueError("No query samples to compare.")

    detail_df = pd.concat(all_detail, ignore_index=True)
    summary_df = summarize_comparison(detail_df)

    detail_out = f"{output_prefix}.detail.tsv"
    summary_out = f"{output_prefix}.summary.tsv"

    detail_df.to_csv(detail_out, sep="\t", index=False)
    summary_df.to_csv(summary_out, sep="\t", index=False)
    print(summary_df)

    print(f"[DONE] detail  -> {detail_out}")
    print(f"[DONE] summary -> {summary_out}")

    return detail_df, summary_df


def main():
    segment_files = glob.glob("/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Results/temp/*25_03*/*/*.segments.new.tsv")
    segment_files += glob.glob("/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Results/temp/*25_04*/*/*.segments.new.tsv")
    print(f"Found segment files: {(segment_files)}")
    compare_all_samples_to_reference(
        segment_files,
        'cbNIPT_25_03_01',
        output_prefix="cnv_reference_comparison",
        min_overlap_fraction=0.5
    )


if __name__ == "__main__":
    main()