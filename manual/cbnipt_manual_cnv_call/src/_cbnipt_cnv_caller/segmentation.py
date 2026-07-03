import numpy as np
import pandas as pd
import ruptures as rpt
import os, sys
 
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

from utils import log

def rolling_micro_cnv_segmentation(
    df,
    signal_col="normalized_count",
    chrom_col="chrom",
    start_col="start",
    end_col="end",
    window=3,
    gain_threshold=1.25,
    loss_threshold=0.75,
    min_bins=5,
    min_segment_bp=500_000,
    max_seg_mad=0.20,
):
    import numpy as np
    import pandas as pd

    segments = []

    for chrom, g in df.groupby(chrom_col, sort=False):
        g = g.sort_values(start_col).reset_index(drop=True).copy()

        x = (
            g[signal_col]
            .replace([np.inf, -np.inf], np.nan)
            .astype(float)
        )

        smooth = x.rolling(
            window=window,
            center=True,
            min_periods=1
        ).median()

        g["smooth_signal"] = smooth
        g["bin_state"] = "NEUT"
        print(g)
        g.loc[g["smooth_signal"] >= gain_threshold, "bin_state"] = "AMP"
        g.loc[g["smooth_signal"] <= loss_threshold, "bin_state"] = "DEL"
        g.loc[g["smooth_signal"].isna(), "bin_state"] = "INVALID"

        states = g["bin_state"].values
        start_i = 0
        print(g["bin_state"].value_counts().to_dict())

        for i in range(1, len(g) + 1):
            if i == len(g) or states[i] != states[start_i]:
                state = states[start_i]
                end_i = i - 1

                if state in ["AMP", "DEL"]:
                    seg = g.iloc[start_i:end_i + 1].copy()

                    seg_start = int(seg[start_col].iloc[0])
                    seg_end = int(seg[end_col].iloc[-1])
                    seg_len = seg_end - seg_start

                    vals = seg["smooth_signal"].dropna().values

                    if len(vals) == 0:
                        start_i = i
                        continue

                    seg_median = float(np.nanmedian(vals))
                    seg_mad = float(np.nanmedian(np.abs(vals - seg_median)))

                    if (
                        len(seg) >= min_bins and
                        seg_len >= min_segment_bp and
                        seg_mad <= max_seg_mad
                    ):
                        seg_mean = float(np.log2(seg_median + 1e-9))
                        cn_signal = float(2.0 * seg_median)

                        segments.append({
                            "chrom": chrom,
                            "start": seg_start,
                            "end": seg_end,
                            "n_bins": int(len(seg)),
                            "n_valid_bins": int(len(vals)),
                            "length": int(seg_len),
                            "seg_mean": seg_mean,
                            "seg_mad": seg_mad,
                            "median_norm": seg_median,
                            "copy_number_signal": cn_signal,
                            "copy_number": int(round(cn_signal)),
                            "cnv_call": state,
                            "method": "rolling_micro_cnv",
                        })

                start_i = i

    return pd.DataFrame(segments)

def segment_one_cell(
    meta,
    df_signals,
    signal_col="log2_chrom_norm",
    method="pelt",
    model="l1",
    penalty=8.0,
    n_bkps=None,
    min_size=4,
    jump=1,
    min_bins=4,
    min_segment_bp=400_000,
    gain_log2=0.35,
    loss_log2=-0.45,
):
    """
    Chromosome별 change-point segmentation.

    핵심 수정:
      1. valid index 기준으로 segmentation 수행
      2. breakpoint를 다시 원래 chromosome-local index로 변환
      3. 변화 포인트마다 segment median 값을 따로 계산
      4. segment별 copy_number_signal / cnv_call까지 같이 생성

    Parameters
    ----------
    signal_col:
        기본값 "log2_chrom_norm"
        정상 diploid 기준 0.0이어야 함.

    method:
        "pelt" 또는 "dynp"

    model:
        ruptures model.
        추천:
          - "l1": WGA noise에 상대적으로 robust
          - "l2": 평균 변화 탐지
          - "rbf": 분포 변화 탐지, 조금 더 민감함

    penalty:
        PELT 사용 시 penalty.

    n_bkps:
        Dynp 사용 시 chromosome당 breakpoint 개수.

    gain_log2:
        AMP threshold.
        log2(1.30) ≈ 0.38

    loss_log2:
        DEL threshold.
        log2(0.70) ≈ -0.51
    """

    if method == "pelt" and penalty is None:
        raise ValueError("Penalty value must be provided for PELT algorithm.")

    if method == "dynp" and n_bkps is None:
        raise ValueError("n_bkps must be provided for Dynp algorithm.")

    segments = []

    chrom_list = meta["chrom"].dropna().unique()

    for chrom in chrom_list:
        chrom_idx = meta.index[meta["chrom"] == chrom]

        y = (
            df_signals.loc[chrom_idx, signal_col]
            .replace([np.inf, -np.inf], np.nan)
            .astype(float)
            .values
            .flatten()
        )

        valid = np.isfinite(y)

        if valid.sum() < max(10, min_size * 2):
            continue

        # valid bin만 segmentation에 사용
        y_valid = y[valid]
        y_v = y_valid.reshape(-1, 1)

        # valid array index → chromosome-local index mapping
        valid_local_idx = np.where(valid)[0]

        # --------------------------------------------------
        # 1. breakpoint detection
        # --------------------------------------------------
        if method == "pelt":
            algo = rpt.Pelt(
                model=model,
                min_size=min_size,
                jump=jump
            ).fit(y_v)

            bkps = algo.predict(pen=penalty)

        elif method == "dynp":
            max_bkps = max(0, len(y_v) // min_size - 1)
            use_n_bkps = min(int(n_bkps), max_bkps)

            if use_n_bkps <= 0:
                bkps = [len(y_v)]
            else:
                algo = rpt.Dynp(
                    model=model,
                    min_size=min_size,
                    jump=jump
                ).fit(y_v)

                bkps = algo.predict(n_bkps=use_n_bkps)

        else:
            raise ValueError("method must be either 'pelt' or 'dynp'.")

        # --------------------------------------------------
        # 2. breakpoint별 segment 생성
        # --------------------------------------------------
        valid_start = 0

        for bp in bkps:
            valid_end = bp - 1

            if valid_end < valid_start:
                valid_start = bp
                continue

            # valid index 기준 segment
            seg_valid_vals = y_valid[valid_start:bp]

            if len(seg_valid_vals) == 0:
                valid_start = bp
                continue

            # 원래 chromosome-local index로 변환
            local_start = valid_local_idx[valid_start]
            local_end = valid_local_idx[valid_end]

            full_seg_idx = chrom_idx[local_start:local_end + 1]

            if len(full_seg_idx) < min_bins:
                valid_start = bp
                continue

            seg_start = int(meta.loc[full_seg_idx[0], "start"])
            seg_end = int(meta.loc[full_seg_idx[-1], "end"])
            seg_len = seg_end - seg_start

            if min_segment_bp and seg_len < min_segment_bp:
                valid_start = bp
                continue

            # segment 값은 valid bin 기준 median 사용
            seg_mean = float(np.nanmedian(seg_valid_vals))
            seg_std = float(np.nanstd(seg_valid_vals))
            seg_mad = float(
                np.nanmedian(np.abs(seg_valid_vals - np.nanmedian(seg_valid_vals)))
            )

            # log2_chrom_norm 기준이면:
            # normalized_count = 2^seg_mean
            # copy_number_signal = 2 * normalized_count
            median_norm = float(np.exp2(seg_mean))
            copy_number_signal = float(2.0 * median_norm)
            copy_number = int(round(copy_number_signal))

            if seg_mean >= gain_log2:
                cnv_call = "AMP"
            elif seg_mean <= loss_log2:
                cnv_call = "DEL"
            else:
                cnv_call = "NEUT"
            print({
                "chrom": chrom,
                "start": seg_start,
                "end": seg_end,
                "n_bins": int(len(full_seg_idx)),
                "n_valid_bins": int(len(seg_valid_vals)),
                "seg_mean": seg_mean,
                "seg_std": seg_std,
                "seg_mad": seg_mad,
                "median_norm": median_norm,
                "copy_number_signal": copy_number_signal,
                "copy_number": copy_number,
                "cnv_call": cnv_call,
                "method": method,
                "model": model,
                "penalty": penalty if method == "pelt" else np.nan,
                "n_bkps": n_bkps if method == "dynp" else np.nan,
            })
            segments.append({
                "chrom": chrom,
                "start": seg_start,
                "end": seg_end,
                "n_bins": int(len(full_seg_idx)),
                "n_valid_bins": int(len(seg_valid_vals)),
                "seg_mean": seg_mean,
                "seg_std": seg_std,
                "seg_mad": seg_mad,
                "median_norm": median_norm,
                "copy_number_signal": copy_number_signal,
                "copy_number": copy_number,
                "cnv_call": cnv_call,
                "method": method,
                "model": model,
                "penalty": penalty if method == "pelt" else np.nan,
                "n_bkps": n_bkps if method == "dynp" else np.nan,
            })

            valid_start = bp

    if not segments:
        try:
            log("Warning: No segments were generated. Check input signal or parameters.")
        except NameError:
            print("Warning: No segments were generated. Check input signal or parameters.")

        return pd.DataFrame(columns=[
            "chrom", "start", "end", "n_bins", "n_valid_bins",
            "seg_mean", "seg_std", "seg_mad", "median_norm",
            "copy_number_signal", "copy_number", "cnv_call",
            "method", "model", "penalty", "n_bkps",
        ])

    return pd.DataFrame(segments)


def assign_cn_state(seg_df, baseline_ploidy=2, sex_tag="UNKNOWN"):
    """
    [SIMPLIFIED & FIXED] 1D seg_mean 기반 상태 할당.
    성별(Sex)에 따라 성염색체(chrX, chrY)의 기본 Ploidy를 동적으로 조절합니다.
    """
    if seg_df is None or seg_df.empty:
        return seg_df
        
    df = seg_df.copy()
    
    # [FIX] seg_median -> seg_mean 컬럼명 통일 확인
    if "seg_mean" not in df.columns:
        raise KeyError("Cannot find 'seg_mean' column in segments dataframe. "
                       "segment_one_cell() output must contain 'seg_mean'.")

    ratios = np.exp2(df["seg_mean"].fillna(0.0))
    
    expected_ploidy = np.full(len(df), baseline_ploidy, dtype=int)
    chrom_array = df["chrom"].values
    
    if sex_tag in ("XY", "Male (XY)", "Male"):
        expected_ploidy[chrom_array == "chrX"] = 1
        expected_ploidy[chrom_array == "chrY"] = 1
    elif sex_tag in ("XX", "Female (XX)", "Female"):
        expected_ploidy[chrom_array == "chrY"] = 0
        
    df["expected_ploidy"] = expected_ploidy
    df["copy_number"] = np.round(ratios * expected_ploidy).astype(int)
    df["copy_number"] = df["copy_number"].clip(lower=0)
    
    conditions = [
        df["copy_number"] > df["expected_ploidy"],
        df["copy_number"] < df["expected_ploidy"]
    ]
    choices = ["AMP", "DEL"]
    df["cnv_call"] = np.select(conditions, choices, default="NEUT")
    
    return df
