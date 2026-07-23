import os
import re
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


# ============================================================
# Configuration
# ============================================================

BASE_DIR = (
    "/storage/home/jhkim/Projects/cbNIPT/"
    "260423-GCX-cbNIPT-ManualMethod/"
    "Results/downsampling/manaual_analysis"
)

RESULTS_DIR = os.path.join(BASE_DIR, "results")

EXPECTED_REPLICATES = ["1234", "1235", "1236"]

DIAGNOSIS_SCORE = {
    "NEG": 0,
    "SUS": 1,
    "POS": 2,
}


# ============================================================
# Utility functions
# ============================================================

def map_marker_diagnosis(value):
    """
    clinical_report.tsv의 개별 마커 DIAGNOSIS를 표준화합니다.

    반환:
        POSITIVE
        SUSPICIOUS
        NEGATIVE
    """
    value = str(value).strip().upper()

    if "POSITIVE" in value:
        return "POSITIVE"

    if "SUSPICIOUS" in value:
        return "SUSPICIOUS"

    return "NEGATIVE"


def parse_sample_name(sample_name):
    """
    샘플 이름을 분리합니다.

    입력 예:
        cbNIPT_24_04_04_0.5X_s1234

    반환:
        BASE_SAMPLE_ID = cbNIPT_24_04_04
        DEPTH          = 0.5X
        REP_ID         = 1234
    """
    sample_name = os.path.basename(sample_name)

    if sample_name.endswith(".clinical_report.tsv"):
        sample_name = sample_name.replace(".clinical_report.tsv", "")

    pattern = re.compile(
        r"^(?P<base_sample>.+?)_"
        r"(?P<depth>[0-9]+(?:\.[0-9]+)?X)_"
        r"s(?P<rep_id>[0-9]+)$",
        flags=re.IGNORECASE,
    )

    match = pattern.match(sample_name)

    if not match:
        raise ValueError(
            f"Invalid sample name: {sample_name}\n"
            "Expected format: cbNIPT_24_04_04_0.5X_s1234"
        )

    return {
        "BASE_SAMPLE_ID": match.group("base_sample"),
        "DEPTH": match.group("depth").upper(),
        "REP_ID": match.group("rep_id"),
    }


def depth_to_float(depth):
    """
    depth 문자열을 정렬용 숫자로 변환합니다.

    예:
        0.1X -> 0.1
        0.5X -> 0.5
        1X   -> 1.0
    """
    try:
        return float(str(depth).upper().replace("X", ""))
    except (TypeError, ValueError):
        return float("inf")


def sort_by_depth(dataframe, additional_columns=None):
    """
    DataFrame을 depth 숫자 기준으로 정렬합니다.
    """
    dataframe = dataframe.copy()
    dataframe["_DEPTH_VALUE"] = dataframe["DEPTH"].apply(depth_to_float)

    sort_columns = ["_DEPTH_VALUE"]

    if additional_columns:
        sort_columns.extend(additional_columns)

    dataframe = (
        dataframe
        .sort_values(sort_columns)
        .drop(columns="_DEPTH_VALUE")
        .reset_index(drop=True)
    )

    return dataframe


# ============================================================
# Data loading
# ============================================================

def load_clinical_reports(search_pattern):
    """
    모든 clinical_report.tsv 파일을 읽고 하나의 long-format DataFrame으로
    합칩니다.
    """
    tsv_files = sorted(
        glob.glob(search_pattern, recursive=True)
    )

    if not tsv_files:
        raise FileNotFoundError(
            f"No TSV files found: {search_pattern}"
        )

    required_columns = [
        "SYNDROME",
        "FEATURE_NAME",
        "FEATURE_TYPE",
        "FEAT_RANK",
        "DIAGNOSIS",
    ]

    all_data = []

    for file_path in tsv_files:
        file_name = os.path.basename(file_path)

        try:
            sample_info = parse_sample_name(file_name)
        except ValueError as error:
            print(f"[SKIP] {file_path}")
            print(f"       {error}")
            continue

        try:
            df = pd.read_csv(
                file_path,
                sep="\t",
                dtype=str,
            )
        except Exception as error:
            print(f"[SKIP] Failed to read: {file_path}")
            print(f"       {error}")
            continue

        missing_columns = [
            column
            for column in required_columns
            if column not in df.columns
        ]

        if missing_columns:
            print(
                f"[SKIP] {file_path}: "
                f"missing columns={missing_columns}"
            )
            continue

        df = df.copy()

        # 성별 판정 행 제외
        df = df[
            df["SYNDROME"].astype(str).str.upper()
            != "NIPT_SEX"
        ].copy()

        if df.empty:
            print(f"[SKIP] No syndrome rows: {file_path}")
            continue

        # FEAT_RANK 숫자 변환
        df["FEAT_RANK"] = pd.to_numeric(
            df["FEAT_RANK"],
            errors="coerce",
        )

        invalid_rank_count = df["FEAT_RANK"].isna().sum()

        if invalid_rank_count > 0:
            print(
                f"[WARN] {file_path}: "
                f"{invalid_rank_count} rows have invalid FEAT_RANK"
            )

        df = df.dropna(
            subset=[
                "SYNDROME",
                "FEATURE_NAME",
                "FEAT_RANK",
            ]
        ).copy()

        if df.empty:
            print(f"[SKIP] No valid marker rows: {file_path}")
            continue

        df["FEAT_RANK"] = df["FEAT_RANK"].astype(int)

        df["MARKER_DIAGNOSIS"] = (
            df["DIAGNOSIS"]
            .apply(map_marker_diagnosis)
        )

        full_sample_id = file_name.replace(
            ".clinical_report.tsv",
            "",
        )

        df["FULL_SAMPLE_ID"] = full_sample_id
        df["BASE_SAMPLE_ID"] = sample_info["BASE_SAMPLE_ID"]
        df["DEPTH"] = sample_info["DEPTH"]
        df["REP_ID"] = sample_info["REP_ID"]
        df["SOURCE_FILE"] = file_path

        all_data.append(df)

    if not all_data:
        raise ValueError(
            "No valid clinical report data was parsed."
        )

    combined_df = pd.concat(
        all_data,
        ignore_index=True,
    )

    combined_df["_DEPTH_VALUE"] = (
        combined_df["DEPTH"].apply(depth_to_float)
    )

    combined_df = (
        combined_df
        .sort_values([
            "BASE_SAMPLE_ID",
            "_DEPTH_VALUE",
            "REP_ID",
            "SYNDROME",
            "FEAT_RANK",
            "FEATURE_NAME",
        ])
        .drop(columns="_DEPTH_VALUE")
        .reset_index(drop=True)
    )

    return combined_df


# ============================================================
# Disease-level classification
# ============================================================

def classify_syndrome(marker_group):
    """
    하나의 샘플 × depth × replicate × 질환에 대해 최종 질환 판정을 합니다.

    판정 규칙
    ------------------------------------------------------------
    1. FEAT_RANK == 1 마커 중 하나라도 POSITIVE
       -> POS

    2. Rank 1 마커가 Positive가 아니면서,
       FEAT_RANK > 1인 하위 마커가 존재하고
       하위 마커가 모두 POSITIVE
       -> SUS

    3. 나머지
       -> NEG
    ------------------------------------------------------------

    개별 마커의 SUSPICIOUS는 하위 마커의 "모두 Positive" 조건에
    포함되지 않습니다.
    """
    marker_group = marker_group.copy()

    # 동일한 마커가 중복되어 있는 경우 가장 심각한 결과를 사용
    marker_priority = {
        "NEGATIVE": 0,
        "SUSPICIOUS": 1,
        "POSITIVE": 2,
    }

    marker_group["MARKER_SCORE"] = (
        marker_group["MARKER_DIAGNOSIS"]
        .map(marker_priority)
        .fillna(0)
        .astype(int)
    )

    marker_group = (
        marker_group
        .sort_values("MARKER_SCORE")
        .drop_duplicates(
            subset=[
                "FEATURE_NAME",
                "FEATURE_TYPE",
                "FEAT_RANK",
            ],
            keep="last",
        )
    )

    rank1 = marker_group[
        marker_group["FEAT_RANK"] == 1
    ].copy()

    lower_rank = marker_group[
        marker_group["FEAT_RANK"] > 1
    ].copy()

    rank1_total_count = len(rank1)
    lower_rank_total_count = len(lower_rank)

    rank1_positive_count = int(
        (
            rank1["MARKER_DIAGNOSIS"]
            == "POSITIVE"
        ).sum()
    )

    lower_rank_positive_count = int(
        (
            lower_rank["MARKER_DIAGNOSIS"]
            == "POSITIVE"
        ).sum()
    )

    # Rank 1 중 하나라도 Positive
    if rank1_positive_count > 0:
        final_diagnosis = "POS"
        decision_reason = (
            "At least one FEAT_RANK=1 marker is positive"
        )

    # Rank 1은 양성이 아니며 하위 마커가 모두 Positive
    elif (
        lower_rank_total_count > 0
        and lower_rank_positive_count
        == lower_rank_total_count
    ):
        final_diagnosis = "SUS"
        decision_reason = (
            "All lower-rank markers are positive "
            "without a positive rank-1 marker"
        )

    else:
        final_diagnosis = "NEG"
        decision_reason = (
            "Rank-1 markers are not positive and "
            "not all lower-rank markers are positive"
        )

    def make_marker_result_text(data):
        if data.empty:
            return ""

        return ";".join(
            data.apply(
                lambda row: (
                    f"{row['FEATURE_NAME']}"
                    f"(rank={row['FEAT_RANK']},"
                    f"type={row['FEATURE_TYPE']})"
                    f":{row['MARKER_DIAGNOSIS']}"
                ),
                axis=1,
            ).tolist()
        )

    return pd.Series({
        "SYNDROME_DIAGNOSIS": final_diagnosis,
        "SYNDROME_SCORE": DIAGNOSIS_SCORE[final_diagnosis],
        "DECISION_REASON": decision_reason,

        "RANK1_MARKER_COUNT": rank1_total_count,
        "RANK1_POSITIVE_COUNT": rank1_positive_count,

        "LOWER_RANK_MARKER_COUNT": lower_rank_total_count,
        "LOWER_RANK_POSITIVE_COUNT": lower_rank_positive_count,

        "RANK1_RESULTS": make_marker_result_text(rank1),
        "LOWER_RANK_RESULTS": make_marker_result_text(lower_rank),
    })


def create_replicate_syndrome_results(combined_df):
    """
    replicate별 질환 최종 판정표를 생성합니다.

    행 단위:
        BASE_SAMPLE_ID × DEPTH × REP_ID × SYNDROME
    """
    group_columns = [
        "BASE_SAMPLE_ID",
        "DEPTH",
        "REP_ID",
        "FULL_SAMPLE_ID",
        "SYNDROME",
    ]

    syndrome_result = (
        combined_df
        .groupby(
            group_columns,
            dropna=False,
            sort=False,
        )
        .apply(
            classify_syndrome,
            include_groups=False,
        )
        .reset_index()
    )

    syndrome_result["_DEPTH_VALUE"] = (
        syndrome_result["DEPTH"]
        .apply(depth_to_float)
    )

    syndrome_result["_REP_VALUE"] = pd.to_numeric(
        syndrome_result["REP_ID"],
        errors="coerce",
    )

    syndrome_result = (
        syndrome_result
        .sort_values([
            "BASE_SAMPLE_ID",
            "_DEPTH_VALUE",
            "_REP_VALUE",
            "SYNDROME",
        ])
        .drop(
            columns=[
                "_DEPTH_VALUE",
                "_REP_VALUE",
            ]
        )
        .reset_index(drop=True)
    )

    return syndrome_result


# ============================================================
# Replicate aggregation
# ============================================================

def summarize_sample_depth_results(replicate_syndrome_df):
    """
    동일한 원본 샘플, depth, 질환에서 replicate 결과를 집계합니다.

    주요 결과:
        POS_REP_COUNT
        SUS_REP_COUNT
        NEG_REP_COUNT
        REP_1234
        REP_1235
        REP_1236
        POS_RESULT
    """
    data = replicate_syndrome_df.copy()

    data["IS_POS"] = (
        data["SYNDROME_DIAGNOSIS"] == "POS"
    ).astype(int)

    data["IS_SUS"] = (
        data["SYNDROME_DIAGNOSIS"] == "SUS"
    ).astype(int)

    data["IS_NEG"] = (
        data["SYNDROME_DIAGNOSIS"] == "NEG"
    ).astype(int)

    group_columns = [
        "BASE_SAMPLE_ID",
        "DEPTH",
        "SYNDROME",
    ]

    summary = (
        data.groupby(
            group_columns,
            as_index=False,
            dropna=False,
        )
        .agg(
            POS_REP_COUNT=("IS_POS", "sum"),
            SUS_REP_COUNT=("IS_SUS", "sum"),
            NEG_REP_COUNT=("IS_NEG", "sum"),
            TOTAL_REP_COUNT=("REP_ID", "nunique"),
            AVAILABLE_REPLICATES=(
                "REP_ID",
                lambda values: ",".join(
                    sorted(
                        set(values.astype(str)),
                        key=lambda value: int(value),
                    )
                ),
            ),
        )
    )

    replicate_wide = data.pivot_table(
        index=group_columns,
        columns="REP_ID",
        values="SYNDROME_DIAGNOSIS",
        aggfunc="first",
    ).reset_index()

    replicate_wide.columns.name = None

    rename_map = {}

    for column in replicate_wide.columns:
        if column not in group_columns:
            rename_map[column] = f"REP_{column}"

    replicate_wide = replicate_wide.rename(
        columns=rename_map
    )

    summary = summary.merge(
        replicate_wide,
        on=group_columns,
        how="left",
    )

    for rep_id in EXPECTED_REPLICATES:
        rep_column = f"REP_{rep_id}"

        if rep_column not in summary.columns:
            summary[rep_column] = "MISSING"
        else:
            summary[rep_column] = (
                summary[rep_column]
                .fillna("MISSING")
            )

    summary["MISSING_REPLICATES"] = (
        summary["AVAILABLE_REPLICATES"]
        .apply(
            lambda value: ",".join(
                sorted(
                    set(EXPECTED_REPLICATES)
                    - set(str(value).split(",")),
                    key=int,
                )
            )
        )
    )

    summary["REP_RESULT"] = summary.apply(
        lambda row: "/".join(
            [
                str(row[f"REP_{rep_id}"])
                for rep_id in EXPECTED_REPLICATES
            ]
        ),
        axis=1,
    )

    summary["POS_RESULT"] = (
        summary["POS_REP_COUNT"].astype(str)
        + "/"
        + summary["TOTAL_REP_COUNT"].astype(str)
    )

    summary["SUS_RESULT"] = (
        summary["SUS_REP_COUNT"].astype(str)
        + "/"
        + summary["TOTAL_REP_COUNT"].astype(str)
    )

    def classify_replicate_consensus(row):
        pos_count = int(row["POS_REP_COUNT"])
        sus_count = int(row["SUS_REP_COUNT"])
        total_count = int(row["TOTAL_REP_COUNT"])

        if total_count == 0:
            return "NO_DATA"

        if pos_count == total_count:
            return "ALL_POS"

        if pos_count >= 2:
            return "MAJORITY_POS"

        if pos_count == 1:
            return "SINGLE_POS"

        if sus_count == total_count:
            return "ALL_SUS"

        if sus_count > 0:
            return "SUS_DETECTED"

        return "ALL_NEG"

    summary["REPLICATE_CONSENSUS"] = summary.apply(
        classify_replicate_consensus,
        axis=1,
    )

    summary["_DEPTH_VALUE"] = (
        summary["DEPTH"].apply(depth_to_float)
    )

    summary = (
        summary
        .sort_values([
            "SYNDROME",
            "_DEPTH_VALUE",
            "BASE_SAMPLE_ID",
        ])
        .drop(columns="_DEPTH_VALUE")
        .reset_index(drop=True)
    )

    return summary


# ============================================================
# Positive sample lists
# ============================================================

def create_depth_positive_sample_detail(
    sample_depth_summary,
    min_positive_replicates=1,
):
    """
    특정 기준 이상의 POS replicate를 가진 샘플만 추출합니다.

    예:
        min_positive_replicates=1
            replicate 중 하나라도 POS

        min_positive_replicates=2
            replicate 중 최소 2개가 POS

        min_positive_replicates=3
            replicate 3개 모두 POS
    """
    positive_df = sample_depth_summary[
        sample_depth_summary["POS_REP_COUNT"]
        >= min_positive_replicates
    ].copy()

    positive_df["_DEPTH_VALUE"] = (
        positive_df["DEPTH"].apply(depth_to_float)
    )

    positive_df = (
        positive_df
        .sort_values([
            "_DEPTH_VALUE",
            "SYNDROME",
            "BASE_SAMPLE_ID",
        ])
        .drop(columns="_DEPTH_VALUE")
        .reset_index(drop=True)
    )

    return positive_df


def create_depth_positive_sample_list(
    sample_depth_summary,
    min_positive_replicates=1,
):
    """
    depth와 질환별 Positive 샘플 목록을 요약합니다.
    """
    positive_df = create_depth_positive_sample_detail(
        sample_depth_summary,
        min_positive_replicates=min_positive_replicates,
    )

    output_columns = [
        "DEPTH",
        "SYNDROME",
        "POSITIVE_SAMPLE_COUNT",
        "POSITIVE_SAMPLE_IDS",
        "SAMPLE_POS_REP_RESULTS",
    ]

    if positive_df.empty:
        return pd.DataFrame(columns=output_columns)

    grouped_rows = []

    for (depth, syndrome), group in positive_df.groupby(
        ["DEPTH", "SYNDROME"],
        sort=False,
    ):
        group = group.sort_values("BASE_SAMPLE_ID")

        sample_ids = sorted(
            group["BASE_SAMPLE_ID"]
            .astype(str)
            .unique()
            .tolist()
        )

        result_text = "; ".join(
            group.apply(
                lambda row: (
                    f"{row['BASE_SAMPLE_ID']}="
                    f"{row['POS_RESULT']}"
                    f"[{row['REP_RESULT']}]"
                ),
                axis=1,
            ).tolist()
        )

        grouped_rows.append({
            "DEPTH": depth,
            "SYNDROME": syndrome,
            "POSITIVE_SAMPLE_COUNT": len(sample_ids),
            "POSITIVE_SAMPLE_IDS": ",".join(sample_ids),
            "SAMPLE_POS_REP_RESULTS": result_text,
        })

    result_df = pd.DataFrame(grouped_rows)

    result_df["_DEPTH_VALUE"] = (
        result_df["DEPTH"].apply(depth_to_float)
    )

    result_df = (
        result_df
        .sort_values([
            "_DEPTH_VALUE",
            "SYNDROME",
        ])
        .drop(columns="_DEPTH_VALUE")
        .reset_index(drop=True)
    )

    return result_df


# ============================================================
# Depth-level disease summary
# ============================================================

def create_depth_syndrome_summary(sample_depth_summary):
    """
    depth와 질환별로 전체 샘플 수 및 POS/SUS/NEG 결과를 요약합니다.

    샘플 단위 분류:
        POS_SAMPLE:
            POS replicate가 1개 이상

        SUS_ONLY_SAMPLE:
            POS replicate는 없고 SUS replicate가 1개 이상

        NEG_SAMPLE:
            POS와 SUS replicate가 모두 없음
    """
    data = sample_depth_summary.copy()

    data["HAS_POS"] = (
        data["POS_REP_COUNT"] > 0
    ).astype(int)

    data["HAS_SUS_ONLY"] = (
        (data["POS_REP_COUNT"] == 0)
        & (data["SUS_REP_COUNT"] > 0)
    ).astype(int)

    data["IS_ALL_NEG"] = (
        (data["POS_REP_COUNT"] == 0)
        & (data["SUS_REP_COUNT"] == 0)
    ).astype(int)

    summary = (
        data.groupby(
            ["DEPTH", "SYNDROME"],
            as_index=False,
        )
        .agg(
            TOTAL_SAMPLE_COUNT=(
                "BASE_SAMPLE_ID",
                "nunique",
            ),
            POSITIVE_SAMPLE_COUNT=(
                "HAS_POS",
                "sum",
            ),
            SUSPICIOUS_ONLY_SAMPLE_COUNT=(
                "HAS_SUS_ONLY",
                "sum",
            ),
            NEGATIVE_SAMPLE_COUNT=(
                "IS_ALL_NEG",
                "sum",
            ),
            TOTAL_POS_REPLICATES=(
                "POS_REP_COUNT",
                "sum",
            ),
            TOTAL_SUS_REPLICATES=(
                "SUS_REP_COUNT",
                "sum",
            ),
            TOTAL_NEG_REPLICATES=(
                "NEG_REP_COUNT",
                "sum",
            ),
        )
    )

    summary["POSITIVE_SAMPLE_RATE"] = (
        summary["POSITIVE_SAMPLE_COUNT"]
        / summary["TOTAL_SAMPLE_COUNT"]
    )

    summary["POSITIVE_SAMPLE_RATE_PERCENT"] = (
        summary["POSITIVE_SAMPLE_RATE"] * 100
    ).round(2)

    summary["_DEPTH_VALUE"] = (
        summary["DEPTH"].apply(depth_to_float)
    )

    summary = (
        summary
        .sort_values([
            "_DEPTH_VALUE",
            "SYNDROME",
        ])
        .drop(columns="_DEPTH_VALUE")
        .reset_index(drop=True)
    )

    return summary


# ============================================================
# Matrices
# ============================================================

def create_disease_positive_count_matrix(sample_depth_summary):
    """
    행:
        BASE_SAMPLE_ID + DEPTH

    열:
        SYNDROME

    값:
        replicate 중 POS 개수
    """
    data = sample_depth_summary.copy()

    data["_DEPTH_VALUE"] = (
        data["DEPTH"].apply(depth_to_float)
    )

    row_order_df = (
        data[
            [
                "BASE_SAMPLE_ID",
                "DEPTH",
                "_DEPTH_VALUE",
            ]
        ]
        .drop_duplicates()
        .sort_values([
            "BASE_SAMPLE_ID",
            "_DEPTH_VALUE",
        ])
    )

    row_order_df["SAMPLE_DEPTH"] = (
        row_order_df["BASE_SAMPLE_ID"]
        + "\n"
        + row_order_df["DEPTH"]
    )

    row_order = row_order_df[
        "SAMPLE_DEPTH"
    ].tolist()

    syndrome_order = sorted(
        data["SYNDROME"]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )

    data["SAMPLE_DEPTH"] = (
        data["BASE_SAMPLE_ID"]
        + "\n"
        + data["DEPTH"]
    )

    matrix = data.pivot_table(
        index="SAMPLE_DEPTH",
        columns="SYNDROME",
        values="POS_REP_COUNT",
        aggfunc="max",
    )

    annotation = data.pivot_table(
        index="SAMPLE_DEPTH",
        columns="SYNDROME",
        values="POS_RESULT",
        aggfunc="first",
    )

    matrix = matrix.reindex(
        index=row_order,
        columns=syndrome_order,
    )

    annotation = annotation.reindex(
        index=row_order,
        columns=syndrome_order,
    ).fillna("-")

    return matrix, annotation


def create_disease_final_status_matrix(sample_depth_summary):
    """
    replicate 집계 결과를 샘플-depth 단위의 최종 질환 상태로 변환합니다.

    최종 상태:
        POS:
            POS replicate가 1개 이상

        SUS:
            POS replicate는 없고 SUS replicate가 1개 이상

        NEG:
            POS와 SUS가 모두 없음
    """
    data = sample_depth_summary.copy()

    def determine_status(row):
        if row["POS_REP_COUNT"] > 0:
            return "POS"

        if row["SUS_REP_COUNT"] > 0:
            return "SUS"

        return "NEG"

    data["FINAL_STATUS"] = data.apply(
        determine_status,
        axis=1,
    )

    data["FINAL_STATUS_SCORE"] = (
        data["FINAL_STATUS"]
        .map(DIAGNOSIS_SCORE)
        .astype(int)
    )

    data["_DEPTH_VALUE"] = (
        data["DEPTH"].apply(depth_to_float)
    )

    row_order_df = (
        data[
            [
                "BASE_SAMPLE_ID",
                "DEPTH",
                "_DEPTH_VALUE",
            ]
        ]
        .drop_duplicates()
        .sort_values([
            "BASE_SAMPLE_ID",
            "_DEPTH_VALUE",
        ])
    )

    row_order_df["SAMPLE_DEPTH"] = (
        row_order_df["BASE_SAMPLE_ID"]
        + "\n"
        + row_order_df["DEPTH"]
    )

    row_order = row_order_df[
        "SAMPLE_DEPTH"
    ].tolist()

    syndrome_order = sorted(
        data["SYNDROME"]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )

    data["SAMPLE_DEPTH"] = (
        data["BASE_SAMPLE_ID"]
        + "\n"
        + data["DEPTH"]
    )

    matrix = data.pivot_table(
        index="SAMPLE_DEPTH",
        columns="SYNDROME",
        values="FINAL_STATUS_SCORE",
        aggfunc="max",
    )

    annotation = data.pivot_table(
        index="SAMPLE_DEPTH",
        columns="SYNDROME",
        values="FINAL_STATUS",
        aggfunc="first",
    )

    matrix = matrix.reindex(
        index=row_order,
        columns=syndrome_order,
    )

    annotation = annotation.reindex(
        index=row_order,
        columns=syndrome_order,
    ).fillna("-")

    return matrix, annotation


# ============================================================
# Plotting
# ============================================================

def plot_disease_positive_count_heatmap(
    matrix,
    annotation,
    output_path,
):
    """
    replicate 3개 중 POS로 판정된 개수를 표시합니다.

    셀:
        0/3
        1/3
        2/3
        3/3
    """
    if matrix.empty:
        print(
            "[WARN] Positive-count matrix is empty. "
            "Heatmap was not created."
        )
        return

    num_diseases = matrix.shape[1]
    num_rows = matrix.shape[0]

    figure_width = max(
        14,
        num_diseases * 0.8,
    )

    figure_height = max(
        8,
        num_rows * 0.55 + 3,
    )

    plt.figure(
        figsize=(
            figure_width,
            figure_height,
        )
    )

    cmap = mcolors.ListedColormap([
        "#f7fbff",
        "#fee391",
        "#fdae6b",
        "#de2d26",
    ])

    bounds = [
        -0.5,
        0.5,
        1.5,
        2.5,
        3.5,
    ]

    norm = mcolors.BoundaryNorm(
        bounds,
        cmap.N,
    )

    ax = sns.heatmap(
        matrix,
        cmap=cmap,
        norm=norm,
        annot=annotation,
        fmt="",
        linewidths=0.5,
        linecolor="lightgray",
        mask=matrix.isna(),
        cbar_kws={
            "ticks": [0, 1, 2, 3],
            "label": "POS replicate count",
            "shrink": 0.6,
        },
    )

    colorbar = ax.collections[0].colorbar
    colorbar.set_ticklabels([
        "0 POS",
        "1 POS",
        "2 POS",
        "3 POS",
    ])

    plt.title(
        "Disease-level Positive Replicate Count by Sample and Depth",
        fontsize=17,
        fontweight="bold",
        pad=20,
    )

    plt.xlabel(
        "Syndrome",
        fontsize=13,
        fontweight="bold",
    )

    plt.ylabel(
        "Sample ID / Depth",
        fontsize=13,
        fontweight="bold",
    )

    plt.xticks(
        rotation=90,
        ha="center",
        fontsize=9,
    )

    plt.yticks(
        rotation=0,
        fontsize=9,
    )

    plt.tight_layout()

    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()

    print(f"[SAVE] {output_path}")


def plot_disease_final_status_heatmap(
    matrix,
    annotation,
    output_path,
):
    """
    샘플-depth별 최종 질환 상태를 표시합니다.

    NEG:
        POS replicate와 SUS replicate 모두 없음

    SUS:
        POS replicate는 없지만 SUS replicate가 존재

    POS:
        POS replicate가 하나 이상 존재
    """
    if matrix.empty:
        print(
            "[WARN] Final-status matrix is empty. "
            "Heatmap was not created."
        )
        return

    num_diseases = matrix.shape[1]
    num_rows = matrix.shape[0]

    figure_width = max(
        14,
        num_diseases * 0.8,
    )

    figure_height = max(
        8,
        num_rows * 0.55 + 3,
    )

    plt.figure(
        figsize=(
            figure_width,
            figure_height,
        )
    )

    cmap = mcolors.ListedColormap([
        "#f7fbff",
        "#ffeda0",
        "#f03b20",
    ])

    bounds = [
        -0.5,
        0.5,
        1.5,
        2.5,
    ]

    norm = mcolors.BoundaryNorm(
        bounds,
        cmap.N,
    )

    ax = sns.heatmap(
        matrix,
        cmap=cmap,
        norm=norm,
        annot=annotation,
        fmt="",
        linewidths=0.5,
        linecolor="lightgray",
        mask=matrix.isna(),
        cbar_kws={
            "ticks": [0, 1, 2],
            "label": "Final disease status",
            "shrink": 0.6,
        },
    )

    colorbar = ax.collections[0].colorbar
    colorbar.set_ticklabels([
        "NEG",
        "SUS",
        "POS",
    ])

    plt.title(
        "Final Disease Status by Sample and Depth",
        fontsize=17,
        fontweight="bold",
        pad=20,
    )

    plt.xlabel(
        "Syndrome",
        fontsize=13,
        fontweight="bold",
    )

    plt.ylabel(
        "Sample ID / Depth",
        fontsize=13,
        fontweight="bold",
    )

    plt.xticks(
        rotation=90,
        ha="center",
        fontsize=9,
    )

    plt.yticks(
        rotation=0,
        fontsize=9,
    )

    plt.tight_layout()

    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()

    print(f"[SAVE] {output_path}")


# ============================================================
# Console summary
# ============================================================

def print_positive_summary(sample_depth_summary):
    """
    콘솔에 depth별 Positive 샘플을 출력합니다.
    """
    positive_df = sample_depth_summary[
        sample_depth_summary["POS_REP_COUNT"] > 0
    ].copy()

    if positive_df.empty:
        print("\nNo positive disease calls found.")
        return

    positive_df["_DEPTH_VALUE"] = (
        positive_df["DEPTH"].apply(depth_to_float)
    )

    positive_df = positive_df.sort_values([
        "_DEPTH_VALUE",
        "SYNDROME",
        "BASE_SAMPLE_ID",
    ])

    print("\n" + "=" * 90)
    print("Positive disease calls by depth")
    print("=" * 90)

    for depth in (
        positive_df[
            ["DEPTH", "_DEPTH_VALUE"]
        ]
        .drop_duplicates()
        .sort_values("_DEPTH_VALUE")["DEPTH"]
    ):
        depth_df = positive_df[
            positive_df["DEPTH"] == depth
        ]

        print(f"\n[Depth: {depth}]")

        for syndrome, syndrome_df in depth_df.groupby(
            "SYNDROME",
            sort=True,
        ):
            print(f"  {syndrome}")

            for _, row in syndrome_df.iterrows():
                print(
                    f"    - {row['BASE_SAMPLE_ID']}: "
                    f"POS {row['POS_RESULT']}, "
                    f"replicates={row['REP_RESULT']}"
                )


# ============================================================
# Main
# ============================================================

def main():
    os.makedirs(
        RESULTS_DIR,
        exist_ok=True,
    )

    search_pattern = os.path.join(
        BASE_DIR,
        "cbNIPT_*_*X_s*",
        "data",
        "*.clinical_report.tsv",
    )

    print("=" * 90)
    print("cbNIPT disease-level replicate analysis")
    print("=" * 90)
    print(f"Search pattern:\n{search_pattern}")

    # --------------------------------------------------------
    # 1. Load all marker-level results
    # --------------------------------------------------------
    print("\n[1/7] Loading clinical reports...")

    combined_df = load_clinical_reports(
        search_pattern
    )

    all_marker_output = os.path.join(
        RESULTS_DIR,
        "01_all_marker_results.tsv",
    )

    combined_df.to_csv(
        all_marker_output,
        sep="\t",
        index=False,
    )

    print(
        f"Loaded marker rows: {len(combined_df):,}"
    )
    print(
        f"Base samples: "
        f"{combined_df['BASE_SAMPLE_ID'].nunique():,}"
    )
    print(
        f"Depths: "
        f"{combined_df['DEPTH'].nunique():,}"
    )
    print(
        f"Replicates: "
        f"{combined_df['FULL_SAMPLE_ID'].nunique():,}"
    )
    print(f"[SAVE] {all_marker_output}")

    # --------------------------------------------------------
    # 2. Determine disease result per replicate
    # --------------------------------------------------------
    print(
        "\n[2/7] Classifying syndrome result "
        "for each replicate..."
    )

    replicate_syndrome_df = (
        create_replicate_syndrome_results(
            combined_df
        )
    )

    replicate_output = os.path.join(
        RESULTS_DIR,
        "02_replicate_syndrome_diagnosis.tsv",
    )

    replicate_syndrome_df.to_csv(
        replicate_output,
        sep="\t",
        index=False,
    )

    print(
        f"Replicate syndrome rows: "
        f"{len(replicate_syndrome_df):,}"
    )
    print(f"[SAVE] {replicate_output}")

    # --------------------------------------------------------
    # 3. Aggregate replicate results by sample and depth
    # --------------------------------------------------------
    print(
        "\n[3/7] Aggregating replicate results "
        "by sample, depth and syndrome..."
    )

    sample_depth_summary = (
        summarize_sample_depth_results(
            replicate_syndrome_df
        )
    )

    sample_depth_output = os.path.join(
        RESULTS_DIR,
        "03_sample_depth_syndrome_summary.tsv",
    )

    sample_depth_summary.to_csv(
        sample_depth_output,
        sep="\t",
        index=False,
    )

    print(f"[SAVE] {sample_depth_output}")

    # --------------------------------------------------------
    # 4. Create positive sample lists
    # --------------------------------------------------------
    print(
        "\n[4/7] Creating positive sample lists..."
    )

    positive_any_detail = (
        create_depth_positive_sample_detail(
            sample_depth_summary,
            min_positive_replicates=1,
        )
    )

    positive_any_detail_output = os.path.join(
        RESULTS_DIR,
        "04_positive_samples_any_rep_detail.tsv",
    )

    positive_any_detail.to_csv(
        positive_any_detail_output,
        sep="\t",
        index=False,
    )

    positive_any_list = (
        create_depth_positive_sample_list(
            sample_depth_summary,
            min_positive_replicates=1,
        )
    )

    positive_any_list_output = os.path.join(
        RESULTS_DIR,
        "05_positive_samples_any_rep_by_depth.tsv",
    )

    positive_any_list.to_csv(
        positive_any_list_output,
        sep="\t",
        index=False,
    )

    positive_2of3_detail = (
        create_depth_positive_sample_detail(
            sample_depth_summary,
            min_positive_replicates=2,
        )
    )

    positive_2of3_detail_output = os.path.join(
        RESULTS_DIR,
        "06_positive_samples_2of3_detail.tsv",
    )

    positive_2of3_detail.to_csv(
        positive_2of3_detail_output,
        sep="\t",
        index=False,
    )

    positive_2of3_list = (
        create_depth_positive_sample_list(
            sample_depth_summary,
            min_positive_replicates=2,
        )
    )

    positive_2of3_list_output = os.path.join(
        RESULTS_DIR,
        "07_positive_samples_2of3_by_depth.tsv",
    )

    positive_2of3_list.to_csv(
        positive_2of3_list_output,
        sep="\t",
        index=False,
    )

    print(f"[SAVE] {positive_any_detail_output}")
    print(f"[SAVE] {positive_any_list_output}")
    print(f"[SAVE] {positive_2of3_detail_output}")
    print(f"[SAVE] {positive_2of3_list_output}")

    # --------------------------------------------------------
    # 5. Create depth-level summary
    # --------------------------------------------------------
    print(
        "\n[5/7] Creating depth-level disease summary..."
    )

    depth_syndrome_summary = (
        create_depth_syndrome_summary(
            sample_depth_summary
        )
    )

    depth_summary_output = os.path.join(
        RESULTS_DIR,
        "08_depth_syndrome_summary.tsv",
    )

    depth_syndrome_summary.to_csv(
        depth_summary_output,
        sep="\t",
        index=False,
    )

    print(f"[SAVE] {depth_summary_output}")

    # --------------------------------------------------------
    # 6. Positive replicate-count heatmap
    # --------------------------------------------------------
    print(
        "\n[6/7] Creating positive replicate-count heatmap..."
    )

    positive_count_matrix, positive_count_annotation = (
        create_disease_positive_count_matrix(
            sample_depth_summary
        )
    )

    positive_count_matrix_output = os.path.join(
        RESULTS_DIR,
        "09_disease_positive_count_matrix.tsv",
    )

    positive_count_matrix.to_csv(
        positive_count_matrix_output,
        sep="\t",
    )

    positive_count_heatmap_output = os.path.join(
        RESULTS_DIR,
        "10_disease_positive_count_heatmap.png",
    )

    plot_disease_positive_count_heatmap(
        matrix=positive_count_matrix,
        annotation=positive_count_annotation,
        output_path=positive_count_heatmap_output,
    )

    print(f"[SAVE] {positive_count_matrix_output}")

    # --------------------------------------------------------
    # 7. Final POS/SUS/NEG heatmap
    # --------------------------------------------------------
    print(
        "\n[7/7] Creating final disease-status heatmap..."
    )

    final_status_matrix, final_status_annotation = (
        create_disease_final_status_matrix(
            sample_depth_summary
        )
    )

    final_status_matrix_output = os.path.join(
        RESULTS_DIR,
        "11_disease_final_status_matrix.tsv",
    )

    final_status_matrix.to_csv(
        final_status_matrix_output,
        sep="\t",
    )

    final_status_heatmap_output = os.path.join(
        RESULTS_DIR,
        "12_disease_final_status_heatmap.png",
    )

    plot_disease_final_status_heatmap(
        matrix=final_status_matrix,
        annotation=final_status_annotation,
        output_path=final_status_heatmap_output,
    )

    print(f"[SAVE] {final_status_matrix_output}")

    # --------------------------------------------------------
    # Console output
    # --------------------------------------------------------
    print_positive_summary(
        sample_depth_summary
    )

    print("\n" + "=" * 90)
    print("Analysis completed")
    print("=" * 90)
    print(f"Results directory:\n{RESULTS_DIR}")


if __name__ == "__main__":
    main()