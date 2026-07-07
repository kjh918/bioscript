import os
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


def map_diagnosis(val):
    """
    텍스트 기반 진단 결과를 Heatmap 컬러 매핑을 위한 숫자 스코어로 변환합니다.
    """
    val = str(val).upper()
    if "POSITIVE" in val: return 2
    if "SUSPICIOUS" in val: return 1
    if "NEGATIVE" in val: return 0
    return 0  # Default fallback


def load_detailed_samples(data_dir):
    """
    모든 샘플의 clinical_report.tsv를 읽어, 증후군의 '세부 마커(Feature)' 단위까지
    펼쳐진 매트릭스(Sample x Feature)를 생성합니다.
    """
    search_pattern = os.path.join(data_dir, "*24_04*/data/*.clinical_report.tsv")
    print(f"Searching for files: {search_pattern}")
    
    tsv_files = glob.glob(search_pattern, recursive=True)
    if not tsv_files:
        raise FileNotFoundError(f"No TSV files found using pattern: {search_pattern}")

    all_data = []

    for file_path in tsv_files:
        sample_id = os.path.basename(file_path).replace(".clinical_report.tsv", "")
        
        df = pd.read_csv(file_path, sep="\t")
        
        # 필수 컬럼 존재 확인
        req_cols = ["SYNDROME", "FEATURE_NAME", "FEATURE_TYPE", "FEAT_RANK", "DIAGNOSIS"]
        if not all(col in df.columns for col in req_cols):
            print(f"Skipping {file_path}: Missing required columns.")
            continue
            
        # 성별 판독 행은 질환(Positive/Negative)이 아니므로 Heatmap에서 제외
        df = df[df["SYNDROME"] != "NIPT_SEX"].copy()
        
        df["SAMPLE_ID"] = sample_id
        df["SCORE"] = df["DIAGNOSIS"].apply(map_diagnosis)
        
        # X축에 표시될 이름 포맷팅: [증후군] 피처명 (피처타입)
        # 예: [DiGeorgy syndrome] TBX1 (CoreGene)
        df["DISPLAY_NAME"] = df.apply(
            lambda r: f"[{r['SYNDROME']}]\n{r['FEATURE_NAME']} ({r['FEATURE_TYPE']})", axis=1
        )
        
        all_data.append(df)
        
    if not all_data:
        raise ValueError("No valid data parsed from TSV files.")

    combined_df = pd.concat(all_data, ignore_index=True)
    
    # [핵심] 컬럼(X축)의 논리적 정렬 순서 보장
    # 1순위: 증후군명, 2순위: 피처의 계층적 랭크(염색체->영역->유전자)
    col_order_df = combined_df[["DISPLAY_NAME", "SYNDROME", "FEAT_RANK"]].drop_duplicates()
    col_order_df = col_order_df.sort_values(["SYNDROME", "FEAT_RANK"])
    ordered_cols = col_order_df["DISPLAY_NAME"].tolist()
    
    # [MODIFIED] 변경 이유: 중복된 마커 이름이 들어올 경우 발생하는 Index Error 방지.
    # pivot 대신 pivot_table을 사용하고, 중복 시 가장 심각한 상태(max)를 반영 (2: POS > 1: SUS > 0: NEG)
    heatmap_data = combined_df.pivot_table(
        index="SAMPLE_ID", 
        columns="DISPLAY_NAME", 
        values="SCORE", 
        aggfunc="max"
    )
    
    # 데이터가 없는 셀은 0(Negative)으로 채우고, 미리 정해둔 순서대로 컬럼 재배치
    heatmap_data = heatmap_data.fillna(0)[ordered_cols]
    
    return heatmap_data


def plot_syndrome_heatmap(heatmap_data, output_path):
    """
    세분화된 Sample-Feature 매트릭스로부터 고해상도 Heatmap을 생성합니다.
    """
    # 컬럼(세부 마커)의 개수가 많으므로 동적으로 가로 길이를 늘림
    num_features = heatmap_data.shape[1]
    num_samples = heatmap_data.shape[0]
    
    fig_width = max(16, num_features * 0.4)
    fig_height = max(8, num_samples * 0.5 + 4)
    
    plt.figure(figsize=(fig_width, fig_height))
    
    # 커스텀 컬러맵: 0=Blue(Negative), 1=Yellow(Suspicious), 2=Red(Positive)
    cmap = mcolors.ListedColormap(["#f7fbff", "#ffeda0", "#f03b20"])
    bounds = [-0.5, 0.5, 1.5, 2.5]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)
    
    ax = sns.heatmap(
        heatmap_data, 
        cmap=cmap, 
        norm=norm,
        linewidths=0.5, 
        linecolor="lightgray",
        cbar_kws={"ticks": [0, 1, 2], "shrink": 0.5}
    )
    
    # 컬러바 설정
    cbar = ax.collections[0].colorbar
    cbar.set_ticklabels(["Negative", "Suspicious", "Positive"])
    
    plt.title("NIPT Detailed Marker Detection Summary Across Samples", fontsize=18, pad=20, fontweight="bold")
    plt.xlabel("Target Disease Features (Chromosome → Region → Gene)", fontsize=14, fontweight="bold")
    plt.ylabel("Sample ID", fontsize=14, fontweight="bold")
    
    # X축 라벨이 길고 많으므로 90도 회전
    plt.xticks(rotation=90, ha="center", fontsize=9)
    plt.yticks(rotation=0, fontsize=11)
    
    # 증후군(Syndrome)이 바뀌는 지점에 세로선을 그어 시각적 분리감 부여
    current_syndrome = ""
    for i, col_name in enumerate(heatmap_data.columns):
        syn = col_name.split("]")[0].strip("[")
        if syn != current_syndrome:
            if i != 0:
                ax.axvline(i, color="black", linewidth=1.5, linestyle="-", zorder=3)
            current_syndrome = syn

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    
    print(f"Heatmap successfully saved to {output_path}")


def main():
    # 작업 경로 설정
    base_dir = '/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Results/temp'
    
    output_heatmap = os.path.join(base_dir, "results", "syndrome_detailed_heatmap.png")
    output_matrix = os.path.join(base_dir, "results", "syndrome_detailed_matrix.tsv")
    
    os.makedirs(os.path.dirname(output_heatmap), exist_ok=True)
    os.makedirs(os.path.dirname(output_matrix), exist_ok=True)
    
    # 데이터 처리 및 매트릭스 생성
    print("Aggregating detailed sample data...")
    heatmap_data = load_detailed_samples(base_dir)
    
    heatmap_data.to_csv(output_matrix, sep="\t")
    print(f"Matrix data saved to {output_matrix}")
    
    # Heatmap 플로팅
    print("Generating heatmap...")
    plot_syndrome_heatmap(heatmap_data, output_heatmap)


if __name__ == "__main__":
    main()