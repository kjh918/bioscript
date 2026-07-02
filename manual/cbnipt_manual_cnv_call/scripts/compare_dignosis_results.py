import os
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


def aggregate_diagnosis(diagnoses):
    """
    Aggregates the diagnosis for a single syndrome based on its features.
    - POSITIVE (2): If ALL features are POSITIVE.
    - SUSPICIOUS (1): If not all are POSITIVE, but at least one is POSITIVE or SUSPICIOUS.
    - NEGATIVE (0): If ALL features are NEGATIVE.
    """
    unique_vals = set(diagnoses.dropna().str.upper())
    
    if len(unique_vals) == 1 and "POSITIVE" in unique_vals:
        return 2  # All POSITIVE
    elif "POSITIVE" in unique_vals or "SUSPICIOUS" in unique_vals:
        return 1  # Next level down (Mixed / Suspicious)
    else:
        return 0  # NEGATIVE


def load_and_summarize_samples(data_dir):
    """
    Loads all TSV files in the directory and summarizes the syndrome status per sample.
    """
    print(os.path.join(data_dir, "*24_04*/data/*.clinical_report.tsv"))
    tsv_files = glob.glob(os.path.join(data_dir, "*24_04*/data/*.clinical_report.tsv"), recursive=True)
    if not tsv_files:
        raise FileNotFoundError(f"No TSV files found in {data_dir}")

    all_summaries = []

    for file_path in tsv_files:
        # Extract sample ID from the filename (e.g., sample_A.tsv -> sample_A)
        sample_id = os.path.basename(file_path).replace(".clinical_report.tsv", "")
        
        df = pd.read_csv(file_path, sep="\t")
        
        # Ensure required columns exist
        if not all(col in df.columns for col in ["SYNDROME", "DIAGNOSIS"]):
            print(f"Skipping {file_path}: Missing required columns.")
            continue
            
        # Group by syndrome and apply the aggregation logic
        summary = df.groupby("SYNDROME")["DIAGNOSIS"].apply(aggregate_diagnosis).reset_index()
        summary["SAMPLE_ID"] = sample_id
        
        all_summaries.append(summary)
        
    # Combine all sample summaries into a single DataFrame
    combined_df = pd.concat(all_summaries, ignore_index=True)
    
    # Pivot to create a matrix: Rows = Samples, Columns = Syndromes
    heatmap_data = combined_df.pivot(index="SAMPLE_ID", columns="SYNDROME", values="DIAGNOSIS")
    
    # Fill any missing syndromes for specific samples with NEGATIVE (0)
    heatmap_data = heatmap_data.fillna(0)
    
    return heatmap_data


def plot_syndrome_heatmap(heatmap_data, output_path):
    """
    Generates and saves a heatmap from the sample-syndrome matrix.
    """
    plt.figure(figsize=(14, 8))
    
    # Define custom colormap: 0=Blue (Negative), 1=Yellow (Suspicious), 2=Red (Positive)
    cmap = mcolors.ListedColormap(["#f7fbff", "#ffeda0", "#f03b20"])
    bounds = [-0.5, 0.5, 1.5, 2.5]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)
    
    ax = sns.heatmap(
        heatmap_data, 
        cmap=cmap, 
        norm=norm,
        linewidths=0.5, 
        linecolor="gray",
        cbar_kws={"ticks": [0, 1, 2]}
    )
    
    # Configure colorbar labels
    cbar = ax.collections[0].colorbar
    cbar.set_ticklabels(["Negative", "Suspicious", "Positive"])
    
    plt.title("NIPT Syndrome Detection Summary Across Samples", fontsize=16, pad=20)
    plt.xlabel("Syndrome", fontsize=12)
    plt.ylabel("Sample ID", fontsize=12)
    
    # Rotate x-axis labels for readability
    plt.xticks(rotation=45, ha="right", fontsize=10)
    plt.yticks(rotation=0, fontsize=10)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    
    print(f"Heatmap successfully saved to {output_path}")


def main():
    # Define paths based on the setup script structure
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    base_dir = '/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Results/temp'
    
    output_heatmap = os.path.join(base_dir, "results", "syndrome_heatmap.png")
    output_matrix = os.path.join(base_dir, "results", "syndrome_matrix.tsv")
    
    os.makedirs(os.path.dirname(output_heatmap), exist_ok=True)
    os.makedirs(os.path.dirname(output_matrix), exist_ok=True   )
    # Process data and generate matrix
    print("Aggregating sample data...")
    heatmap_data = load_and_summarize_samples(base_dir)
    
    heatmap_data.to_csv(output_matrix, sep="\t")
    print(f"Matrix data saved to {output_matrix}")
    
    # Plot and save the heatmap
    print("Generating heatmap...")
    plot_syndrome_heatmap(heatmap_data, output_heatmap)


if __name__ == "__main__":
    main()