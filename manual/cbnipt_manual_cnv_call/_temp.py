import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr

df = pd.read_csv('temp.csv', sep='\t')
#df['Pos_Ratio'] = df['CopiesPer20uLWell'] / ( df['AcceptedDroplets'])
df['Pos_Ratio'] = df['CopiesPer20uLWell'] /(df['% cfDNA'] /100)
# seaborn style
sns.set(style="whitegrid")

# figure 생성
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Target 별 plot
for ax, target in zip(axes, df["Target"].unique()):

    subset = df[df["Target"] == target]

    # Spearman correlation
    rho, pval = pearsonr(
        subset["% cfDNA"],
        subset["Pos_Ratio"]
    )

    # scatter + regression line
    sns.regplot(
        data=subset,
        x="% cfDNA",
        y="Pos_Ratio",
        ax=ax,
        scatter_kws={"s": 80}
    )

    ax.set_title(
        f"{target}\nPearson r = {rho:.2f}, p = {pval:.3f}"
    )

    ax.set_xlabel("% cfDNA")
    ax.set_ylabel("Pos_Ratio")

plt.tight_layout()
plt.savefig("correlation_plots.ratio.png")