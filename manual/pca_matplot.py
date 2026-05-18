# pca_plot.py

# ==============================================================================
# 1. COMPATIBILITY PATCH (Must remain at the top)
# ==============================================================================
import inspect
import collections

# Patch for Python 3.11+ compatibility with older Keras/UMAP libraries
if not hasattr(inspect, 'ArgSpec'):
    inspect.ArgSpec = collections.namedtuple('ArgSpec', ['args', 'varargs', 'keywords', 'defaults'])

# ==============================================================================
# 2. IMPORTS
# ==============================================================================
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
import umap.umap_ as umap

# [MODIFIED] Removed Plotly imports. Kept Matplotlib and Seaborn for static plotting.
import matplotlib
matplotlib.use('Agg') # Essential to prevent display errors in GUI-less server environments
import matplotlib.pyplot as plt
import seaborn as sns

# ==============================================================================
# 3. DIMENSIONALITY REDUCTION FUNCTIONS
# ==============================================================================

def compute_pca(expression_matrix: pd.DataFrame, n_components: int = 2) -> pd.DataFrame:
    """
    Computes Principal Component Analysis (PCA) on the expression data.
    """
    # Standardize the data (mean=0, variance=1)
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(expression_matrix)
    
    # Run PCA
    pca = PCA(n_components=n_components, random_state=42)
    pca_results = pca.fit_transform(scaled_data)
    
    # Format as DataFrame
    columns = [f'PC{i+1}' for i in range(n_components)]
    return pd.DataFrame(pca_results, columns=columns, index=expression_matrix.index)

def compute_tsne(expression_matrix: pd.DataFrame, n_components: int = 2) -> pd.DataFrame:
    """
    Computes t-Distributed Stochastic Neighbor Embedding (t-SNE).
    """
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(expression_matrix)
    
    # init='pca' is highly recommended for data with very high feature counts (like RNA-seq)
    tsne = TSNE(
        n_components=n_components, 
        init='pca', 
        learning_rate='auto', 
        random_state=42
    )
    tsne_results = tsne.fit_transform(scaled_data)
    
    columns = [f'tSNE{i+1}' for i in range(n_components)]
    return pd.DataFrame(tsne_results, columns=columns, index=expression_matrix.index)

def compute_umap(expression_matrix: pd.DataFrame, n_components: int = 2) -> pd.DataFrame:
    """
    Computes Uniform Manifold Approximation and Projection (UMAP).
    """
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(expression_matrix)
    
    # Run UMAP
    reducer = umap.UMAP(n_components=n_components, random_state=42)
    umap_results = reducer.fit_transform(scaled_data)
    
    columns = [f'UMAP{i+1}' for i in range(n_components)]
    return pd.DataFrame(umap_results, columns=columns, index=expression_matrix.index)


# ==============================================================================
# 4. VISUALIZATION FUNCTIONS (MATPLOTLIB)
# ==============================================================================

# [MODIFIED] Completely refactored to return a Matplotlib Figure using Seaborn for mapping
def _create_scatter_plot(df: pd.DataFrame, x_col: str, y_col: str, color_col: str, symbol_col: str, title: str) -> plt.Figure:
    """
    Internal helper function to maintain consistent styling across all Matplotlib figures.
    """
    fig, ax = plt.subplots(figsize=(10, 8))
    
    sns.scatterplot(
        data=df, 
        x=x_col, 
        y=y_col, 
        hue=color_col,     # Maps the 'class' to the point color
        style=symbol_col,  # Maps the 'group' (batch) to the point shape
        s=100,             # Marker size
        alpha=0.8,
        palette='Set1',    # A distinct color palette
        ax=ax
    )
    
    ax.set_title(title, fontsize=16, fontweight='bold')
    ax.set_xlabel(x_col, fontsize=14)
    ax.set_ylabel(y_col, fontsize=14)
    
    # Move the legend outside the plot area to prevent it from covering data points
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., title=f"{color_col} & {symbol_col}")
    
    # Add a subtle grid
    ax.grid(True, linestyle='--', alpha=0.6)
    
    # Adjust layout to ensure the external legend is not cut off
    fig.tight_layout()
    
    return fig

# [MODIFIED] Changed return type to plt.Figure
def plot_pca(df: pd.DataFrame, class_col: str = 'class', group_col: str = 'group') -> plt.Figure:
    """
    Generates a static PCA scatter plot to visualize batch effects and classes.
    Requires 'PC1' and 'PC2' columns in the DataFrame.
    """
    return _create_scatter_plot(df, x_col='PC1', y_col='PC2', color_col=class_col, symbol_col=group_col, title='PCA: RNA-seq Expression')

# [MODIFIED] Changed return type to plt.Figure
def plot_tsne(df: pd.DataFrame, class_col: str = 'class', group_col: str = 'group') -> plt.Figure:
    """
    Generates a static t-SNE scatter plot to visualize batch effects and classes.
    Requires 'tSNE1' and 'tSNE2' columns in the DataFrame.
    """
    return _create_scatter_plot(df, x_col='tSNE1', y_col='tSNE2', color_col=class_col, symbol_col=group_col, title='t-SNE: RNA-seq Expression')

# [MODIFIED] Changed return type to plt.Figure
def plot_umap(df: pd.DataFrame, class_col: str = 'class', group_col: str = 'group') -> plt.Figure:
    """
    Generates a static UMAP scatter plot to visualize batch effects and classes.
    Requires 'UMAP1' and 'UMAP2' columns in the DataFrame.
    """
    return _create_scatter_plot(df, x_col='UMAP1', y_col='UMAP2', color_col=class_col, symbol_col=group_col, title='UMAP: RNA-seq Expression')