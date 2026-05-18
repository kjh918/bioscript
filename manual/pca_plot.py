# pca_plot.py

# ==============================================================================
# 1. COMPATIBILITY PATCH (Must remain at the top)
# ==============================================================================
import inspect
import collections
# ==============================================================================
# 2. IMPORTS
# ==============================================================================
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
import umap.umap_ as umap
import plotly.express as px
import plotly.graph_objects as go

import matplotlib
matplotlib.use('Agg') # 필수: 리눅스 서버 등 GUI가 없는 환경에서 에러 방지
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
# 4. VISUALIZATION FUNCTIONS
# ==============================================================================

def _create_scatter_plot(df: pd.DataFrame, x_col: str, y_col: str, color_col: str, symbol_col: str, title: str) -> go.Figure:
    """
    Internal helper function to maintain consistent styling across all Plotly figures.
    """
    fig = px.scatter(
        df, 
        x=x_col, 
        y=y_col, 
        color=color_col,     # Maps the 'class' to the point color
        symbol=symbol_col,   # Maps the 'group' (batch) to the point shape
        title=title,
        hover_data=df.columns, 
        opacity=0.8
    )
    
    fig.update_layout(
        template="plotly_white",
        title_font_size=20,
        legend_title_text=f"{color_col.capitalize()} & {symbol_col.capitalize()}", 
        legend=dict(itemsizing='constant')
    )
    
    # Increased marker size to make shapes (symbols) easier to distinguish
    fig.update_traces(marker=dict(size=8, line=dict(width=0.5, color='DarkSlateGrey')))
    return fig


def plot_pca(df: pd.DataFrame, class_col: str = 'class', group_col: str = 'group') -> go.Figure:
    """
    Generates a PCA scatter plot to visualize batch effects and classes.
    Requires 'PC1' and 'PC2' columns in the DataFrame.
    """
    return _create_scatter_plot(df, x_col='PC1', y_col='PC2', color_col=class_col, symbol_col=group_col, title='PCA: RNA-seq Expression')


def plot_tsne(df: pd.DataFrame, class_col: str = 'class', group_col: str = 'group') -> go.Figure:
    """
    Generates a t-SNE scatter plot to visualize batch effects and classes.
    Requires 'tSNE1' and 'tSNE2' columns in the DataFrame.
    """
    return _create_scatter_plot(df, x_col='tSNE1', y_col='tSNE2', color_col=class_col, symbol_col=group_col, title='t-SNE: RNA-seq Expression')


def plot_umap(df: pd.DataFrame, class_col: str = 'class', group_col: str = 'group') -> go.Figure:
    """
    Generates a UMAP scatter plot to visualize batch effects and classes.
    Requires 'UMAP1' and 'UMAP2' columns in the DataFrame.
    """
    return _create_scatter_plot(df, x_col='UMAP1', y_col='UMAP2', color_col=class_col, symbol_col=group_col, title='UMAP: RNA-seq Expression')

