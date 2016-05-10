
import seaborn as sns
from matplotlib import pyplot as plt


def plot_set_expression(expr, set_genes, sign_genes=None, cmap_kwargs=None):
    """Plots clustermap of expression of (significant) genes in gene set."""

    # Subset to significant genes if given.
    if sign_genes is not None:
        set_genes = sign_genes & set_genes

    # Get expression for gene set.
    set_expr = expr.ix[set_genes].dropna()

    # Draw clustermap.
    g = sns.clustermap(set_expr, **cmap_kwargs)
    g.ax_heatmap.set_xticks([])

    # Correct label orientation (bug in Seaborn).
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)

    return g
