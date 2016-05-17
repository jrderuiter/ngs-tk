from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

# noinspection PyUnresolvedReferences
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

import numpy as np
import seaborn as sns
from natsort import natsorted
from matplotlib import pyplot as plt
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage

from .resample import resample

from ngs_tk.plotting.genomic import plot_genomic, plot_genomic_segments


def plot_profile(data, **kwargs):
    """Plots a CNV profile for a single sample."""

    # Plot data points.
    ax = plot_genomic(data, **kwargs)

    # Draw line for at y = 0.
    ax.axhline(0, color='grey', lw=0.5, zorder=-1)

    return ax


def plot_segments(data, **kwargs):
    """ Plots CNV segments for a single sample."""

    # Plot data points.
    ax = plot_genomic_segments(data, **kwargs)

    # Draw line for at y = 0.
    ax.axhline(0, color='grey', lw=0.5, zorder=-1)

    return ax


def plot_heatmap(data, reference=None, chromosomes=None,
                 metric='correlation', method='complete',
                 plot_bin_size=None, **kwargs):
    """Plots a (clustered) CNV heatmap for multiple samples."""

    if chromosomes is None:
        chromosomes = natsorted(data.index.levels[0])

    # Sort data by position.
    data = data.sort_index()
    data = data.ix[chromosomes]

    # Drop infinite values and NAs.
    data = data.replace([np.inf, -np.inf], np.nan)
    data = data.dropna()

    # Do clustering.
    dist = pdist(data.T, metric=metric)
    z = linkage(dist, method=method)

    # Resample data if needed for plotting.
    if plot_bin_size is not None:
        if reference is None:
            raise ValueError('Reference is required for resampling')
        data = resample(data, reference=reference, bin_size=plot_bin_size)

    # Plot heatmap.
    g = sns.clustermap(data.T, linewidths=0, col_cluster=False,
                       row_linkage=z, **kwargs)
    g.ax_heatmap.set_xticks([])

    # Plot chromosome breaks.
    breaks = data.groupby(level=0).size()[chromosomes].cumsum()
    for loc in breaks[:-1]:
        g.ax_heatmap.axvline(loc, color='black', lw=1)

    # Add chromosome labels.
    label_pos = np.concatenate([[0], breaks])
    label_pos = (label_pos[:-1] + label_pos[1:]) / 2

    g.ax_heatmap.set_xticks(label_pos)
    g.ax_heatmap.set_xticklabels(chromosomes, rotation=0)

    # Rotate sample labels.
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)

    # Label axes.
    g.ax_heatmap.set_xlabel('Chromosome')

    return g
