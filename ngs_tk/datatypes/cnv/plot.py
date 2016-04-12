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


def plot_profile(data, reference, ax=None, color=None, chromosomes=None):
    """Plots a CNV profile for a single sample."""

    if ax is None:
        _, ax = plt.subplots()

    # Lookup chromosomes.
    if chromosomes is None:
        chromosomes = natsorted(data.index.levels[0])

    # Lookup chromosome lengths/offsets.
    chrom_lengths = [len(reference[k]) for k in chromosomes]
    chrom_cumsums = np.concatenate([[0], np.cumsum(chrom_lengths)])

    chrom_offset = dict(zip(chromosomes, chrom_cumsums))

    # Plot data points for each chromosome individually.
    groups = data.groupby(level=0)
    for chrom in chromosomes:
        grp = groups.get_group(chrom)
        grp = grp.sort_index()

        x = grp.index.get_level_values(1) + chrom_offset[chrom]
        ax.plot(x, grp.values, '.', color=color)

    # Draw dividers and x-tick-labels.
    for loc in chrom_cumsums[1:-1]:
        ax.axvline(loc, color='grey', lw=0.5, zorder=5)

    ax.set_xticks((chrom_cumsums[:-1] + chrom_cumsums[1:]) / 2)
    ax.set_xticklabels(chromosomes)

    # Add axis labels and title.
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('Value')

    # ax.set_title(y)

    return ax


def plot_profile_segments(data, y, reference, chrom='chrom',
                          start='start', end='end', ax=None):
    """ Plots segments on a drawn CNV axis. """

    if ax is None:
        _, ax = plt.subplots()

    # Lookup chromosomes.
    try:
        # Try to get order from categorical.
        chrom_ids = data[chrom].cat.categories
    except AttributeError:
        # If that fails, just fall back on unique.
        chrom_ids = data[chrom].unique()

    # Lookup chromosome lengths/offsets.
    chrom_lengths = [len(reference[k]) for k in chrom_ids]
    chrom_cumsums = np.concatenate([[0], np.cumsum(chrom_lengths)])

    chrom_offset = dict(zip(chrom_ids, chrom_cumsums))

    for _, row in data.iterrows():
        seg_start = row[start] + chrom_offset[row[chrom]]
        seg_end = row[end] + chrom_offset[row[chrom]]
        ax.plot([seg_start, seg_end], [row[y]] * 2, color='red')


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
