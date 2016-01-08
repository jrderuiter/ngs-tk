from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

# noinspection PyUnresolvedReferences
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)


import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt


def plot_profile(data, y, reference, chrom='chrom',
                 position='position', ax=None, color=None,
                 chromosomes=None):
    """Plots a CNV profile for a single sample."""

    if ax is None:
        _, ax = plt.subplots()

    # Lookup chromosomes.
    if chromosomes is None:
        try:
            # Try to get order from categorical.
            chrom_ids = data[chrom].cat.categories
        except AttributeError:
            # If that fails, just fall back on unique.
            chrom_ids = data[chrom].unique()
    else:
        chrom_ids = chromosomes

    # Lookup chromosome lengths/offsets.
    chrom_lengths = [len(reference[k]) for k in chrom_ids]
    chrom_cumsums = np.concatenate([[0], np.cumsum(chrom_lengths)])

    chrom_offset = dict(zip(chrom_ids, chrom_cumsums))

    # Plot data points for each chromosome individually.
    groups = data.groupby(chrom)
    for chrom_id in chrom_ids:
        grp = groups.get_group(chrom_id)
        grp = grp.sort_values(position, ascending=True)

        x = grp[position] + chrom_offset[chrom_id]
        ax.plot(x, grp[y], '.', color=color)

    # Draw dividers and x-tick-labels.
    for loc in chrom_cumsums[1:-1]:
        ax.axvline(loc, color='grey', lw=0.5, zorder=5)

    ax.set_xticks((chrom_cumsums[:-1] + chrom_cumsums[1:]) / 2)
    ax.set_xticklabels(chrom_ids)

    # Add axis labels and title.
    ax.set_xlabel(chrom)
    ax.set_ylabel(y)

    ax.set_title(y)

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


def plot_heatmap(data,  columns=None, chrom='chrom',
                 position='position', vline_color='black', **kwargs):
    """Plots a (clustered) CNV heatmap for multiple samples."""

    # Select all columns by default.
    if columns is None:
        columns = [c for c in data if c not in {chrom, position}]

    # Sort data by position.
    data = data.sort([chrom, position], ascending=True)

    # Plot heatmap.
    g = sns.clustermap(data[columns].T, linewidths=0,
                       col_cluster=False, **kwargs)
    g.ax_heatmap.set_xticks([])

    # Plot chromosome breaks.
    breaks = np.where(~data[chrom].duplicated(take_last=True))[0]
    breaks += 1

    for loc in breaks[:-1]:
        g.ax_heatmap.axvline(loc, color=vline_color)

    # Add chromosome labels.
    label_pos = np.concatenate([[0], breaks])
    label_pos = (label_pos[:-1] + label_pos[1:]) / 2

    g.ax_heatmap.set_xticks(label_pos)
    g.ax_heatmap.set_xticklabels(
        data[chrom].unique(), rotation=0)

    # Label axes.
    g.ax_heatmap.set_xlabel(chrom)

    return g
