import numpy as np
from matplotlib import pyplot as plt


def plot_cnv(data, y, reference,
             chrom='chrom', position='position',
             ax=None, color=None):
    """Plots a CNV profile for a single sample."""

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

    # Plot data points for each chromosome individually.
    groups = data.groupby(chrom)
    for chrom_id in chrom_ids:
        grp = groups.get_group(chrom_id)
        grp = grp.sort_values(position, ascending=True)

        x = grp[position] + chrom_offset[chrom_id]
        ax.plot(x, grp[y], '.', color=color)

    # Draw dividers and xticklabels.
    for loc in chrom_cumsums[1:-1]:
        ax.axvline(loc, color='grey', lw=0.5, zorder=5)

    ax.set_xticks((chrom_cumsums[:-1] + chrom_cumsums[1:]) / 2)
    ax.set_xticklabels(chrom_ids)

    # Add axis labels and title.
    ax.set_xlabel(chrom)
    ax.set_ylabel(y)

    ax.set_title(y)

    return ax


def plot_segments(data, y, reference, chrom='chrom',
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
