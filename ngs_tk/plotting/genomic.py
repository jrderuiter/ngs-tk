
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from natsort import natsorted


def plot_genomic(data, reference, chrom=None, position=None, y=None, hue=None,
                 ax=None, color=None, chromosomes=None, ymin=None,ymax=None, **kwargs):
    """Plots a genomic data profile for a single sample."""

    # Pre-process data.
    data, chrom, position, y = _preprocess_data(data, chrom, position, y)

    # Default axes.
    if ax is None:
        _, ax = plt.subplots()

    # Lookup chromosomes.
    if chromosomes is None:
        chromosomes = natsorted(data.index.levels[0])
    else:
        chromosomes = list(chromosomes)

    # Lookup chromosome lengths/offsets.
    chrom_lengths = [len(reference[k]) for k in chromosomes]
    chrom_cumsums = np.concatenate([[0], np.cumsum(chrom_lengths)])

    chrom_offset = dict(zip(chromosomes, chrom_cumsums))

    # Build frame with plotting data.
    plot_data = pd.DataFrame({
        'chrom': pd.Categorical(data[chrom], chromosomes),
        'x': data[position] + data[chrom].map(chrom_offset),
        'y': data[y]})

    plot_data.dropna(subset=['chrom'], inplace=True)

    # Do the actual plotting.
    if hue is not None:
        # Plot by hue group.
        plot_data['hue'] = data[hue]

        for label, grp in plot_data.groupby('hue'):
            ax.plot(grp['x'], grp['y'], '.', label=label, **kwargs)

        ax.legend(frameon=True, title=hue)
    else:
        # Plot by chromosome.
        for _, grp in plot_data.groupby('chrom'):
            ax.plot(grp['x'], grp['y'], '.', color=color, **kwargs)

    # Draw dividers and x-tick-labels.
    for loc in chrom_cumsums[1:-1]:
        ax.axvline(loc, color='grey', lw=0.5, zorder=5)

    ax.set_xticks((chrom_cumsums[:-1] + chrom_cumsums[1:]) / 2)
    ax.set_xticklabels(chromosomes)

    # Add axis labels and title.
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('Value')

    ax.set_title(y)

    # Set limits.
    ax.set_ylim(ymin, ymax)

    return ax


def _preprocess_data(data, chrom=None, position=None, y=None):
    """Pre-processes Series/Frame inputs into common frame format."""

    # Get chrom/position from index if not given.
    if chrom is None:
        chrom = data.index.names[0]

    if position is None:
        position = data.index.names[1]

    # Get y if not given.
    if isinstance(data, pd.Series):
        y = data.name
    else:
        y = data.columns[0]

    # Convert to data frame.
    data = data.reset_index()

    return data, chrom, position, y
