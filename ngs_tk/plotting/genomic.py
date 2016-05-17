
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from natsort import natsorted


def plot_genomic(data, reference, chrom=None, position=None, y=None,
                 hue=None, ax=None, color=None, chromosomes=None,
                 ymin=None, ymax=None, **kwargs):
    """Plots a genomic data profile for a single sample."""

    # Default axes.
    if ax is None:
        _, ax = plt.subplots()

    # Pre-process data.
    data = _normalize_data(data, chrom, position, y)

    # Lookup chromosomes if not given.
    if chromosomes is None:
        chromosomes = natsorted(set(data['chrom']))
    else:
        chromosomes = list(chromosomes)

    # Lookup chromosome lengths/offsets.
    chrom_offsets = _chromosome_offsets(chromosomes, reference)

    # Build frame with plotting data.
    plot_data = pd.DataFrame({
        'chrom': pd.Categorical(data['chrom'], chromosomes),
        'x': data['position'] + data['chrom'].map(chrom_offsets),
        'y': data['y']})

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

    # Style axes.
    _draw_dividers(ax, chromosomes, chrom_offsets)
    _style_axes(ax, title=y, ymin=ymin, ymax=ymax)

    return ax


def _normalize_data(data, chrom=None, position=None, y=None):
    """Pre-processes Series/DataFrame inputs into common DataFrame format."""

    # Get chrom/position from index if not given.
    if chrom is None:
        chrom = data.index.names[0]

    if position is None:
        position = data.index.names[1]

    # Get y if not given.
    if y is None:
        if isinstance(data, pd.Series):
            y = data.name or 0
    else:
        y = data.columns[0]

    # Convert to data frame and rename columns.
    data = data.reset_index()
    data = data.rename(columns={chrom: 'chrom', position: 'position', y: 'y'})

    return data


def _chromosome_offsets(chromosomes, reference):
    chrom_lengths = [len(reference[k]) for k in chromosomes]
    chrom_cumsums = np.concatenate([[0], np.cumsum(chrom_lengths)])

    chrom_offset = dict(zip(chromosomes, chrom_cumsums))

    return chrom_offset


def _draw_dividers(ax, chromosomes, chrom_offsets):
    positions = np.sort(list(chrom_offsets.values()))

    # Draw dividers.
    for loc in positions[1:-1]:
        ax.axvline(loc, color='grey', lw=0.5, zorder=5)

    # Draw xtick labels.
    ax.set_xticks((positions[:-1] + positions[1:]) / 2)
    ax.set_xticklabels(chromosomes)

    return ax


def _style_axes(ax, title=None, ymin=None, ymax=None):
    # Add axis labels and title.
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('Value')

    if title is not None:
        ax.set_title(title)

    # Set limits.
    ax.set_ylim(ymin, ymax)

    return ax


def plot_genomic_segments(data, reference, y, chrom='chrom',
                          start='start', end='end', ax=None,
                          chromosomes=None, ymin=None, ymax=None,
                          style=True):
    # Default axes.
    if ax is None:
        _, ax = plt.subplots()

    # Lookup chromosomes.
    if chromosomes is None:
        chromosomes = natsorted(set(data[chrom]))
    else:
        chromosomes = list(chromosomes)

    # Lookup chromosome lengths/offsets.
    chrom_offsets = _chromosome_offsets(chromosomes, reference)

    # Plot segments.
    groups = data.groupby(chrom)
    for chrom_name in chromosomes:
        grp = groups.get_group(chrom_name)
        offset = chrom_offsets[chrom_name]

        for _, row in grp.iterrows():
            seg_start = row[start] + offset
            seg_end = row[end] + offset
            ax.plot([seg_start, seg_end], [row[y]] * 2, color='red')

    # Style axes.
    if style:
        _draw_dividers(ax, chromosomes, chrom_offsets)
        _style_axes(ax, title=y, ymin=ymin, ymax=ymax)

    return ax
