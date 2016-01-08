from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

# noinspection PyUnresolvedReferences
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)


from collections import OrderedDict

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection


ONCOPRINT_ALTERATIONS = OrderedDict(
    [('cna', ['homdel', 'amp']),
     ('mutation', ['trunc', 'frameshift', 'missense', 'inframe'])])

ONCOPRINT_COLORS = {
    'cna': {
        'homdel': 'blue',
        'amp': 'red'
    },
    'mutation': {
        'trunc': 'black',
        'frameshift': 'black',
        'inframe': 'peru',
        'missense': 'green'
    }
}

ONCOPRINT_TYPES = {'cna': 'fill',
                   'mutation': 'inset'}


def oncoprint(data, **kwargs):
    """Creates oncoplot in cBioPortals oncoprint style."""
    return oncoplot(data, alterations=ONCOPRINT_ALTERATIONS,
                    colors=ONCOPRINT_COLORS,
                    plot_types=ONCOPRINT_TYPES, **kwargs)


def oncoplot(data, alterations, colors, plot_types, fig_kws=None):
    """Creates oncoprint style plots of alterations."""

    type_funcs = {
        'fill': _fill,
        'inset': _inset_square
    }

    # Fill default args.
    fig_kws = fig_kws or {}

    # Sort genes/samples before plotting.
    data = sort_alterations(data, alterations)

    # Create figure and plot empty background.
    fig, ax = plt.subplots(**fig_kws)
    _plot_background(data, ax)

    # Plot per alteration type.
    for type_, grp in data.groupby('type'):
        type_colors = colors[type_]
        plot_type = plot_types[type_]

        plot_func = type_funcs[plot_type]
        plot_func(grp, type_colors, ax=ax)

    return fig, ax


def _plot_background(data, ax):
    # Plot empty 'null' background.
    null_data = pd.DataFrame(0, index=data['gene'].cat.categories,
                             columns=data['sample'].cat.categories)
    sns.heatmap(null_data, cmap=discrete_cmap(['lightgrey']),
                ax=ax, linewidths=0.2, cbar=False)


def _fill(data, colors, ax):
    if len(data) == 0:
        return

    data = data.copy()

    # Subset data and colors.
    data = data.ix[data['alteration'].isin(colors)]
    colors = {k: v for k, v in colors.items()
              if k in set(data['alteration'])}

    # Build number/color maps.
    num_map, colors_ord = {}, []
    for i, (type_, color) in enumerate(colors.items()):
        num_map[type_] = i
        colors_ord.append(color)

    # Map to numeric and pivot.
    data['value'] = data['alteration'].map(num_map)
    mat = _pivot(data, values='value')

    # Mask na-values.
    mask = pd.isnull(mat)

    # Draw fills using heatmap.
    cmap = discrete_cmap(colors_ord)
    sns.heatmap(mat, cmap=cmap, ax=ax, mask=mask,
                cbar=False, linewidths=0.5)


def _inset_square(data, colors, ax):
    # Create lookups for gene/sample indices.
    sample_lookup = _cat_lookup(data['sample'])
    gene_lookup = _cat_lookup(data['gene'])

    # Mirror gene lookup
    n_genes = len(gene_lookup)
    gene_lookup = {k: (n_genes - 1) - v
                   for k, v in gene_lookup.items()}

    # Plot mutations.
    for type_, grp in data.groupby('alteration'):
        if len(grp) > 0:
            patches = (Rectangle(xy=(sample_lookup[row['sample']] + 0.2,
                                     gene_lookup[row['gene']] + 0.33),
                                 width=0.6, height=0.33)
                       for _, row in grp.iterrows())

            patches = PatchCollection(patches, facecolor=colors[type_])
            ax.add_collection(patches)


def sort_alterations(data, alterations, inplace=False):
    if not inplace:
        data = data.copy()

    # Order genes by alteration frequency.
    gene_ord = (data.groupby('gene')['sample'].nunique()
                .sort_values(ascending=False))

    data['gene'] = pd.Categorical(
        data['gene'], categories=list(gene_ord.index))

    # Split priorities for alterations.
    type_order = []
    alt_order = []
    for type_, type_alts in alterations.items():
        type_order.append(type_)
        alt_order += type_alts

    # Apply alteration type order.
    data['type'] = pd.Categorical(data['type'], type_order, ordered=True)
    data['alteration'] = pd.Categorical(data['alteration'],
                                        alt_order, ordered=True)

    # Now finally sort to obtain sample order.

    # We use the numeric codes of the categorical for sorting,
    # as sorting on categoricals after pivoting seems to sort on
    # the string values, not on the actual categorical.
    data['alt_num'] = data['alteration'].cat.codes
    data.ix[data['alt_num'] == -1, 'alt_num'] = np.nan

    # Here we pivot the data into matrix format, which we then sort
    # by its columns (in gene order) to get the desired sample ordering.
    mat = pd.pivot_table(data, index='sample', columns='gene',
                         values='alt_num', aggfunc=min)
    mat.sort_values(by=list(mat.columns), ascending=True, inplace=True)

    # Actually 're-order' samples in categorical.
    data['sample'] = _reorder_category(data['sample'], order=list(mat.index))

    return data


def _reorder_category(column, order):
    if column.dtype.name == 'category':
        if len(order) < len(column.cat.categories):
            order += [c for c in column.cat.categories
                      if c not in set(order)]
        new_column = column.cat.reorder_categories(order, ordered=True)
    else:
        new_column = pd.Categorical(column, categories=order, ordered=True)

    return new_column


def _pivot(data, values, aggfunc=min):
    return pd.pivot_table(data, columns='sample', index='gene',
                          values=values, aggfunc=aggfunc)


def _cat_lookup(series):
    categories = series.cat.categories
    return dict(zip(categories, range(len(categories))))


def discrete_cmap(colors):
    """Create an N-bin discrete colormap from specified colors."""

    if len(colors) == 0:
        return None
    else:
        if len(colors) == 1:
            # Repeat single colors as a work-around
            # for single color issue.
            colors = colors * 2

        n = len(colors)
        base = plt.cm.get_cmap()
        cmap_name = base.name + str(n)
        return base.from_list(cmap_name, colors, n)
