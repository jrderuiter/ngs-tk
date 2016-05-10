from collections import OrderedDict

import numpy as np
import pandas as pd

import seaborn as sns

from itertools import chain

from matplotlib import patches as mpatches
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize, rgb2hex


def clustermap(data, *args,
               col_annotation=None, col_palette=None,
               row_annotation=None, row_palette=None,
               legend_position=(1.2, 1),
               legend_offset=0.2, **kwargs):

    # Color annotation if given.
    color_maps = []

    if col_annotation is not None:
        kwargs['col_colors'], col_color_map = \
            color_annotation(col_annotation, col_palette)
        color_maps.append(col_color_map)

    if row_annotation is not None:
        kwargs['row_colors'], row_color_map = \
            color_annotation(row_annotation, row_palette)
        color_maps.append(row_color_map)

    # Draw clustermap.
    g = sns.clustermap(data, *args, **kwargs)

    # Draw legends.
    draw_legends(color_maps, ax=g.ax_heatmap,
                 position=legend_position, offset=legend_offset)

    return g


def draw_legends(color_maps, ax, position, offset,
                 loc='upper right', **kwargs):
    # Flatten list of dicts to items.
    color_maps = chain.from_iterable(
        (cmap.items() for cmap in color_maps))

    # Draw legends.
    for i, ((name, color_map), not_last) in enumerate(lookahead(color_maps)):
        bbox_position = (position[0], position[1] - (offset * i))
        legend = draw_legend(
            color_map, ax=ax, name=name, loc=loc,
            bbox_to_anchor=bbox_position, **kwargs)

        if not_last:
            ax.add_artist(legend)


def draw_legend(color_map, ax, name=None, **kwargs):
    patches = [mpatches.Patch(color=color, label=label)
               for label, color in color_map.items()]
    legend = ax.legend(handles=patches,
                       title=name, **kwargs)
    return legend


def color_annotation(df, colors, bg_color='#ffffff'):
    """Converts a data frame of annotations to colors."""

    # TODO: Fix boolean/nan case.

    colored_cols = OrderedDict()
    color_maps = OrderedDict()

    for (col_name, values), color in zip(df.items(), colors):
        colored, color_map = _color_column(values, color, bg_color)
        colored_cols[col_name] = colored
        color_maps[col_name] = color_map

    return pd.DataFrame(colored_cols), color_maps


def _color_column(series, color, bg_color):
    """Converts an annotation column to colors."""
    if series.dtype == bool:
        # Boolean case.
        return _color_bool_column(series, color, bg_color)
    elif series.dtype.kind in 'biufc':
        # Numeric case.
        return _color_numeric_column(series, color, bg_color)
    elif series.dtype.kind == 'O':
        # Strings and categorical case..
        return _color_cat_column(series, color, bg_color)
    else:
        raise ValueError('Unsupported dtype: {}'.format(series.dtype))


def _color_bool_column(series, color, bg_color):
    """Converts a boolean annotation column to colors."""
    return _color_numeric_column(series.astype(int), color, bg_color)


def _color_numeric_column(series, color, bg_color, n=200):
    """Converts a numeric annotation column to colors."""

    cmap = LinearSegmentedColormap.from_list(
        name='custom_cmap', colors=[bg_color, color], N=n)
    norm = Normalize(vmin=series.min(), vmax=series.max())
    mappable = ScalarMappable(norm=norm, cmap=cmap)

    rgba_colors = mappable.to_rgba(series.values)
    hex_colors = (rgb2hex(rgba) for rgba in rgba_colors)

    mapped = pd.Series(hex_colors, index=series.index, name=series.name)

    return mapped, color


def _color_cat_column(series, colors, bg_color):
    """Converts a categorical annotation column to colors."""

    if series.dtype.name == 'category':
        str_values = series.cat.categories
    else:
        str_values = set(series.dropna())

    if isinstance(colors, list):
        colors = [_rgb_to_hex(c, normalized=True) for c in colors]
        color_map = dict(zip(str_values, colors))
        mapped = series.map(color_map).fillna(bg_color)
    else:
        numeric_values = np.linspace(0, 1, len(str_values) + 1)[1:]
        numeric_map = dict(zip(str_values, numeric_values))

        mapped = _color_numeric_column(series.map(numeric_map),
                                       colors, bg_color=bg_color)[0]
        color_map = dict(zip(series.values, mapped.values))

    return mapped, color_map


def _rgb_to_hex(rgb, normalized=True):
    if normalized:
        rgb = tuple(map(lambda x: int(x * 255), rgb))
    return '#%02x%02x%02x' % rgb


def lookahead(iterable):
    """Pass through all values from the given iterable, augmented by the
    information if there are more values to come after the current one
    (True), or if it is the last value (False).
    """
    # Get an iterator and pull the first value.
    it = iter(iterable)
    last = next(it)
    # Run the iterator to exhaustion (starting from the second value).
    for val in it:
        # Report the *previous* value (more to come).
        yield last, True
        last = val
    # Report the last value.
    yield last, False
