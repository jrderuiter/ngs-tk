from collections import OrderedDict

import numpy as np
import pandas as pd

from matplotlib import patches as mpatches
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize, rgb2hex


def color_annotation(df, colors, bg_color='#ffffff'):
    """Converts a data frame of annotations to colors."""

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

    if series.dtype.name == 'category' and series.cat.ordered:
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
                                       colors, bg_color=bg_color)
        color_map = dict(zip(series, mapped.values))

    return mapped, color_map


def _rgb_to_hex(rgb, normalized=True):
    if normalized:
        rgb = tuple(map(lambda x: int(x * 255), rgb))
    return '#%02x%02x%02x' % rgb


def draw_legends(color_maps, ax, **kwargs):
    for name, color_map in color_maps.items():
        draw_legend(color_map, ax=ax, name=name, **kwargs)


def draw_legend(color_map, ax, name, **kwargs):
    patches = [mpatches.Patch(color=color, label=label)
               for label, color in color_map.items()]
    legend = ax.legend(handles=patches, frameon=True,
                       title=name, **kwargs)
    return legend
