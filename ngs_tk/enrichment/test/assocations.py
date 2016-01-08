from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

import itertools

import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


ALT_MAP = {
    'mutex': 'less',
    'co-occ': 'greater',
    'two-sided': 'two-sided'
}


def convert_to_bool(frame):
    bool_cols = [_convert_column_to_bool(values)
                 for _, values in frame.items()]
    return pd.concat(bool_cols, axis=1)


def _convert_column_to_bool(series):
    if series.dtype == bool:
        return series
    elif series.dtype == object and type(series[0]) == str:
        return _convert_str_column_to_bool(series)
    else:
        raise ValueError('Unsupported dtype: {}'.format(series.dtype))


def _convert_str_column_to_bool(series):
    col_names = ('{}_{}'.format(series.name, v) for v in series.unique())
    col_values = list(series == v for v in series.unique())
    return pd.DataFrame(dict(zip(col_names, col_values)))


def contingency_table(data):
    if len(data.columns) != 2:
        raise ValueError('DataFrame does not have two columns')

    # Ensure data is boolean.
    data = convert_to_bool(data)

    # Create 2D contingency table.
    table = pd.pivot_table(data, index=data.columns[0],
                           columns=data.columns[1], aggfunc=len)
    table = table.fillna(0)

    return table


def test_association(data, test_type='two-sided'):
    # Try to look-up test type.
    try:
        alt = ALT_MAP[test_type]
    except KeyError:
        alt = test_type

    # Test using fishers_exact.
    table = contingency_table(data)
    p_value = fisher_exact(table, alternative=alt)[1]

    return p_value


def test_associations(data, test_types=('two-sided',), threshold=None,
                      corr_method='fdr_bh', associations=None):
    if associations is None:
        associations = itertools.combinations(data.columns, 2)

    row_gen = ((a, b, test_type,
                test_association(data[[a, b]], test_type=test_type))
               for a, b in associations
               for test_type in test_types)

    frame = pd.DataFrame(row_gen, columns=['a', 'b', 'test_type', 'p_value'])
    frame['p_value_adj'] = multipletests(frame['p_value'],
                                         method=corr_method)[1]
    frame.sort_values(by='p_value_adj', inplace=True)

    if threshold is not None:
        frame = frame.query('p_value_adj <= {}'.format(threshold))

    return frame
