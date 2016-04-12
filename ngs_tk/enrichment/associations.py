import itertools

import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


def contingency_table(a, b):
    """Calculates a 2x2 contigency table for two boolean Series.

    Args:
        a (pandas.Series): Boolean series a.
        b (pandas.Series): Boolean series b.

    Returns:
        pd.DataFrame: A 2x2 contingency table for a and b.

    """

    if a.name == b.name:
        raise ValueError('a and b should not have the same name ({})'
                         .format(a.name))

    # Concatenate both series (should take care of matching indices).
    # Note we drop any NAs to only keep matching entries.
    data = pd.concat([a, b], axis=1).dropna()

    # Create 2D contingency table.
    table = pd.pivot_table(data, index=data.columns[0],
                           columns=data.columns[1], aggfunc=len)
    table = table.fillna(0)

    return table


def test_association(a, b, alternative='two-sided'):
    """Tests for an association between series a and b using fishers exact.

    Args:
        a (pandas.Series): Boolean series a.
        b (pandas.Series): Boolean series b.
        alternative (str): Alternative hypothesis to test.

    Returns:
        float: p-value from the fishers exact test, reflecting the degree
            of association between series a and b.

    """
    table = contingency_table(a, b)
    return fisher_exact(table, alternative=alternative)


def test_associations(a, b, alternative='two-sided',
                      corr_method='fdr_bh', labels=('a', 'b')):
    """Tests for significant associations between columns in dataframes a and b.

    Args:
        a (pandas.DataFrame): DataFrame a. Should either be boolean or contain
            categorical values that can be expanded to a boolean matrix.
        b (pandas.DataFrame): DataFrame b, same format as a.
        alternative (str): Alternative hypothesis to test.
        corr_method (str): Correction method to use for to correct p-values
            for multiple testing.
        labels (str): Labels to use to in result.

    Returns:
        pd.DataFrame: Result dataframe containing associations and
            corresponding (corrected) p-values.

    """

    # Ensure a and b are boolean.
    a = convert_to_boolean(a)
    b = convert_to_boolean(b)

    # Combinations to try.
    combinations = itertools.product(a.iteritems(), b.iteritems())

    # Generate result frame from combinations.
    rows = ((a, b, test_association(a_values, b_values,
                                    alternative=alternative)[1])
            for (a, a_values), (b, b_values) in combinations)

    col_names = labels + ('p_value',)
    result = pd.DataFrame.from_records(rows, columns=col_names)

    # Apply multiple testing correction.
    result['p_value_corr'] = multipletests(
        result['p_value'], method=corr_method)[1]

    return result


def convert_to_boolean(data):
    """Converts categorical series or dataframe to boolean."""

    if isinstance(data, pd.DataFrame):
        return _convert_frame_to_boolean(data)
    elif isinstance(data, pd.Series):
        return _convert_series_to_boolean(data)
    else:
        raise ValueError('Unknown data type ({})'.format(type(data)))


def _convert_series_to_boolean(series):
    """Converts categorical series to boolean matrix."""

    if series.dtype.name == 'bool':
        return series
    else:
        expanded = pd.pivot_table(series.reset_index(), index=series.index.name,
                                  columns=series.name, aggfunc=len)
        return expanded > 0


def _convert_frame_to_boolean(frame):
    """Converts categorical DataFrame to boolean matrix."""

    converted = []
    for col_name, col_values in frame.items():
        if col_values.dtype.name != 'bool':
            conv = convert_to_boolean(col_values)
            conv.columns = ['{}_{}'.format(col_name, c)
                            for c in conv.columns]
            converted.append(conv)
        else:
            converted.append(col_values)

    return pd.concat(converted, axis=1)
