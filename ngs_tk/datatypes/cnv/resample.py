from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

# noinspection PyUnresolvedReferences
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

import numpy as np
import pandas as pd


def resample(data, reference, bin_size):
    """Resamples frame values to bins."""

    # Drop infinite values.
    data = data.replace([np.inf, -np.inf], np.nan)

    # Resample per chromosome and concatenate results.
    resampled = (_resample_chrom(chrom, grp, reference, bin_size)
                 for chrom, grp in data.groupby(level=0))
    resampled = pd.concat(resampled)

    # Drop any resulting NA rows.
    resampled.dropna(inplace=True)

    return resampled


def _resample_chrom(chrom, grp, reference, bin_size):
    chrom_length = len(reference[chrom])

    # Ensure group is sorted.
    grp = grp.sort_index()

    # Divide frame into position bins.
    bins = range(0, chrom_length + bin_size, bin_size)
    binned = pd.cut(grp.index.get_level_values(1), bins=bins,
                    include_lowest=True, right=False)

    # Calculate bin locations and means.
    binned_locs = np.arange(0, chrom_length, bin_size) + (bin_size / 2)
    binned_means = grp.groupby(binned).mean()

    # Add chrom and position to means.
    binned_means.index = pd.MultiIndex.from_arrays(
        [[chrom] * len(binned_locs), binned_locs],
        names=grp.index.names)

    return binned_means
