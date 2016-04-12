from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

# noinspection PyUnresolvedReferences
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

import numpy as np


def normalize_counts(counts):
    size_factors = estimate_size_factors(counts)
    return counts / size_factors


def estimate_size_factors(counts):
    log_geo_means = np.mean(np.log(counts), axis=1)
    sf = np.apply_along_axis(_estimate_size_factors_col, axis=0,
                             arr=counts, log_geo_means=log_geo_means)
    return sf


def _estimate_size_factors_col(counts, log_geo_means):
    log_counts = np.log(counts)
    mask = np.isfinite(log_geo_means) & (counts > 0)
    return np.exp(np.median((log_counts - log_geo_means)[mask]))
