import numpy as np
from math import sqrt

import numba


def preranked_running_sum(ranking, sets, **kwargs):
    raise NotImplementedError()


def running_sum(vec, n_perm=None):
    es, rs = _running_sum(vec)

    if n_perm is not None:
        perm_es = []

        high = len(vec)
        size = vec.sum()

        for _ in range(n_perm):
            rand_ind = np.random.randint(high, size=size)
            rand_vec = np.bincount(rand_ind, minlength=high)
            rand_es, rand_rs = _running_sum(rand_vec)
            perm_es.append(rand_es)

        perm_es = np.array(perm_es)

        if es >= 0:
            p_val = np.sum(perm_es >= es) / n_perm
        else:
            p_val = np.sum(perm_es <= es) / n_perm

        return es, rs, p_val, perm_es
    else:
        return es, rs


@numba.jit(nopython=True)
def _running_sum(vec):
    rs = np.zeros_like(vec, dtype=np.float64)

    Ns = np.sum(vec)
    N = len(vec)

    increment = sqrt((N - Ns) / Ns)
    decrement = sqrt(Ns / (N - Ns))

    for i in range(1, len(rs)):
        if vec[i] > 0:
            rs[i] = rs[i-1] + (vec[i] * increment)
        else:
            rs[i] = rs[i-1] - decrement

    # ES score is maximum deviation from zero.
    es_max = rs.max()
    es_min = rs.min()

    es = es_min if abs(es_min) > es_max else es_max

    return es, rs
