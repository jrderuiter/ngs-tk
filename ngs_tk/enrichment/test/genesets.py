from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

from itertools import chain
from toolz import curry

import pandas as pd
from scipy.stats import hypergeom, fisher_exact
from statsmodels.stats.multitest import multipletests as p_adj


def test_sets(selected_genes, gene_sets, all_genes=None, direction=None,
              corr_method=None, threshold=None, sort=True):
    if all_genes is None:
        all_genes = set(chain.from_iterable(gene_sets.values()))

    # Curry test function for convenience.
    test_func = curry(_test_sets, gene_sets=gene_sets, all_genes=all_genes,
                      corr_method=corr_method, sort=sort, threshold=threshold)

    if direction is None:
        # Do normal test.
        result = test_func(selected_genes)
    else:
        # Do test per-direction.
        genes_up = [g for g, d in zip(selected_genes, direction) if d == 1]
        genes_dn = [g for g, d in zip(selected_genes, direction) if d == -1]

        sets_up = test_func(genes_up).assign(direction='up')
        sets_down = test_func(genes_dn).assign(direction='down')

        result = pd.concat([sets_up, sets_down], ignore_index=True)

    return result


def _test_sets(selected_genes, gene_sets, all_genes,
               corr_method=None, threshold=None, sort=True):
    # Calculate p-values for each of the gene sets.
    rows = ((set_name, test_set_hypergeom(
             selected_genes, all_genes, set_genes))
            for set_name, set_genes in gene_sets.items())
    result = pd.DataFrame(rows, columns=['set_name', 'p_value'])

    # Sort if required. Note corrected p-value should be monotonic
    # with uncorrected, so we can sort just on the p-value.
    if sort:
        result.sort_values('p_value', ascending=True, inplace=True)

    # Calculate a corrected p-value if required.
    if corr_method is not None:
        result['p_value_adj'] = p_adj(result['p_value'], method=corr_method)[1]

    # Return only significant cases if threshold is specified.
    if threshold is not None:
        threshold_col = 'p_value_adj' if corr_method is not None else 'p_value'
        result = result.ix[result[threshold_col] < threshold]

    return result


def test_set_hypergeom(selected_genes, all_genes, set_genes):
    # Reduce the gene_set to the universe of all_genes,
    # as we can only sample from this set.
    all_genes = set(all_genes)
    set_genes = set(set_genes).intersection(all_genes)
    selected_genes = set(selected_genes).intersection(all_genes)

    # Calculate overlap of the gene_set with selected genes.
    selected_set_genes = set_genes.intersection(selected_genes)

    # Calculate p-value using the hyper-geometric test.
    p_val = hypergeom.sf(
        M=len(all_genes), n=len(set_genes),
        N=len(selected_genes), k=len(selected_set_genes), loc=1)

    return p_val


def test_set_fisher(selected_genes, all_genes, set_genes):
    # Define our sets to analyze.
    all_genes = set(all_genes)
    set_genes = set(set_genes).intersection(all_genes)
    selected_genes = set(selected_genes).intersection(all_genes)

    # Calculate intersections/differences.
    both = len(selected_genes.intersection(set_genes))
    sel_only = len(selected_genes.difference(set_genes))
    set_only = len(set_genes.difference(selected_genes))
    neither = len(all_genes.difference(selected_genes).difference(set_genes))

    # Calculate p-value using fishers exact test.
    odds, p_val = fisher_exact([[both, sel_only],
                                [set_only, neither]],
                               alternative='greater')

    return p_val
