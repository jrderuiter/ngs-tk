from itertools import chain

import pandas as pd

from .annovar import summarize_annovar, get_annovar_fields
from .snpeff import summarize_snpeff, get_snpeff_fields
from . import stats


def summarize_vcf(records, summ_funcs):
    """Builds summary frames of given VCF records."""

    # Prepend index function to summary funcs.
    summ_funcs = [_index_func] + summ_funcs

    # Process each record.
    rows = (tuple(func(rec) for func in summ_funcs)
            for rec in records)
    index, *summaries = zip(*rows)

    # Build index as MultiIndex.
    index = pd.MultiIndex.from_tuples(
        index, names=['chrom', 'pos', 'ref', 'alt'])

    # Convert to dataframe(s).
    frames = []
    for summ in summaries:
        if isinstance(summ[0], (list, tuple)):
            # Handle nested list case.
            with_index = (([i] * len(s), s)
                          for i, s in enumerate(summ)
                          if s is not None)
            row_indices, row_summs = zip(*with_index)

            # Flatten nested values.
            row_indices = chain.from_iterable(row_indices)
            row_summs = chain.from_iterable(row_summs)
        else:
            # Handle basic entry case.
            with_index = ((i, s) for i, s in enumerate(summ)
                          if s is not None)
            row_indices, row_summs = zip(*with_index)

        # Create frame using summaries and index entries.
        frame = pd.DataFrame.from_records(
            row_summs, index=index[list(row_indices)])
        frames.append(frame)

    return frames


def _index_func(rec):
    """Builds index row for summarize_vcf."""
    return (rec.CHROM, rec.POS, rec.REF, str(rec.ALT[0]))


def summarize_ad(rec):
    """Summarizes sample ADs for record."""
    return {s.sample: stats.sample_ad(s) for s in rec.samples}


def summarize_af(rec):
    """Summarizes sample AFs for record."""
    return {s.sample: stats.sample_af(s) for s in rec.samples}


def summarize_dp(rec):
    """Summarizes sample depths (DP) for record."""
    return {s.sample: stats.sample_dp(s) for s in rec.samples}
