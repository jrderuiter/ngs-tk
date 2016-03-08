
import pandas as pd

from .annovar import summarize_annovar
from . import stats


def summarize_vcf(records, summ_funcs):
    """Builds summary frames of given VCF records."""

    # Prepend index function to summary funcs.
    summ_funcs = [_index_func] + summ_funcs

    # Process each record.
    # TODO: Change to generator.
    rows = []
    for rec in records:
        rows.append(tuple(func(rec) for func in summ_funcs))

    # Split index and summary records.
    index, *summs = zip(*rows)

    # Build index as MultiIndex.
    index = pd.MultiIndex.from_tuples(
        index, names=['chrom', 'pos', 'ref', 'alt'])

    # Convert to dataframe(s).
    frames = (pd.DataFrame.from_records(summ, index=index)
              for summ in summs)

    return tuple(frames)


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
