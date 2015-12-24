import numpy as np
import pandas as pd


def resample(frame, reference, bin_size,
             chrom='chrom', position='position', columns=None):
    """Resamples frame values to bins."""

    def _resample_chrom(chrom_id, grp, columns):
        chrom_length = len(reference[chrom_id])

        # Divide frame into position bins.
        bins = range(0, chrom_length + bin_size, bin_size)
        binned = pd.cut(grp[position], bins=bins,
                        include_lowest=True, right=False)

        # Calculate bin locations and means.
        binned_locs = np.arange(0, chrom_length, bin_size) + (bin_size / 2)
        binned_means = grp.groupby(binned)[columns].mean()

        # Add chrom and position to means.
        binned_means[chrom] = chrom_id
        binned_means[position] = binned_locs

        return binned_means

    # Select all columns by default.
    if columns is None:
        columns = [c for c in frame if c not in {chrom, position}]

    # Resample per chromosome and concatenate results.
    resampled = pd.concat((_resample_chrom(chrom_id, grp, columns)
                           for chrom_id, grp in frame.groupby(chrom)),
                          ignore_index=True)

    # Convert 'chrom' to category if needed.
    if frame[chrom].dtype == 'category':
        order = frame[chrom].cat.categories
        resampled[chrom] = pd.Categorical(resampled[chrom], categories=order)

    # Reorder columns.
    resampled = resampled[[chrom, position] + columns]

    return resampled
