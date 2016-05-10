
import numpy as np
from matplotlib import pyplot as plt


def volcano_plot(result, min_fc=2, max_pval=0.05, ax=None,
                 fc_col='log2FoldChange', pval_col='padj'):
    if ax is None:
        _, ax = plt.subplots()

    is_sign = ((result[fc_col].abs() > min_fc) &
               (result[pval_col] < max_pval))

    ax.plot(result.ix[is_sign, fc_col],
            -np.log2(result.ix[is_sign, pval_col]),
            'o', alpha=0.6)

    ax.plot(result.ix[~is_sign, fc_col],
            -np.log2(result.ix[~is_sign, pval_col]),
            'o', color='grey', alpha=0.6)

    xlim = np.ceil(result[fc_col].abs().max())
    ax.set_xlim(-xlim, xlim)

    ax.set_xlabel('Fold change')
    ax.set_ylabel('-log10(p-value)')

    return ax
