
def annotate_counts_barplot(data, x, y, ax, **kwargs):
    """Annotates an existing barplot with counts above the bars."""

    # Update default annotate keywords with any kwargs.
    ann_kws = {'xytext': (0, 3),
               'ha': 'center',
               'textcoords': 'offset points'}
    ann_kws.update(kwargs)

    # Create lookup for x position of label.
    x_pos = {label.get_text(): tick for tick, label in
             zip(ax.get_xticks(), ax.get_xticklabels())}

    # Annotate with counts!
    for x_val, y_val in zip(data[x], data[y]):
        ax.annotate(xy=(x_pos[x_val], y_val), s=y_val, **ann_kws)

    return ax
