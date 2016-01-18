from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

# noinspection PyUnresolvedReferences
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)


def parse_barcode(barcode):
    """Parse a tcga barcode into its constituents."""

    # Define fields and their locations in the barcode.
    fields = {
        'project': (0, ),
        'tss': (1, ),
        'participant': (2, ),
        'sample': (3, (0, 2)),
        'vial': (3, (2, )),
        'portion': (4, (0, 2)),
        'analyte': (4, (2, )),
        'plate': (5, ),
        'center': (6, )
    }

    # Split barcode and extract fields.
    split = barcode.split('-')

    values = {field: _extract_field(split, indices)
              for field, indices in fields.items()}

    return values


def _extract_field(split, indices):
    try:
        if len(indices) == 1:
            value = split[indices[0]]
        elif len(indices) == 2:
            ind1, ind2 = indices
            if len(ind2) == 1:
                value = split[ind1][ind2[0]]
            else:
                value = split[ind1][ind2[0]:ind2[1]]
        return value
    except IndexError:
        return None


def get_sample(barcode):
    parsed = parse_barcode(barcode)
    return parsed['sample']


def get_particpant(barcode):
    parsed = parse_barcode(barcode)
    return '-'.join([parsed['project'], parsed['tss'], parsed['participant']])
