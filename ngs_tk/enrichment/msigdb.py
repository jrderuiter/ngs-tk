from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

# noinspection PyUnresolvedReferences
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)


def read_gmt(file_path):
    with open(file_path, 'r') as file_:
        gene_sets = dict(map(_parse_gmt_line, file_))
    return gene_sets


def _parse_gmt_line(line):
    split = line.strip().split()

    set_name = split[0]
    set_genes = set(split[2:])

    return set_name, set_genes
