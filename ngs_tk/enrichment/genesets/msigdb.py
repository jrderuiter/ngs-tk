# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401


def read_gmt(file_path):
    """Reads an msigdb file containing genesets."""
    with open(file_path, 'r') as file_:
        return dict(map(_parse_gmt_line, file_))

def _parse_gmt_line(line):
    split = line.strip().split()

    set_name = split[0]
    set_genes = set(split[2:])

    return set_name, set_genes
