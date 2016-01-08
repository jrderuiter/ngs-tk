from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

from contextlib import contextmanager

from rpy2 import robjects
from rpy2.robjects import pandas2ri
from rpy2.rinterface import RNULLType, RRuntimeError

from rpy2.robjects.packages import importr as rpy2_importr

pandas2ri.activate()

base = rpy2_importr('base')


class RDependencyError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


def pandas_to_dataframe(pd_frame, check_names=False):
    r_frame = pandas2ri.py2ri_pandasdataframe(pd_frame)
    r_frame.colnames = pd_frame.columns

    if not check_names:
        r_frame.rownames = pd_frame.index

    return r_frame


def dataframe_to_pandas(r_frame):
    pd_frame = pandas2ri.ri2py_dataframe(r_frame)

    # Extract column names if possible.
    col_names = robjects.r.colnames(r_frame)
    if not type(col_names) == RNULLType:
        pd_frame.columns = col_names

    # Extract row names if possible.
    index = robjects.r.rownames(r_frame)
    if not type(index) == RNULLType:
        pd_frame.index = index

    return pd_frame


def importr(pkg_name):
    try:
        pkg = rpy2_importr(pkg_name)
    except RRuntimeError:
        raise RDependencyError(
            'Required R package {} is not installed'.format(pkg_name))
    else:
        return pkg


@contextmanager
def change_dir(dir_path):
    curr_wd = base.getwd()
    base.setwd(dir_path)
    yield
    base.setwd(curr_wd)
