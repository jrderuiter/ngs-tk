# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

import os
import shutil
import tempfile
import warnings

try:
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr

    pandas2ri.activate()
except ImportError:
    importr = None


def pathview(values, species, pathway_id, output_file=None, **kwargs):
    """Plots data values on genes in given KEGG pathway using pathview.

    Args:
        values (pandas.Series): Values to be plotted on the KEGG pathway. Should
            be indexed using the KEGG entrez ids of the corresponding genes.
        species (str): Name of the KEGG species.
        pathway_id (str): ID of the desired KEGG pathway.
        output_file (str): Output path for the drawn pathway.

    Returns:
        bytes or str: Path to output file if output_file is given, else
            a bytes object containing the PNG data of the plot.

    """
    if importr is None:
        raise ValueError('Rpy2 must be installed for pathview')

    # Load pathview in R.
    r_pathview = importr('pathview')

    # Run pathview within temp_dir and with disabled warnings.
    with tempfile.TemporaryDirectory() as tmp_dir:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            r_pathview.pathview(values, species=species, pathway_id=pathway_id,
                                kegg_dir=tmp_dir, **kwargs)

    # Get path to result.
    pathview_path = '{}{}.pathview.png'.format(species, pathway_id)

    # Handle returning of output.
    if output_file is None:
        # Read and return PNG if no file is given.
        with open(pathview_path, 'rb') as file_:
            data = file_.read()
            os.unlink(pathview_path)
            return data
    else:
        # Move figure to given location.
        shutil.move(pathview_path, output_file)
        return output_file
