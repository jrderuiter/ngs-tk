
import numpy as np

from .stats import sample_af


def filter_normals(records, normals, max_af=0.2):
    for rec in records:
        afs = np.array([sample_af(rec.genotype(s)) for s in normals])
        if not afs.max() > max_af:
            yield rec
