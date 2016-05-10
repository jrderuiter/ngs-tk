import numpy as np


def sample_ad(sample):
    return sample.data.AD[1]


def sample_af(sample):
    try:
        return sample.data.AD[1] / sample.data.DP
    except (ZeroDivisionError, TypeError):
        return np.nan


def sample_dp(sample):
    return sample.data.DP
