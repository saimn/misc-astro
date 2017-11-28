# -*- coding: utf-8 -*-

import numpy as np


def generate_moffat_im(center=(12, 12), shape=(25, 25), alpha=2, beta=2.5):
    """Generate an image with a Moffat profile."""
    xx, yy = np.indices(shape)
    res = np.sqrt(((xx - center[0])**2 + (yy - center[1])**2))
    res = (1 + (res / alpha)**2)**(-beta)
    res /= np.sum(res)
    return res
