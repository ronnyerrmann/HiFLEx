import os
import numpy as np


def nanminmax(data, axis=None):
    """
    returns the minimum and maximum of a list or array
    """
    return [np.nanmin(data, axis=axis), np.nanmax(data, axis=axis)]

def round_sig(value, sig=5):
    """
    Rounds to a significant number, e.g. to have 5 digits
    :param value: number
    :param sig: int, number of dgits allowed
    :return rounded: float, rounded to the number of sig digits
    """
    if np.isnan(value) or np.isinf(value):
        return value
    if value == 0.0:
        return value
    rounded = round(value, sig-int(np.floor(np.log10(abs(value))))-1)
    
    return rounded
