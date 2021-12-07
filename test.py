from Utils import *
import numba
import numpy as np
import time

def generate_identifier_fast(grams, forward, backward):
    grams_size = len(grams)
    return_vec = np.zeros(2 * grams_size, dtype=bool)
    for index, g in enumerate(grams):
        if g in forward:
            return_vec[index] = 1
        if g in backward:
            return_vec[index + grams_size] = 1
    return return_vec

@
