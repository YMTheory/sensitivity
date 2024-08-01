import numpy as np

def fm_to_MeVminus1(l):
    return l * 1 / 197.3

def day_to_s(t):
    return t * 24 * 3600.

def MeV_to_cm(E):
    return 0.197e-10 # cm