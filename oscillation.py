import numpy as np
import matplotlib.pyplot as plt
import scipy as sp


def electron_neutrino_survival_probability(dm2, sin2theta_square, E, L):
    survival_prob = 1 - sin2theta_square * np.sin(1.27*dm2*L/E)**2
    return survival_prob

    
def maximum_oscillation_energy_with_fixed_baseline(dm2, L, n):
    E = (1.27*dm2*L) / (np.pi /2. + n*np.pi)
    return E 

def maximum_oscillation_energy_with_fixed_energy(dm2, E, n):
    L = (np.pi/2+n*np.pi) / (1.27*dm2/E)
    return L

def oscillation_length(dm2, E):
    L = np.pi / (1.27 * dm2 /E)
    return L