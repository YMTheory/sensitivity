from oscillation import *
from unit_conversion import *

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import nquad
from scipy.stats import norm


class signal_calculation:

    def __init__(self, gen, det, energy_bins, dm2, sin2theta_square) -> None:
        self.gen = gen
        self.det = det
        self.n_events = 12518
        self.energy_bins = energy_bins
        self.dm2 = dm2
        self.sin2theta_square = sin2theta_square
        self.Enu = 0.75


    def _set_dm2(self, dm2):
        self.dm2 = dm2
        
    def _set_sin2theta_square(self, sin2theta_square):
        self.sin2theta_square = sin2theta_square
        
    def _set_energy_bins(self, energy_bins):
        self.energy_bins = energy_bins

    def total_events(self):
        self.n_events = day_to_s(self.det.livetime ) * self.gen.activity * (self.det.radius*2 * self.det.height) / (4*np.pi*self.det.baseline)


    def integrand(self, z, r, theta, x0, Ev):
        x = x0 + r * np.cos(theta)
        y = r * np.sin(theta)
        bl = np.sqrt(x**2+y**2+z**2)
        sur_prob = electron_neutrino_survival_probability(self.dm2, self.sin2theta_square, Ev, bl)
        return sur_prob * r

    def normalized_oscillated_electron_energy_spectrum(self):
        energy_spec = np.zeros(len(self.energy_bins)-1)
        width = self.energy_bins[1] - self.energy_bins[0]
        cents = (self.energy_bins[:-1] + self.energy_bins[1:]) / 2.
        probs, specs = [], []
        theta_limits = [0, 2 * np.pi]
        r_limits = [0, self.det.radius]
        z_limits = [-self.det.height/2., self.det.height/2.]

        ## Here we are considering the source as ideal mono-energetic neutrino source with discrete neutrino energies.
        for R, Enu in zip(self.gen.ratios, self.gen.energies):
            res, err = nquad(self.integrand, [z_limits, r_limits, theta_limits], args=(self.det.baseline, Enu))
            probs.append(res / self.det.volume * R )
            resol = self.det.empirical_combined_energy_resolution(Enu)
            tmp_spec = probs[-1] * norm.pdf(cents, loc=Enu, scale=resol*Enu) * width
            specs.append(tmp_spec)
        
            energy_spec += tmp_spec
        
        return probs, specs, energy_spec
            

        