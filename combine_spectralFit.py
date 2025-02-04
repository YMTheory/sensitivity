import numpy as np
import pandas as pd
import pickle
import histlite as hl
from iminuit import Minuit
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import sys
sys.path.append("/fs/ddn/sdf/group/nexo/users/miaoyu/Sterile_nu/sensitivity")

from spectralFit import spectralFit
from oscillation import *

class combine_spectralFit:
    def __init__(self) -> None :
        self.fitters = []

        # systematic budget:
        ## relative percentage
        self.sigma_flux                 = 0.02
        self.sigma_xsec_CC              = 0.13
        self.sigma_xsec_ES              = 0.003
        self.sigma_efficiency           = 0.01
        self.sigma_background           = 0.005
        self.alpha_flux                 = 0
        self.alpha_xsec_CC              = 0
        self.alpha_xsec_ES              = 0
        self.alpha_efficiency           = 0
        self.alpha_background           = 0 

        self.all_parameters             = ['alpha_flux', 'alpha_xsec_CC', 'alpha_xsec_ES', 'alpha_efficiency', 'alpha_background']
        self.fixed_parameters           = []
        
        ## Fitting results:
        self.minimizer                  = None
        self.fitted_values              = None
        self.fitted_errors              = None
        self.fval                       = 0.
        
    def _empty_fitter_list(self):
        self.fitters = []
        
    def _add_fitter(self, fitter):
        self.fitters.append( fitter )
        
    def _set_fixed_parameters(self, params):
        self.fixed_parameters = params
    
        
    def combine_chi_square(self, alpha_flux, alpha_xsec_CC, alpha_xsec_ES, alpha_efficiency, alpha_background):
        dchi2 = 0.
        for fitter in self.fitters:
            if fitter.channel == 'ES':
                dchi2 += fitter.chi_square(alpha_flux, alpha_xsec_ES, alpha_efficiency, alpha_background)
            elif fitter.channel == 'CC':
                dchi2 += fitter.chi_square(alpha_flux, alpha_xsec_CC, alpha_efficiency, alpha_background)
        # The shared uncertainties shoud be subtracted from it.
        dchi2 -= alpha_flux**2/self.sigma_flux**2
        dchi2 -= alpha_efficiency**2 / self.sigma_efficiency**2
        return dchi2
    
    def minimize_combine_chi_square(self, alpha_flux0=0.0, alpha_xsec_CC0=0.0, alpha_xsec_ES0=0.0, alpha_efficiency0=0.0, alpha_background0=0.0):
        is_fit_valid = False
        n_fit, n_fit_max0 = 0, 10
        
        while (not is_fit_valid) and (n_fit < n_fit_max0):
            # initial values of fitting parameters
            initials = { \
                        'alpha_flux' :       np.random.normal(alpha_flux0, 0.001), \
                        'alpha_xsec_CC' :    np.random.normal(alpha_xsec_CC0, 0.005), \
                        'alpha_xsec_ES' :    np.random.normal(alpha_xsec_ES0, 0.001), \
                        'alpha_efficiency' : np.random.normal(alpha_efficiency0, 0.001), \
                        'alpha_background' : np.random.normal(alpha_background0, 0.001) \
                        }
        
            for k in self.fixed_parameters:
                initials[k]     = 0.0

            m = Minuit(self.combine_chi_square, **initials)
            for k in self.all_parameters:
                if k in self.fixed_parameters:
                    m.fixed[k]      = True
                else:
                    m.limits[k] = (-10.0, 10.0)
        
            
            m.migrad()
            is_fit_valid = m.valid
            n_fit += 1

        if not is_fit_valid:
            print(f'Fitting fails after {n_fit} trys, will suspend it for now.')
    
        self.alpha_flux         = m.values['alpha_flux']
        self.alpha_xsec_CC      = m.values['alpha_xsec_CC']
        self.alpha_xsec_ES      = m.values['alpha_xsec_ES']
        self.alpha_efficiency   = m.values['alpha_efficiency']
        self.alpha_background   = m.values['alpha_background']

        self.minimizer          = m
        self.fval               = m.fval
        self.fitted_values      = m.values
        self.fitted_errors      = m.errors
    

        


