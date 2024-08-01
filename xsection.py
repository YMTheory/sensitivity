# We're gonna use nu_e + {^A_54} Xe -> {^A_55} *Cs + e- to measure the electron neutrinos

import numpy as np
import scipy as sp
import scipy.special
from unit_conversion import *

from scipy import integrate

class xsection:
    def __init__(self):
        self.GF = 1.1663787e-5 * 1e-6 # MeV^-2
        self.theta_c = 13.2/180*np.pi
        self.Zi =54 # initial state atomic number
        self.Zf = 55 # final state atomic number
        self.alpha = 0.0072973525643 # fine structure constant
        self.me = 0.510998910 # electron mass in MeV
        self.e = 0.303 #1.602176565e-19 # electron charge in Coulomb
        self.A = 136 # atomic number
        self.Rf = fm_to_MeVminus1(1.2 * np.power(self.A, 1/3.) ) # fm

        self.load_flag = False
        
    ## I am gonna calculate the differential and total cross sections by my own later
    ## Currently, I am stuck with how to calculate the multipole operators in Coulomb longitudinal and transverse components.
    def Fermi_function(self, El):
        # Here, all momentum and energy are replaced by corresponding effective values.
        VC0 = - self.Zf * self.alpha / 2. / self.Rf
        Eeff = El - VC0
        gamma = np.sqrt(1 - (self.Zf * self.alpha)**2)
        peff = np.sqrt(Eeff**2 - self.me**2)
        y = self.alpha * self.Zf * Eeff / peff
        
        F = 2 * (1+gamma) * np.power(2*peff*self.Rf, 2*(gamma-1)) * np.exp(np.pi*y) * abs(scipy.special.gamma(gamma+y*1j))**2 / abs(scipy.special.gamma(2*gamma+1))**2
        
        return F
      
        
    def transfer_momentum(self, Enu, El, theta):
        pl = np.sqrt(El**2 - self.me**2)
        q = np.sqrt( (Enu-pl)**2 + 4*Enu * pl * np.sin(theta / 2.)**2 )
        return q

    def longitudinal_component(self, Enu, omega, J):
        pass
    
    def transverse_component(self, Enu, omega, J):
        pass
    

    ### Load the data points from Brian's paper and interpolate to get the total cross section
    def load_csv(self):
        self.Ev1, self.xsec1 = np.loadtxt('/p/lustre1/yu47/Sterile_Neutrino/sensitivity/data/Cs136_1+590keV_totXsec.csv', delimiter=',', unpack=True)
        self.Ev2, self.xsec2 = np.loadtxt('/p/lustre1/yu47/Sterile_Neutrino/sensitivity/data/Cs136_1+850keV_totXsec.csv', delimiter=',', unpack=True)
        self.load_flag = True

    def interp_xsec_state1(self, E):
        if not self.load_flag:
            self.load_csv()
        return np.interp(E, self.Ev1, self.xsec1)
    
    def interp_xsec_state2(self, E):
        if not self.load_flag:
            self.load_csv()
        return np.interp(E, self.Ev2, self.xsec2)

        
        
    ### Directly use the total cross section from Brian's paper
    def total_xsec(self, Ev):
        pass
    
    
    
    def differential_xsec_nuescatter(self, Enu, Te):
        # neutrino energy in unit of MeV
        Q_plus = 0.231
        Q_minus = 0.5 + Q_plus
        dsigmadTe = 2 * self.GF**2 * self.me / np.pi * (Q_minus**2 + Q_plus**2 * (1 - Te/Enu)**2 - Q_minus*Q_plus*self.me*Te/Enu**2 )
        return dsigmadTe * self.Zi
    
    
    def total_xsec_nuescatter(self, Enu):
        Te_max = 2*Enu**2 / (2*Enu+self.me)
        f = lambda Te, E: self.differential_xsec_nuescatter(E, Te)
        res, err = integrate.quad(f, 0, Te_max, args=(Enu))
        res = MeV_to_cm(1.0)**2 * res
        return res