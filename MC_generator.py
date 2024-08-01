import numpy as np
import hist
from hist import Hist

from signal_spectrum import signal_calculation
from detector import detector
from neutrino_source import neutrino_source
from oscillation import *
from xsection import xsection

from unit_conversion import *

class MC_generator:
    def __init__(self, source, det, energy_bins=[], dm2=1.0, sin2theta_square=0.1):
        self.energy_bins = energy_bins
        self.source = source
        self.det = det
        self.xsec = xsection()
        
        # neutrino mixing parameters:
        self.dm2 = dm2
        self.sin2theta_square = sin2theta_square

        # expected un-oscillated neutrino number:
        #self.source.scale_counts(self.det.height. self.det.volume, self.det.baseline, self.det.run_time )
        self.n_events = self.scale_counts()


    def scale_counts(self):
        # Scale the number of neutrinos generated by the source from the LZ paper (JHEP11(2014) 042)
        # Approximately take the average arrived neutrinos scaled by the square of distances from the source to the detector center (L+H/2.)**2
        # The LZ detector is hard-coded here.
        H0 = 1.38 # unit: m, LZ detector height
        R0 = H0 / 2. # unit: m, LZ detector radius
        V0 = np.pi * R0**2 * H0 # unit: m^3, LZ detector volume
        L0 = 1 # unit: m, below the fiducial volume, along the central axis of the cylinder. Distance from the source to the bottom plane of the cylinder.
        N0 = 12518 # 100 days run
        T0 = 100 # unit: days, 100 days run
        activity0 = 5e6 # unit: Ci, 1 Ci = 3.7e10 dacays per second
        
        H = self.det.height
        V = self.det.volume
        #L = self.det.baseline - self.det.height/2.
        L = np.sqrt((self.source.position[0]-self.det.position[0])**2 + (self.source.position[1]-self.det.position[1])**2 + (self.source.position[2]-self.det.position[2])**2 ) - self.det.height/2.
        T = self.det.run_time
        
        Enu = self.source.energies[0]
        br  = self.source.ratios[0]
        xsec_nue = self.xsec.total_xsec_nuescatter(Enu)
        xsec_cc = self.xsec.interp_xsec_state1(Enu)

        self.n_events = (N0 / V0 / T0 * (L0+ H0/2.)**2 ) * V * T / (L+H/2.)**2 * (self.source.activity / activity0) * br * xsec_cc / xsec_nue
        return self.n_events
    


    def _get_source(self):
        return self.gen
    
    def _get_detector(self):
        return self.det
    
    def _get_nevents(self):
        return self.n_events

    def generate_in_cubic(self, x0, y0, z0, H0, L0, N0):
        x = np.random.uniform(-L0/2+x0, L0/2+x0, size=N0)
        y = np.random.uniform(-L0/2+y0, L0/2+y0, size=N0)
        z = np.random.uniform(-H0/2+z0, H0/2+z0, size=N0)
        return x, y, z
    
    def Is_in_detector(self, x, y, z):
        if ( np.sqrt((x-self.det.position[0])**2 + (y-self.det.position[1])**2) < self.det.radius) and (-self.det.height/2. <= z-self.det.position[2] < self.det.height/2.) :
            return True
        else:
            return False

    def smear_position(self, pos, res):
        x_smeared = np.random.normal(loc=pos[:,0], scale=res, size=len(pos))
        y_smeared = np.random.normal(loc=pos[:,1], scale=res, size=len(pos))
        z_smeared = np.random.normal(loc=pos[:,2], scale=res, size=len(pos))
        pos_smeared = np.array([x_smeared, y_smeared, z_smeared]).T
        bl = np.sqrt( (x_smeared - self.source.position[0])**2 + (y_smeared - self.source.position[1])**2 + (z_smeared - self.source.position[2])**2 )
        return pos_smeared, bl


    def generate_nonoscillate_Asimov_dataset(self, smear=False):
        nu_pos = []
        bl = []
        n_event0 = int(2 * self.n_events)
        x0, y0, z0 = self.generate_in_cubic(*self.det.position, self.det.height, 2*self.det.radius, n_event0)
        for x, y, z in zip(x0, y0, z0) :
            if self.Is_in_detector(x, y, z) :
                nu_pos.append( (x, y, z) )
                bl.append( np.sqrt((x - self.source.position[0])**2 + (y - self.source.position[1])**2 + (z - self.source.position[2])**2) )
            else:
                continue
        if not smear:
            return np.array(nu_pos[0:int(self.n_events)]), np.array(bl[0:int(self.n_events)])
        else:
            pos_smeared, bl_smeared = self.smear_position(np.array(nu_pos[0:int(self.n_events)]), self.det.spatial_resolution)
            return pos_smeared, bl_smeared
        
        
        
    def generate_oscillate_Asimov_dataset(self, Enu, smear=False):
        nu_pos0, bl0 = self.generate_nonoscillate_Asimov_dataset()
        sur_p = electron_neutrino_survival_probability(self.dm2, self.sin2theta_square, Enu, bl0)
        p0 = np.random.uniform(0, 1, size=len(sur_p))
        sur_id = np.where(p0 < sur_p)[0]
        nu_pos, bl = [], []
        for id in sur_id:
            nu_pos.append( nu_pos0[id] ) 
            bl.append(bl0[id])
        
        if not smear:
            return np.array(nu_pos), np.array(bl)
        else:
            pos_smeared, bl_smeared = self.smear_position(np.array(nu_pos), self.det.spatial_resolution)
            return pos_smeared, bl_smeared
    
    
    
        

            