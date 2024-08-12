import os
import numpy as np
from iminuit import cost, Minuit
import histlite as hl
import sys
import h5py as h5

from MC_generator import MC_generator


class fitter:
    def __init__(self, source, detector) -> None:
        self.source = source
        self.detector = detector

        self.MC_gen = MC_generator(self.source, self.detector)

        source_z = self.source.position[2]
        detector_center_z = self.detector.position[2]
        self.Lmin, self.Lmax, self.Lstep = detector_center_z - source_z - self.detector.height/2. ,  detector_center_z - source_z + self.detector.height/2., 0.03
        self.n_Lbins = int((self.Lmax - self.Lmin) / self.Lstep) + 1
        
        self.dm2 = 0.
        self.sin2theta_square = 0.

        self.scale_count_flag = False
        self.MC_gen.scale_counts()

        self.dataset = []
        self.PDF  = None
        self.data_filename = None
        self.load_data_flag = False
        self.pdf_filename = None
        self.load_pdf_flag = False
        

        self.data_dm2 = 0.0
        self.data_sin2 = 0.0
        self.fit_dm2 = 0.0
        self.fit_sin2 = 0.0

    ## Setters:
    def _set_data_dm2(self, val):
        self.data_dm2 = val
    
    def _set_data_sin2(self, val):
        self.data_sin2 = val
        
    def _set_fit_dm2(self, val):
        self.fit_dm2 = val
        
    def _set_fit_sin2(self, val):
        self.fit_sin2 = val

    def _set_data_filename(self, name):
        self.data_filename = name
        
    def _set_pdf_filename(self, name):
        self.pdf_filename = name
   ######################################################################## 
        
    
    ## Gettters:
    def _get_data_dm2(self):
        return self.data_dm2
    
    def _get_data_sin2(self):
        return self.data_sin2
    
    def _get_fit_dm2(self):
        return self.fit_dm2
    
    def _get_fit_sin2(self):
        return self.fit_sin2

   ######################################################################## 
        
    def load_data_file(self, evno=[]):
        if not os.path.exists(self.data_filename):
            print(f'{self.data_filename} does not exist :(')
            sys.exit(-1)

        with h5.File(self.data_filename, 'r') as f:

            if len(evno) == []:
                dataset_names = f.keys()
                for name in dataset_names:
                    self.dataset.append( f[name][:]) 
            
            else:
                for i in evno:
                    name = f'event{i}'
                    if name not in f.keys():
                        continue
                    self.dataset.append( f[name][:])
            
                    
                    
    
    
    def calculate_rateonly_dchi2_asimov(self):
        ## There is actually no only fitting, just calculating rate-only delta_chi2 for (sin2, dm2) pairs.
        # Fitting histogram
        h_fit = self.MC_gen.generate_asimov_dataset(data_dm2=self.fit_dm2, data_sin2=self.fit_sin2, )
        fit_nevt = np.sum( h_fit.values )
        
        # Data asimov histogram
        h_data = self.MC_gen.generate_asimov_dataset(data_dm2=self.data_dm2, data_sin2=self.data_sin2, )
        data_nevt = np.sum( h_data.values )

        # For now, only consider Poisson statistics fluctuation.
        sigma_data = np.sqrt( data_nevt )
        
        chi2 = (fit_nevt - data_nevt)**2 / sigma_data**2

        return chi2
        

    def calculate_shape_only_dchi2_asimov_fixedRate(self):
        ## There is actually no only fitting, just calculating rate-only delta_chi2 for (sin2, dm2) pairs.
        # Fitting histogram
        h_fit = self.MC_gen.generate_asimov_dataset(data_dm2=self.fit_dm2, data_sin2=self.fit_sin2, )
        fit_spec = h_fit.values
        
        # Data asimov histogram
        h_data = self.MC_gen.generate_asimov_dataset(data_dm2=self.data_dm2, data_sin2=self.data_sin2, )
        data_spec = h_data.values

        # For now, only consider Poisson statistics fluctuation.
        sigma_data = np.sqrt( data_spec )
        
        mask_id = np.where(sigma_data != 0)
        chi2 = np.sum( (fit_spec[mask_id] - data_spec[mask_id])**2 / sigma_data[mask_id]**2 ) 

        return chi2

    def calculate_shape_only_dchi2_asimov_scanRate(self, coarse_scan_step=0.01, coarse_scan=11, fine_scan_step=0.001, fine_scan=11):
        h_fit = self.MC_gen.generate_asimov_dataset(data_dm2=self.fit_dm2, data_sin2=self.fit_sin2, )
        fit_spec = h_fit.values
        
        # Data asimov histogram
        h_data = self.MC_gen.generate_asimov_dataset(data_dm2=self.data_dm2, data_sin2=self.data_sin2, )
        data_spec = h_data.values

        sigma_data = np.sqrt( data_spec )
        
        mask_id = np.where( sigma_data != 0 )

        ## Coarse scanning
        coarse_dchi2 = np.zeros(coarse_scan)
        coarse_R     = np.zeros(coarse_scan)
        N_oneSide = int((coarse_scan/2) -1)
        low = 1 - coarse_scan_step * N_oneSide
        for i in range(coarse_scan):
            Ri = low + coarse_scan_step * i
            coarse_R[i] =  Ri 
            dchi2 =  np.sum( (data_spec[mask_id] - Ri * fit_spec[mask_id] )**2 / sigma_data[mask_id]**2 ) 
            coarse_dchi2[i] = dchi2 

        coarse_min_R = coarse_R[np.argmin(coarse_dchi2) ]

        # Fine scanning
        fine_dchi2 = np.zeros(fine_scan)
        fine_R     = np.zeros(fine_scan)
        N_oneSide = int((fine_scan/2) -1)
        low = 1 - fine_scan_step * N_oneSide
        for i in range(fine_scan):
            Ri = low + fine_scan_step * i
            fine_R[i] =  Ri 
            dchi2 =  np.sum( (data_spec[mask_id] - Ri * fit_spec[mask_id] )**2 / sigma_data[mask_id]**2 ) 
            fine_dchi2[i] = ( dchi2 )

        fine_min_R = fine_R[np.argmin(fine_dchi2) ]
        fine_min_dchi2 = np.min( fine_dchi2 )
        
        
        return coarse_R, coarse_dchi2, fine_R, fine_dchi2