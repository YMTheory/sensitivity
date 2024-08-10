import os
import numpy as np
from iminuit import cost, Minuit
import histlite as hl
import pickle

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
        self.pdf_filename = None
        

    def _set_datafilename(self, filename):
        self.data_filename = filename
        
    def _set_pdffilename(self, filename):
        self.pdf_filename = filename

        
    def _set_dm2(self, dm2):
        self.MC_gen.dm2 = dm2
        self.dm2 = dm2
    
    def _set_sin2theta_square(self, sin2):
        self.MC_gen.sin2theta_square = sin2
        self.sin2theta_square = sin2

        
        
    def generate_PDF(self, dm2, sin2theta_square, Enu):
        self.MC_gen.dm2 = dm2
        self.MC_gen.sin2theta_square = sin2theta_square
        
        self.MC_gen.n_events = 1e6 # Hard-coded for now, enough statistics.
        self.scale_count_flag = False
        _, bl = self.MC_gen.generate_oscillate_Asimov_dataset(Enu)
        
        h = hl.hist(bl, bins=self.n_Lbins, range=(self.Lmin, self.Lmax))
        h_norm = h.normalize()
        return h_norm
    
    
    def write_PDFs_intoFile(self, PDFs, filename):
        write_file = True
        if os.path.exists(filename):
            write_file = False
            user_input = input(f"{filename} already exists, do you want to overwrite it?").strip().lower()
            if user_input == 'yes':
                write_file = True        
        if write_file:
            with open(filename, 'wb') as file:
                pickle.dump(PDFs, file)


    def load_PDFs_fromFile(self):
        filename = self.pdf_filename
        try:
            with open(filename, 'rb') as file:
                self.PDF = pickle.load(file)
            return self.PDF
        except:
            print(f"{filename} does not exists.")
    
    
    def generate_dataset(self, dm2, sin2theta_square, Enu):
        self.MC_gen.dm2 = dm2
        self.MC_gen.sin2theta_square = sin2theta_square
        if not self.scale_count_flag:
            self.MC_gen.scale_counts()
            self.scale_count_flag = True
        
        _, bl = self.MC_gen.generate_oscillate_Asimov_dataset(Enu)
        
        return bl

    def write_data_histo_intoFiles(self, bls, filename):
        write_file = True
        if os.path.exists(filename):
            write_file = False
            user_input = input(f"{filename} already exists, do you want to overwrite it?").strip().lower()
            if user_input == 'yes':
                write_file = True        
        if write_file:
            hists = []
            for bl in bls:
                h = hl.hist(bl, bins=self.n_Lbins, range=(self.Lmin, self.Lmax))
                hists.append(h)
            with open(filename , 'wb') as file:
                pickle.dump(hists, file)
            
    
    def load_data_fromFile(self):
        filename = self.data_filename
        try:
            with open(filename, 'rb') as file:
                self.dataset = pickle.load(file)
            return self.dataset
        except:
            print(f"{filename} does not exists.")
            

    def expected_rate(self, bl, N):
        if self.PDF is None:
            self.load_PDFs_fromFile()
            return N * self.PDF.get_value(bl)

            
    def fit_statiscs_only_oneEvent(self, evno):
        if len(self.dataset) == 0:
            self.load_data_fromFile()

        if evno > len(self.dataset):
            evno = len(self.dataset) - 1
            
        # The raw dataset is array, need to fit it into histograms
        h_data = self.dataset[evno]
        
        c = cost.BinnedNLL(h_data.values, h_data.bins[0], self.expected_rate)
        m = Minuit(c, N=4000,)
        m.limits['N'] = (3000, 5000)
        m.migrad()
        
        return m
            

    def fit_statiscs_only_oneEvent_manually(self, evno, Nmin=300, Nmax=500):
        if len(self.dataset) == 0:
            self.load_data_fromFile()
            
        if evno > len(self.dataset):
            evno = len(self.dataset) - 1
            
        h_data = self.dataset[evno]

        def get_chi2(N):
            chi2 = 0
            for i in range(h_data.n_bins[0]):
                obs = h_data.values[i]
                err = h_data.errors[i]
                pred = N * self.PDF.values[i]
                if err == 0:
                    continue
                chi2 += (obs - pred)**2 / err**2

            return chi2

        # Coarse scanning:
        step = 0.1
        coarse_step, coarse_chi2 = np.arange(Nmin, Nmax, step), []
        for Ni in coarse_step:
            tmp_chi2 = get_chi2(Ni)
            coarse_chi2.append(tmp_chi2)
        coarse_chi2 = np.array(coarse_chi2)
        N0 = coarse_step[np.argmin(coarse_chi2)]

        # Fine scanning:
        step = 0.01
        fine_step, fine_chi2 = np.arange(N0-5, N0+5, step), []
        for Ni in fine_step:
            tmp_chi2 = get_chi2(Ni)
            fine_chi2.append(tmp_chi2)
        fine_chi2 = np.array(fine_chi2)

        idx = np.argmin(fine_chi2)
        return fine_step[idx], fine_chi2[idx],coarse_step, coarse_chi2, fine_step, fine_chi2

        