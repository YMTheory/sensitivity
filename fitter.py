import numpy as np
import iminuit
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
        self.Lmin, self.Lmax, Lstep = detector_center_z - source_z - self.detector.height/2. ,  detector_center_z - source_z + self.detector.height/2., 0.03
        self.n_Lbins = int((self.Lmax - self.Lmin) / Lstep) + 1

        self.MC_gen.scale_counts()
        
    def generate_PDF(self, dm2, sin2theta_square, Enu):
        self.MC_gen.dm2 = dm2
        self.MC_gen.sin2theta_square = sin2theta_square
        
        self.MC_gen.n_events = 1e6 # Hard-coded for now, enough statistics.
        _, bl = self.MC_gen.generate_oscillate_Asimov_dataset(Enu)
        
        h = hl.hist(bl, bins=self.n_Lbins, range=(self.Lmin, self.Lmax))
        h_norm = h.normalize()
        return h_norm
    
    
    def write_PDFs_intoFile(self, PDFs, filename):
        with open(filename, 'wb') as file:
            pickle.dump(PDFs, file)


    def load_PDFs_fromFile(self, filename):
        with open(filename, 'rb') as file:
            PDFs = pickle.load(file)
        return PDFs
    
    
    def generate_dataset(self, dm2, sin2theta_square, Enu):
        self.MC_gen.dm2 = dm2
        self.MC_gen.sin2theta_square = sin2theta_square
        
        _, bl = self.MC_gen.generate_oscillate_Asimov_dataset(Enu)
        
        return bl
        