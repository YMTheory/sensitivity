import numpy as np

class neutrino_source:
    def __init__(self, name='Cr51', n_events=100000, energies=['0.75'], ratios=[1.0], time=100, activity=5e6, height=0.14, diameter=0.14):
        # one needs to specify the name of the source, total number of events, energies of the neutrinos and their branch ratios.
        print(f"A {name} hot neutrino source with nominal {n_events} in xx years with neutrino energies {energies} and branch ratios {ratios} is created.")
        
        self.name = name
        self.halflife = 27.7 # unit day
        if name == "Cr51":
            self.halflife = 27.7
        self.n_events = n_events
        self.energies = energies
        self.ratios = ratios
        self.time = time
        
        self.height = height # unit: m
        self.diameter = diameter # unit: m

        self.position = (0, 0, 0)

        self.activity = activity # unit: MCi, 1 Ci = 3.7e10 dacays per second


    def set_source_positin(self, x, y, z):
        self.position = (x, y, -1)
        
    def sample_vertice(self, nevent):
        # TODO: consider the geometry of the source, which will influence the baseline calculation.
        pass
    
    def flux_activity(self, t):
        return self.activity * np.exp(-t*(0.693/self.halflife))
    
    def flux_integral(self):
        lbd = 0.693 / self.halflife
        return -1/lbd*(self.flux_activity(self.time) - self.flux_activity(0.))        
        
    def flux_integral_nominal(self):
        lbd = 0.693 / 27.7
        return -1/lbd*(self.flux_activity(100) - self.flux_activity(0.))        
