import numpy as np

class neutrino_source:
    def __init__(self, name, n_events, energies, ratios):
        # one needs to specify the name of the source, total number of events, energies of the neutrinos and their branch ratios.
        print(f"A {name} hot neutrino source with nominal {n_events} in xx years with neutrino energies {energies} and branch ratios {ratios} is created.")
        
        self.name = name
        self.n_events = n_events
        self.energies = energies
        self.ratios = ratios
        
        self.height = 0.14 # unit: m
        self.diameter = 0.14 # unit: m

        self.position = (0, 0, 0)

        self.activity = 5e6 # unit: Ci, 1 Ci = 3.7e10 dacays per second


    def set_source_positin(self, x, y, z):
        self.position = (x, y, -1)
        
    def sample_vertice(self, nevent):
        # TODO: consider the geometry of the source, which will influence the baseline calculation.
        pass
    
