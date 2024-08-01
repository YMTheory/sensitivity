import numpy as np

class detector:
    def __init__(self, name):
        self.name = name

        self.livetime = 40  # unit: day
        
        # Detection efficiency
        self.efficiency = 1.0
        
        # Detection resolution
        self.energy_resolution = 0.01
        self.spatial_resolution = 0.01 # unit: m

        self.comb_res_a = 0.317 # sqrt(keV)
        self.comb_res_b = 0.015
        
        self.q_res_a = 0.75 # sqrt(keV)
        self.q_res_b = 0.035
        
        # position of detector center
        self.position = ( 0, 0, 0)
        
        # Detector-source distance
        self.baseline = 1# unit: m
        # which can also be calculated by detector center position and source position

        # Running time
        self.run_time = 100 # unit: day

        # Detector mass and geometry
        self.FV_mass = 3281 # kg
        self.Xe136_mass = 136 # g / mol
        self.NA = 6.022e23 # mol^-1
        self.N_Xe136 = self.FV_mass / self.Xe136_mass * self.NA 
        self.height = 1.183 # m
        self.radius = 0.5665 # m
        self.volume = np.pi * self.radius**2 * self.height # m^3

        #TODO: add Xe enrichment?
        
        
    # Two empirical energy resolution functions from Brain's solar nu paper
    def empirical_combined_energy_resolution(self, E):
        E = E * 1000. # input unit: MeV, convert to keV
        return E * (self.comb_res_a / np.sqrt(E) + self.comb_res_b) / 1000. # output unit: MeV
    
    def empirical_charge_energy_resolution(self, E):
        E = E * 1000. # input unit: MeV, convert to keV
        return E*(self.q_res_a / np.sqrt(E) + self.q_res_b)/1000. # output unit: MeV
        