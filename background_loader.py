import pickle
import histlite
import numpy as np
import matplotlib.pyplot as plt
import sys

class BkgLoader:
    
    def __init__(self):
        self.bkg_filename = "/fs/ddn/sdf/group/nexo/users/miaoyu/Sterile_nu/sensitivity/data/backgrounds/background_ES_PRD.p"
    
        self.bkg_dataset = None
    
    
    def load_file(self):
        try:
            with open(self.bkg_filename, 'rb') as f:
                self.bkg_dataset = pickle.load(f)
        except Exception as e:
            print(e)
            sys.exit(-1)

    def spectrum(self, name, E):
        if name not in self.bkg_dataset.keys():
            print(f'{name} does not exist, must be in {self.bkg_dataset.keys()}..')
            return 0.
        else:
            h = self.bkg_dataset[name]
            return h.get_value(E)
    
    def plot(self, name):
        fig, ax = plt.subplots(figsize=(10, 8))
        
        if isinstance(name, list):
            for n in name:
                if n not in self.bkg_dataset.keys():
                    print(f'{n} does not exist, must be in {self.bkg_dataset.keys()}..')
                else:
                    h = self.bkg_dataset[n]
                    histlite.plot1d(ax, h, label=n)
        else:
            if name not in self.bkg_dataset.keys():
                print(f'{name} does not exist, must be in {self.bkg_dataset.keys()}..')
            else:
                h = self.bkg_dataset[name]
                histlite.plot1d(ax, h, label=name)
        ax.set_xlabel('Electron recoil kinetic energy [MeV]', fontsize=14)
        ax.set_ylabel('Counts per {int((h.bins[0][1] - h.bins[0][0])*1000)} keV', fontsize=14)
        ax.tick_params(labelsize=13)
        ax.legend(fontsize=14)
        ax.loglog()
        plt.tight_layout()
        plt.show()
        return fig

            

                
