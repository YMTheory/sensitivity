import numpy as np
import pandas as pd
import scipy.integrate as integrate
import matplotlib.pyplot as plt
plt.style.use('fast')

import sys
sys.path.append('/p/lustre1/yu47/Sterile_Neutrino/sensitivity')
import xsection

class bkg:
    def __init__(self):
        self.dataset = {}

        self.time = 365 # unit: d
        self.mass = 1 # unit: ton
        
        self.load_flag = False

    def load_data(self, name):
        path = "/p/lustre1/yu47/Sterile_Neutrino/sensitivity/data/backgrounds/" 
        df = pd.read_csv(path+name+'.csv', sep=',')
        self.dataset[name] = df
        self.dataset[name]['energy'] = self.dataset[name]['energy'] / 1000.
        self.load_flag = True


    def plot_data(self):
        if not self.load_flag:
            self.load_data('2nbb')
            self.load_data('solar-pp')
            self.load_data('Kr85')
            self.load_data('Rn222')
        fig, ax = plt.subplots(figsize=(8, 6))
        for k, v in self.dataset.items():
            lb = k
            E = v['energy']
            C = v['events']
            ax.plot(E, C, label=lb)
        ax.set_xlabel('electron energy [keV]')
        ax.set_ylabel(r'events [keV$^{-1}\cdot$ ton $^{-1}\cdot$] yr $^{-1}$')
        ax.legend()
        ax.loglog()
        plt.tight_layout()
        plt.show()


    def background_counts(self, det, time, bkg, Ethr=0.0, Ehig=1000):
        if not self.load_flag:
            self.load_data('2nbb')
            self.load_data('solar-pp')
            self.load_data('Kr85')
            self.load_data('Rn222')
            
        if bkg not in self.dataset.keys():
            print(f'No {bkg} in our background list!!')
            return 0.0
            
        v = self.dataset[bkg]
        x, y = v['energy'], v['events']
        result = integrate.quad(lambda e: np.interp(e, x, y), Ethr, Ehig)
        itg = result[0] 
        n_events = itg * (det.FV_mass*1e-3) * (time/365.)
        #print(f'{lb} background count in {det.name} detector within {time} days is: {n_events:.3e}')
        return n_events
            

    def background_spectrum(self, det, time, bkg, E):
        if not self.load_flag:
            self.load_data('2nbb')
            self.load_data('solar-pp')
            self.load_data('Kr85')
            self.load_data('Rn222')
        
        if bkg not in self.dataset.keys() and bkg != 'total':
            print(f'No {bkg} in our background list!!')
            return 0.0

        elif bkg in self.dataset.keys():
            v = self.dataset[bkg] 
            Emax = np.max(v['energy'])
            x, y = v['energy'], v['events']*(det.FV_mass*1e-3) *(time/365.)

            res = np.interp(E, x, y)
            return np.where(E>Emax, 0, res)

        elif bkg == 'total':
            res = 0
            for v in self.dataset.values():
                Emax = np.max(v['energy'])
                x, y = v['energy'], v['events']*(det.FV_mass*1e-3) *(time/365.)
                tmp_res = np.where(E>Emax, 0, np.interp(E, x, y))
                res += tmp_res
            return res

                


    def generate_bkg(self, det, source):
        """
        The unit from the loading file is: events / keV / ton / yr
        """
        path = "/p/lustre1/yu47/Sterile_Neutrino/sensitivity/data/backgrounds/" 
        for k, v in self.dataset.items():
            filename = path + k + '_' + det.name + '.csv'

            factor = (det.FV_mass*1e-3) * (source.time / 365)
            




'''
obj = bkg()
obj.load_data("2nbb")
obj.load_data("Kr85")
obj.load_data("Rn222")
obj.load_data("solar-pp")


obj.plot_data()

'''