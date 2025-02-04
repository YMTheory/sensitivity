import sys
sys.path.append("/p/lustre1/yu47/Sterile_Neutrino/sensitivity")

import numpy as np
from iminuit import Minuit

import numpy as np
from fitting import fitting

class combining:
    def __init__(self) -> None:
        self.fitter_list = {}
        self.syst_err = {}
        self.nuisance = {}

        ## Fitting results
        self.minimizer = None
        self.fit_value = 0.

    def add_fitter(self, name, fitter):
        self.fitter_list[name] = fitter

    def add_syst_err(self, name, sigma):
        self.syst_err[name] = sigma

    def add_nuisance(self, name, alpha):
        self.nuisance[name] = alpha
        
    def fluctuate_nuisance(self):
        for k, v in self.nuisance.items():
            self.nuisance[k] = np.random.uniform(v, 0.005)
        
    def load(self):
        for f in self.fitter_list.values():
            f.reload()

    def chi_square(self, alpha_init_flux=0.0, alpha_xsec_CC=0.0, alpha_xsec_ES=0.0, alpha_det_eff=0.0):
    #def chi_square(self, **alphas):
        dchi2 = 0.
        for k, f in self.fitter_list.items():
            if 'ES' in k:
                dchi2 += f.chi_square(alpha_xsec_ES, alpha_init_flux, alpha_det_eff, 0.)
            elif 'CC' in k:
                dchi2 += f.chi_square(alpha_xsec_CC, alpha_init_flux, alpha_det_eff, 0.)
            else:
                print(f'Unknow process in {k}')
        dchi2 += alpha_init_flux**2 / self.syst_err['init_flux']**2
        dchi2 += alpha_xsec_CC**2 / self.syst_err['xsec_CC']**2
        dchi2 += alpha_xsec_ES**2 / self.syst_err['xsec_ES']**2
        dchi2 += alpha_det_eff**2 / self.syst_err['det_eff']**2
        return dchi2

    def minimize_combined_chi_square(self):
        fit_valid = False
        N_fit, N_fit_max = 0, 10
        
        while (not fit_valid) and (N_fit < N_fit_max):
            self.fluctuate_nuisance()
            initials = self.nuisance
            m = Minuit(self.chi_square, **initials)
            for k in self.nuisance.keys():
                m.limits[k] = (-10., 10.)
            m.migrad()
            fit_valid = m.valid
            N_fit += 1
        
        if not fit_valid:
            print(f"Fitting fails after {N_fit} trys, will suspend it now.")
        
        for k in self.nuisance.keys():
            self.nuisance[k] = m.values[k]
        self.minimizer = m
        self.fit_value = m.fval
        
        return m, m.values, m.errors        

        
'''
cbf = combining()        

f = fitting()
f._set_fitting_mode('shape')
f._set_constraint('fixed')
filename_data = '/p/lustre1/yu47/Sterile_Neutrino/jobs/MC/dmsquare0.00000eV2/MC_dmsquare0.00000sinsquare0.00000_source14cm_smear0.01m_pre_smeared_dist10cm_scaled.p'
f._set_dataset_asimov_filename(filename_data)
filename_fit = '/p/lustre1/yu47/Sterile_Neutrino/jobs/MC/dmsquare3.05386eV2/MC_dmsquare3.05386sinsquare0.86975_source4cm_smear0.01m_pre_smeared_scaled.p'
f._set_PDF_filename(filename_fit)
f.reload()

cbf.add_fitter('nEXO_10cm_CC', f)

f1 = fitting()
f1._set_fitting_mode('shape')
f1._set_constraint('fixed')
filename_data = '/p/lustre1/yu47/Sterile_Neutrino/jobs/MC/dmsquare0.00000eV2/MC_dmsquare0.00000sinsquare0.00000_source14cm_smear0.01m_pre_smeared_dist10cm_scaled_EStot.p'
f1._set_dataset_asimov_filename(filename_data)
filename_fit = '/p/lustre1/yu47/Sterile_Neutrino/jobs/MC/dmsquare3.05386eV2/MC_dmsquare3.05386sinsquare0.86975_source14cm_smear0.01m_pre_smeared_nEXO_dist10cm_ES.p'
f1._set_PDF_filename(filename_fit)
f1.reload()

cbf.add_fitter('nEXO_10cm_ES', f1)

cbf.add_syst_err('init_flux', 0.02)
cbf.add_syst_err('xsec_CC', 0.11)
cbf.add_syst_err('xsec_ES', 0.11)
cbf.add_syst_err('det_eff', 0.01)

cbf.add_nuisance('alpha_init_flux', 0.0)
cbf.add_nuisance('alpha_xsec_CC', 0.0)
cbf.add_nuisance('alpha_xsec_ES', 0.0)
cbf.add_nuisance('alpha_det_eff', 0.0)

cbf.load()
m, vals, errs = cbf.minimize_combined_chi_square()
print(m)
print(cbf.fit_value)
#cbf.chi_square(0, 0)


'''
            
     


