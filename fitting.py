import os
import pickle
import re
import histlite as hl
import matplotlib.pyplot as plt
plt.style.use('fivethirtyeight')
#plt.style.use('ggplot')
import numpy as np

from iminuit import cost, Minuit


class fitting:
    def __init__(self) -> None:

        self.dataset_asimov = None
        self.dataset_asimov_filename = None
        self.dataset_asimov_flag = False
        
        self.dataset_toyMC = None
        self.dataset_toyMC_filename = None
        self.dataset_toyMC_flag = False
        
        self.data_dm2 = 0.0
        self.data_sin2 = 0.0
        
        self.PDF = None
        self.PDF_filename = None
        self.PDF_flag = False
        self.PDF_fitted = None

        self.fit_dm2 = 0.0
        self.fit_sin2 = 0.0
        
        self.asimov_fit = True
        self.fit_flag = False
        
        self.other_exp_files = []
        self.other_exp_sens = []

        # systematic budget:
        ## relative percentage
        self.sigma_init_flux = 0.02
        self.sigma_xsec = 0.11
        self.sigma_det_eff = 0.03

        self.alpha_init_flux = 0.0
        self.alpha_xsec = 0.0
        self.alpha_det_eff = 0.0
        

    ## Setters:
    def _set_data_dm2(self, val):
        self.data_dm2 = val
        
    def _set_data_sin2(self, val):
        self.data_sin2 = val
        
    def _set_fit_dm2(self, val):
        self.fit_dm2 = val
        
    def _set_fit_sin2(self, val):
        self.fit_sin2 = val

    def _set_dataset_asimov_filename(self, name):
        self.dataset_asimov_filename = name
    
    def _set_dataset_toyMC_filename(self, name):
        self.dataset_toyMC_filename = name
        
    def _set_PDF_filename(self, name):
        self.PDF_filename = name
            
    ## Getters:
    def _get_data_dm2(self):
        return self.data_dm2
    
    def _get_data_sin2(self):
        return self.data_sin2
    
    def _get_fit_dm2(self):
        return self.fit_dm2
    
    def _get_fit_sin2(self):
        return self.fit_sin2
    
        
    ### Loaders:
    def load_dataset_asimov(self):
        try:
            with open(self.dataset_asimov_filename, 'rb') as f:
                self.dataset_asimov = pickle.load(f)
            self.dataset_asimov_flag = True
        except Exception as e:
            print(f"Error occurs -> {e}")
            
    def load_dataset_toyMC(self):
        try:
            with open(self.dataset_toyMC_filename, 'rb') as f:
                self.dataset_toyMC = pickle.load(f)
            self.dataset_toyMC_flag = True
        except Exception as e:
            print(f"Error occurs -> {e}")

    def load_PDF(self):
        try:
            with open(self.PDF_filename, 'rb') as f:
                self.PDF = pickle.load(f)
            self.PDF_flag = True
            self.parse_parameters(self.PDF_filename)
        except Exception as e:
            print(f"Error occurs -> {e}")

    def reload(self):
        if not (self.dataset_asimov_filename is None):
            self.load_dataset_asimov()
        if not (self.dataset_toyMC_filename is None):
            self.load_dataset_toyMC()
        if not (self.PDF_filename is None):
            self.load_PDF()

        
    def parse_parameters(self, filename):
        dm2_pattern = r'dmsquare\s*([+-]?\d*\.\d+|\d+)'
        sin2_pattern = r'sinsquare\s*([+-]?\d*\.\d+|\d+)'
        dm2_match = re.search(dm2_pattern, filename)
        sin2_match = re.search(sin2_pattern, filename)
        dm2_value = dm2_match.group(1) if dm2_match else None
        sin2_value = sin2_match.group(1) if sin2_match else None

        dm2_value = float(dm2_value)
        sin2_value = float(sin2_value)
        self._set_fit_dm2(dm2_value)
        self._set_fit_sin2(sin2_value)
        
        label = r'$\Delta m^2=$' + f'{dm2_value:.5f}' + r' eV$^2, \sin^2(2\theta) = $' + f'{sin2_value:.5f}'
        return label
            
    def draw_histograms(self, data_asimov=True, data_toyMC=False, pdf=False, pdf_fitted=False):
        fig, ax = plt.subplots(figsize=(8, 6))
        if data_asimov:
            if not self.dataset_asimov_flag:
                self.load_dataset_asimov()
            h_data_asimov_error = hl.Hist(self.dataset_asimov.bins[0], \
                                          self.dataset_asimov.values, \
                                          errors = np.sqrt(self.dataset_asimov.values))
            data_asimov_linestyle = hl.LineStyle(line=True, errorbands=True, lw=2)
            hl.plot1d(ax, h_data_asimov_error, style=data_asimov_linestyle, label=f"Asimov dataset ({self.parse_parameters(self.dataset_asimov_filename)})")

        if data_toyMC:
            if not self.dataset_toyMC_flag:
                self.load_dataset_toyMC()
            h_data_toyMC_error = hl.Hist(self.dataset_toyMC.bins[0], \
                                         self.dataset_toyMC.values, \
                                         errors = np.sqrt(self.dataset_toyMC.values))
            data_toyMC_linestyle = hl.LineStyle(line=True, errorbands=True, marker='o', markers=True, lw=2)
            hl.plot1d(ax, h_data_toyMC_error, style=data_toyMC_linestyle, label=f"ToyMC dataset ({self.parse_parameters(self.dataset_asimov_filename)})")
            
        if pdf:
            if not self.PDF_flag:
                self.load_PDF()
            h_PDF_error = hl.Hist(self.PDF.bins[0], \
                                  self.PDF.values, \
                                  errors = np.sqrt(self.PDF.values**2 *(self.sigma_det_eff**2 + self.sigma_init_flux**2 + self.sigma_xsec**2)))
            pdf_linestyle = hl.LineStyle(line=True, errorbars=True, errorcaps=True, capsize=8, lw=2)
            hl.plot1d(ax, h_PDF_error, style=pdf_linestyle, label=f'PDF (before fitting, {self.parse_parameters(self.PDF_filename)})')
            
        if pdf_fitted:
            if not self.fit_flag:
                self.minimize_chi_square()
            if not self.fit_flag:
                print('Can not draw fitted PDF as fitting fails!')
            else:
                self.get_fitted_PDf()
                
            h_PDF_fitted_error = hl.Hist(self.PDF_fitted.bins[0], \
                                         self.PDF_fitted.values, \
                                         errors = np.sqrt(self.PDF_fitted.values**2 *(self.sigma_det_eff**2 + self.sigma_init_flux**2 + self.sigma_xsec**2)))
            pdf_fitted_linestyle = hl.LineStyle(line=True, errorbars=True, errorcaps=True, capsize=8, lw=2)
            hl.plot1d(ax, h_PDF_fitted_error, style=pdf_fitted_linestyle, label=f'PDF (after fitting, {self.parse_parameters(self.PDF_filename)})')
            
        ax.set_xlabel('Baseline [m]', fontsize=14)
        ax.set_ylabel('Signal count', fontsize=14)
        ax.legend(fontsize=13, bbox_to_anchor=(1.0, 1.3))
        ax.tick_params(labelsize=13)
        plt.tight_layout()
        plt.show()
        return fig

        
   
    def chi_square(self, alpha_xsec, alpha_init_flux, alpha_det_eff):
        if not self.PDF_flag:
            self.load_PDF()
        if self.asimov_fit :
            if not self.dataset_asimov_flag:
                self.load_dataset_asimov()
        else:
            if not self.dataset_toyMC_flag:
                self.load_dataset_toyMC()
                
        measured = self.dataset_asimov.values
        stat_err2 = measured
        sys_err2 = measured**2 * (self.sigma_det_eff**2 + self.sigma_init_flux**2 + self.sigma_xsec**2)
        tot_err = np.sqrt(stat_err2 + sys_err2)
        
        predicted = self.PDF.values * (1 + alpha_xsec + alpha_init_flux + alpha_det_eff)
        
        dchi2 = 0
        nbin = self.dataset_asimov.n_bins[0]
        for i in range(nbin):
            if tot_err[i] != 0:
                dchi2 += (measured[i] - predicted[i])**2 / tot_err[i]**2
                
        dchi2 += (alpha_det_eff/self.sigma_det_eff)**2 + (alpha_init_flux/self.sigma_init_flux)**2 + (alpha_xsec/self.sigma_xsec)**2

        return dchi2
        
        
    def minimize_chi_square(self, alpha_init_flux0 = 0.1, alpha_xsec0 = 0.02, alpha_det_eff0 = 0.02):
        fit_valid = False
        N_fit, N_fit_max = 0, 10
        while (not fit_valid) and (N_fit < N_fit_max):
            initials = {'alpha_xsec': np.random.normal(alpha_xsec0, 0.005), 'alpha_init_flux': np.random.normal(alpha_init_flux0, 0.005), 'alpha_det_eff': np.random.normal(alpha_det_eff0, 0.005)}
            m = Minuit(self.chi_square, **initials)
            m.limits['alpha_init_flux'] = (-1.0, 1.0)
            m.limits['alpha_det_eff']   = (-1.0, 1.0)
            m.limits['alpha_xsec']      = (-1.0, 1.0)
            
            m.migrad()
            fit_valid = m.valid
            N_fit += 1
            
        if not fit_valid:
            print(f'Fitting fails after {N_fit} trys, will suspend it for now.')
            return 0.0

        self.alpha_init_flux = m.values['alpha_init_flux']
        self.alpha_xsec = m.values['alpha_xsec']
        self.alpha_det_eff = m.values['alpha_det_eff']

        self.fit_flag = True
        return m.values, m.errors, m
        
        
    def get_fitted_PDf(self):
        values = self.PDF.values
        bins = self.PDF.bins[0]
        fitted_values = values * (1 + self.alpha_det_eff + self.alpha_xsec + self.alpha_xsec)
        self.PDF_fitted = hl.Hist(bins, fitted_values)

        

        
        