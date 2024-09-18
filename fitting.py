import os
import pickle
import re
import histlite as hl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
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
        self.minimizer = None
        
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
        
        ## Besides do fitting, another functionailty will be loading the pre-fitted results and plotting
        self.fitting_filename = None
        self.fitting_results = None
        self.fitting_load_flag = False
        

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

    def _set_fitting_filename(self, name):
        self.fitting_filename = name
        
    def _set_other_experimental_filenames(self, names):
        for name in names:
            self.other_exp_files.append(name)

    def _set_sigma_det_eff(self, val):
        self.sigma_det_eff = val

    def _set_sigma_xsec(self, val):
        self.sigma_xsec = val

    def _set_sigma_init_flux(self, val):
        self.sigma_init_flux = val
            
    ## Getters:
    def _get_data_dm2(self):
        return self.data_dm2
    
    def _get_data_sin2(self):
        return self.data_sin2
    
    def _get_fit_dm2(self):
        return self.fit_dm2
    
    def _get_fit_sin2(self):
        return self.fit_sin2

    def _get_sigma_det_eff(self):
        return self.sigma_det_eff

    def _get_sigma_init_flux(self):
        return self.sigma_init_flux
    
    def _get_sigma_xsec(self):
        return self.sigma_xsec
        
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
        # The fitting flag is set as fault here as new histograms has been loaded.
        self.fit_flag = False

    def load_other_experimental_sensitivities(self):
        for file in self.other_exp_files:
            if not os.path.exists(file):
                print(f"{file} does not exist !")
                continue
            else:
                arr = np.loadtxt(file)
            self.other_exp_sens.append( arr )

    def load_fitting_results(self):
        try:
            self.fitting_results = np.loadtxt(self.fitting_filename)
            self.fitting_load_flag = True
        except Exception as e:
            print(f"Error occurs -> {e}")

        
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

    def draw_fitting_results(self, variable, contour=False, others=False, N_dm2=100, dm2_start=-2, dm2_stop=1,  N_sin2=100, sin2_start=-2, sin2_stop=0):
        if not self.fitting_load_flag:
            self.load_fitting_results()
        if variable == 'alpha_xsec':
            arr = self.fitting_results[:, 0]
            lb = r"$\alpha_{xsec}$"
        elif variable == 'alpha_init_flux':
            arr = self.fitting_results[:, 1]
            lb = r"$\alpha_{init-flux}$"
        elif variable == 'alpha_det_eff':
            lb = r"$\alpha_{eff}$"
            arr = self.fitting_results[:, 2]
        elif variable == 'chi_square':
            arr = self.fitting_results[:, 3]
            lb = r"$\Delta\chi^2$"
        elif variable == 'total_rate_scaling':
            arr = self.fitting_results[:, 0] + self.fitting_results[:, 1] + self.fitting_results[:, 2] + 1
            lb = r"$1+\alpha_{xsec}+\alpha_{eff}+\alpha_{init-flux}$"
        else:
            print('No such variable, can only be in [alpha_xsec, alpha_init_flux, alpha_det_eff, chi_square, total_rate_scaling]')
        arr = np.reshape(arr, (N_dm2, N_sin2))
        
        
        dm2_edges = np.logspace(dm2_start, dm2_stop, N_dm2+1)
        sin2_edges = np.logspace(sin2_start, sin2_stop, N_sin2+1)
        dm2_cents = np.logspace(dm2_start, dm2_stop, N_dm2)
        sin2_cents = np.logspace(sin2_start, sin2_stop, N_sin2)

        fig, ax = plt.subplots(figsize=(9, 6))
        if variable == 'chi_square':
            im = ax.pcolormesh(sin2_edges, dm2_edges, arr, shading='flat', cmap='viridis', norm=colors.LogNorm(vmin=1e-5, vmax=np.max(arr)))
            if contour:
                X, Y = np.meshgrid(sin2_cents, dm2_cents)
                CS = ax.contour(X, Y, arr, levels=[4.605], cmap='Spectral',)
            if others:
                self.load_other_experimental_sensitivities()
                for data, name in zip(self.other_exp_sens, self.other_exp_files):
                    if 'reactor' in name:
                        ax.plot(data[:,0], data[:,1], '--', lw=3, color='darkviolet')
                    elif 'LZ' in name:
                        ax.plot(data[:,0], data[:,1], '-.', lw=3, color='blue')
        else:
            im = ax.pcolormesh(sin2_edges, dm2_edges, arr, shading='flat', cmap='viridis', )

        ax.loglog()
        cb = plt.colorbar(im, ax=ax,)
        cb.set_label(lb, fontsize=14, rotation=270, labelpad=25)
        cb.ax.tick_params(labelsize=13)

        ax.set_xlabel(r'$\sin^2(2\theta)$', fontsize=14)
        ax.set_ylabel(r'$\Delta m^2$', fontsize=14)
        ax.tick_params(labelsize=13)
        ax.set_xlim(sin2_edges[0], sin2_edges[-1])
        ax.set_ylim(dm2_edges[0], dm2_edges[-1])
        #fig.savefig("../plots/dchi2_MCalg_1cmSmear_3cmBinWidth_14cmSource.png")

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
        self.minimizer = m

        self.fit_flag = True
        return m.values, m.errors, m
        
        
    def get_fitted_PDf(self):
        values = self.PDF.values
        bins = self.PDF.bins[0]
        fitted_values = values * (1 + self.alpha_det_eff + self.alpha_xsec + self.alpha_xsec)
        self.PDF_fitted = hl.Hist(bins, fitted_values)

        
    def print_fitting(self):
        if not self.fit_flag:
            self.minimize_chi_square()
            if not self.fit_flag:
                print('The minimization failed :( ')
            else:
                print( self.minimizer )
        else:
            print( self.minimizer )

        
    def total_rate_scaling(self):
        return 1 + self.alpha_init_flux + self.alpha_xsec + self.alpha_det_eff