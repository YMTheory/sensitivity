import os
import pickle
import re
import pandas as pd
import histlite as hl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
plt.style.use('fivethirtyeight')
#plt.style.use('ggplot')
import numpy as np

from iminuit import cost, Minuit

import sys
sys.path.append("/p/lustre1/yu47/Sterile_Neutrino/sensitivity/")
from oscillation import *

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
        self.alpha  = 0.0
        
        ## Besides do fitting, another functionailty will be loading the pre-fitted results and plotting
        self.fitting_filename = None
        self.fitting_results = None
        self.fitting_load_flag = False
        
        self.fitting_mode = 'shape' # options in: [rate, shape]
        self.constraint   = 'constraint' # options in: [free, constraint, fixed] 

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
            
    def _set_fitting_mode(self, mode):
        self.fitting_mode = mode

    def _set_constraint(self, mode):
        self.constraint = mode
            
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
        
    def _get_fitting_mode(self):
        return self.fitting_mode
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
            #self.fitting_results = np.loadtxt(self.fitting_filename)
            self.fitting_results = pd.read_csv(self.fitting_filename)
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
            
    def draw_histograms(self, data_asimov=True, data_toyMC=False, pdf=False, pdf_fitted=False,):
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

    def draw_fitting_results(self, variable, contour=False, others=False, N_dm2=100, dm2_start=-2, dm2_stop=1,  N_sin2=100, sin2_start=-2, sin2_stop=0, drawLength=False, Enu=0.75):
        if not self.fitting_load_flag:
            self.load_fitting_results()

            
        arr = self.fitting_results[variable].to_numpy()
        if variable == 'alpha_xsec':
            lb = r"$\alpha_{xsec}$"
        elif variable == 'alpha_xsec_err':
            lb = r"$\epsilon_{\alpha_{xsec}}$"
        elif variable == 'alpha_init_flux':
            lb = r"$\alpha_{init-flux}$"
        elif variable == 'alpha_init_flux_err':
            lb = r"$\epsilon_{\alpha_{init_flux}}$"
        elif variable == 'alpha_det_eff':
            lb = r"$\alpha_{eff}$"
        elif variable == 'alpha_det_eff_err':
            lb = r"$\epsilon_{\alpha_{det_eff}}$"
        elif variable == 'delta_chisquare':
            lb = r"$\Delta\chi^2$"
        elif variable == 'total_rate_scaling':
            if self.constraint == 'constraint':
                lb = r"$1+\alpha_{xsec}+\alpha_{eff}+\alpha_{init-flux}$"
            elif self.constraint == "free":
                lb = r"$1+\alpha$"
        else:
            print(f'No such variable, can only be in {self.fitting_results.columns}')
        
        arr = np.reshape(arr, (N_dm2, N_sin2))
        
        dm2_edges   = np.logspace(dm2_start, dm2_stop, N_dm2+1)
        sin2_edges  = np.logspace(sin2_start, sin2_stop, N_sin2+1)
        dm2_cents   = np.logspace(dm2_start, dm2_stop, N_dm2)
        sin2_cents  = np.logspace(sin2_start, sin2_stop, N_sin2)

        if drawLength:
            length_edges = oscillation_length(dm2_edges, Enu)
            length_cents = oscillation_length(dm2_cents, Enu)

        fig, ax = plt.subplots(figsize=(9, 6))
        if variable == 'delta_chisquare':
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
        ax.tick_params(labelsize=13)
        ax.set_xlim(sin2_edges[0], sin2_edges[-1])
        ax.set_ylim(dm2_edges[0], dm2_edges[-1])
        if drawLength:
            ax.set_ylabel('oscillation length [m]', fontsize=14)
            #locs, _ = plt.yticks()
            #ax.set_yticks(locs, [f'{elem:.4f}' for elem in oscillation_length(locs, Enu)])
            locs = np.array([1e-2, 1e-1, 1, 10])
            length_ticks = oscillation_length(locs, Enu)
            ax.set_yticks(locs, [f'{el:.3f}' for el in length_ticks])
        else:
            ax.set_ylabel(r'$\Delta m^2\ [\mathrm{eV}^2]$', fontsize=14)
        #fig.savefig("../plots/dchi2_MCalg_1cmSmear_3cmBinWidth_14cmSource.png")

        plt.tight_layout()
        plt.show()
        if contour:
            contour_points = []
            for coll in CS.collections:
                for path in coll.get_paths():
                    vertices = path.vertices
                    x_contour, y_contour = vertices[:, 0], vertices[:, 1]
                    contour_points.append((x_contour, y_contour))
            return fig, contour_points
        return fig 
            
    def decompose_chi_square(self):
        if not self.fit_flag:
            self.minimize_chi_square()
        if self.asimov_fit:
            if not self.dataset_asimov_flag:
                self.load_dataset_asimov()
            measured = self.dataset_asimov.values
        else:
            if not self.dataset_toyMC_flag:
                self.load_dataset_toyMC()
            measured = self.dataset_toyMC.values
        
        predicted = self.PDF.values
        predicted = predicted * (1 + self.alpha_det_eff + self.alpha_init_flux + self.alpha_xsec + self.alpha)

        dchi2_bins = np.zeros(len(predicted))
        if self.fitting_mode == 'rate':
            return dchi2_bins
        
        if self.fitting_mode == 'shape':
            stat_shape_err2 = measured
            sys_shape_err2 = 0.
            #if self.constraint == 'constraint':
            #    sys_shape_err2 = measured**2 * (self.sigma_det_eff**2 + self.sigma_init_flux**2 + self.sigma_xsec**2)
            tot_shape_err2 = stat_shape_err2 + sys_shape_err2 

            nbin = self.dataset_asimov.n_bins[0]
            for i in range(nbin):
                if tot_shape_err2[i] != 0:
                    tmp_dchi2 = (measured[i] - predicted[i])**2 / tot_shape_err2[i]
                    dchi2_bins[i] = tmp_dchi2
            
            return dchi2_bins
    
    def chi_square(self, alpha_xsec, alpha_init_flux, alpha_det_eff, alpha):
        # Check PDF and dataset loading first.
        if not self.PDF_flag:
            self.load_PDF()
        if self.asimov_fit:
            if not self.dataset_asimov_flag:
                self.load_dataset_asimov()
            measured = self.dataset_asimov.values
        else:
            if not self.dataset_toyMC_flag:
                self.load_dataset_toyMC()
            measured = self.dataset_toyMC.values

        # Re-prediction:
        ## TO BE NOTED: one need to fixed the unused nuisance parameters all as 0. 
        predicted = self.PDF.values
        predicted = predicted * (1 + alpha_det_eff + alpha_init_flux + alpha_xsec + alpha) 

        # calculating delta_chisquare
        dchi2 = 0
        if self.fitting_mode == 'rate':
            measured_rate = np.sum(measured)
            predicted_rate = np.sum(predicted)
            stat_rate_err2 = measured_rate
            sys_rate_err2 = 0.
            #if self.constraint:
            #    sys_rate_err2 = measured_rate**2 * (self.sigma_det_eff**2 + self.sigma_init_flux**2 + self.sigma_xsec**2)
            tot_rate_err2 = stat_rate_err2 + sys_rate_err2

            dchi2 += (measured_rate - predicted_rate)**2 / tot_rate_err2
            dchi2 += (alpha_det_eff)**2/self.sigma_det_eff**2 + alpha_init_flux**2/self.sigma_init_flux**2 + alpha_xsec**2/self.sigma_xsec**2
            return dchi2

        elif self.fitting_mode == 'shape':
            stat_shape_err2 = measured
            sys_shape_err2 = 0.
            #if self.constraint == 'constraint':
            #    sys_shape_err2 = measured**2 * (self.sigma_det_eff**2 + self.sigma_init_flux**2 + self.sigma_xsec**2)
            tot_shape_err2 = stat_shape_err2 + sys_shape_err2 

            nbin = self.dataset_asimov.n_bins[0]
            for i in range(nbin):
                if tot_shape_err2[i] != 0:
                    dchi2 += (measured[i] - predicted[i])**2 / tot_shape_err2[i]
            dchi2 += (alpha_det_eff)**2/self.sigma_det_eff**2 + alpha_init_flux**2/self.sigma_init_flux**2 + alpha_xsec**2/self.sigma_xsec**2
            return dchi2 
   
    #def chi_square(self, alpha_xsec, alpha_init_flux, alpha_det_eff):
    #    if not self.PDF_flag:
    #        self.load_PDF()
    #    if self.asimov_fit :
    #        if not self.dataset_asimov_flag:
    #            self.load_dataset_asimov()
    #    else:
    #        if not self.dataset_toyMC_flag:
    #            self.load_dataset_toyMC()
    #            
    #    measured = self.dataset_asimov.values
    #    stat_err2 = measured
    #    sys_err2 = measured**2 * (self.sigma_det_eff**2 + self.sigma_init_flux**2 + self.sigma_xsec**2)
    #    tot_err = np.sqrt(stat_err2 + sys_err2)
    #    
    #    predicted = self.PDF.values * (1 + alpha_xsec + alpha_init_flux + alpha_det_eff)
    #    
    #    dchi2 = 0
    #    nbin = self.dataset_asimov.n_bins[0]
    #    for i in range(nbin):
    #        if tot_err[i] != 0:
    #            dchi2 += (measured[i] - predicted[i])**2 / tot_err[i]**2
    #            
    #    dchi2 += (alpha_det_eff/self.sigma_det_eff)**2 + (alpha_init_flux/self.sigma_init_flux)**2 + (alpha_xsec/self.sigma_xsec)**2

    #    return dchi2
        
        
    def minimize_chi_square(self, alpha_init_flux0 = 0.1, alpha_xsec0 = 0.02, alpha_det_eff0 = 0.02, alpha0=0.0):
        
        fit_valid = False
        N_fit, N_fit_max = 0, 10
        if self.constraint != 'fixed':
            while (not fit_valid) and (N_fit < N_fit_max):
                if self.constraint == 'constraint':
                    initials = {'alpha_xsec': np.random.normal(alpha_xsec0, 0.005), 'alpha_init_flux': np.random.normal(alpha_init_flux0, 0.005), 'alpha_det_eff': np.random.normal(alpha_det_eff0, 0.005), 'alpha': 0}
                    m = Minuit(self.chi_square, **initials)

                    m.limits['alpha_init_flux'] = (-10.0, 10.0)
                    m.limits['alpha_det_eff']   = (-10.0, 10.0)
                    m.limits['alpha_xsec']      = (-10.0, 10.0)
                    m.fixed['alpha']            = True
            
                elif self.constraint == 'free':
                    initials = {'alpha_xsec': 0.0, 'alpha_init_flux': 0., 'alpha_det_eff': 0.,  'alpha': np.random.normal(alpha0, 0.005)}
                    m = Minuit(self.chi_square, **initials)
                    m.fixed['alpha_init_flux']  = True
                    m.fixed['alpha_det_eff']    = True
                    m.fixed['alpha_xsec']       = True
                    m.limits['alpha']           = (-10.0, 10.0)
            
                m.migrad()
                fit_valid = m.valid
                N_fit += 1
            
            if not fit_valid:
                print(f'Fitting fails after {N_fit} trys, will suspend it for now.')

            ## Allocate nuisance parameter fitting results
            self.alpha_init_flux = m.values['alpha_init_flux']
            self.alpha_xsec = m.values['alpha_xsec']
            self.alpha_det_eff = m.values['alpha_det_eff']
            self.alpha = m.values['alpha']
            self.minimizer = m
            self.fit_flag = True
            return m, m.values, m.errors

        elif self.constraint == 'fixed':
            # If everything is fixed, we do not need to minimize dchi_square. We only need to calculate it.
            dchi2 = self.chi_square(0.0, 0.0, 0.0, 0.0)
            self.alpha_det_eff = 0.0
            self.alpha_init_flux = 0.0
            self.alpha_xsec = 0.0
            self.alpha = 0.0
            return dchi2, [0, 0, 0, 0,], [0, 0, 0, 0,]
        
        
    def get_fitted_PDf(self):
        values = self.PDF.values
        bins = self.PDF.bins[0]
        fitted_values = values * (1 + self.alpha_det_eff + self.alpha_xsec + self.alpha_init_flux + self.alpha)
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
        return 1 + self.alpha_init_flux + self.alpha_xsec + self.alpha_det_eff + self.alpha

    def get_dataset_stats(self):
        return np.sum(self.dataset_asimov.values)

    def get_pdf_stats(self):
        return np.sum( self.PDF.values )

    def get_fitted_pdf_stats(self):
        return np.sum( self.PDF_fitted.values)