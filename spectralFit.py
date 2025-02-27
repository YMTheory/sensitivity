import numpy as np
import pandas as pd
import pickle
import histlite as hl
from iminuit import Minuit
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import sys
sys.path.append("/fs/ddn/sdf/group/nexo/users/miaoyu/Sterile_nu/sensitivity")
from oscillation import *

class spectralFit:
    
    def __init__(self) -> None:
        
        ## Fitting data and PDF
        self.asimov_dataset_filename    = None
        self.asimov_dataset             = None
        self.asimov_dataset_loadflag    = False
        self.asimov_dataset_projected   = None
        
        self.signal_PDF_filename        = None
        self.signal_PDF                 = None
        self.signal_PDF_projected       = None
        self.signal_PDF_loadflag        = False
        self.background_PDF_filename    = None
        self.background_PDF             = None
        self.background_PDF_projected   = None
        self.background_PDF_loadflag    = False
        self.PDF                        = None
        self.PDF_projected              = None

        ## Fitting parameters        
        self.dm_square_fit              = 0.0
        self.sin2theta_square_fit       = 0.0
        self.dm_square_data             = 0.0
        self.sin2theta_square_data      = 0.0

        ## Fitting channel
        self.channel                    = 'CC'
        self.fit_dimension              = 2

        ## Systematic uncertainties
        self.sigma_flux                 = 0.02
        self.sigma_xsec                 = 0.13
        self.sigma_efficiency           = 0.01
        self.sigma_background           = 0.005
        self.alpha_flux                 = 0
        self.alpha_xsec                 = 0
        self.alpha_efficiency           = 0
        self.alpha_background           = 0

        ## Fitting-related:
        self.all_parameters             = ['alpha_flux', 'alpha_xsec', 'alpha_efficiency', 'alpha_background']
        self.fixed_parameters           = []
        self.minimizer                  = None
        self.fitted_values              = None
        self.fitted_errors              = None
        self.fval                       = 0

        self.Ethr                       = 0.5
        self.Ethr_index                 = 0
        
        ## Sensitivity files
        self.other_exp_filename         = []
        self.other_exp_sensitivity      = []
        self.fitting_filename           = None
        self.fitting_results            = None
        self.fitting_loadflag           = False
        
    
    ## Files and loading
    def _set_asimov_dataset_filename(self, filename):
        self.asimov_dataset_filename = filename

    
    def _load_asimov_dataset(self):
        try:
            with open(self.asimov_dataset_filename, 'rb') as f:
                self.asimov_dataset = pickle.load(f)
            self.asimov_dataset_loadflag = True
            self.asimov_dataset_projected = self.asimov_dataset.project(axes=[0])
        except Exception as e:
            print(f'ERROR occurs at load_asimov_dataset -> {e}')

    def _set_signal_PDF_filename(self, filename):
        self.signal_PDF_filename = filename
        
    def _load_signal_PDF(self):
        try:
            with open(self.signal_PDF_filename, 'rb') as f:
                self.signal_PDF = pickle.load(f)
            self.signal_PDF_loadflag = True
            self.signal_PDF_projected = self.signal_PDF.project(axes=[0])
        except Exception as e:
            print(f'ERROR occurs load_signal_PDF -> {e}')

    def _set_background_PDF_filename(self, filename):
        self.background_PDF_filename = filename
        
    def _load_background_PDF(self):
        try:
            with open(self.background_PDF_filename, 'rb') as f:
                self.background_PDF = pickle.load(f)
            self.background_PDF_loadflag = True
            self.background_PDF_projected = self.background_PDF.project(axes=[0])
        except Exception as e:
            print(f'ERROR occurs load_background_PDF -> {e}')

    def _get_PDF(self, weight=False):
        if (not self.asimov_dataset_loadflag):
            self._load_asimov_dataset()
        if not self.signal_PDF_loadflag:
            self._load_signal_PDF()
        if self.background_PDF_filename is None:
            # In this case (CC channel), no background, PDF = signal PDF
            self.PDF = self.signal_PDF 
            if weight:
                self.PDF = hl.Hist(self.signal_PDF.bins, self.signal_PDF.values*(1+self.alpha_flux+self.alpha_xsec+self.alpha_efficiency))
        else:
            if not self.background_PDF_loadflag:
                self._load_background_PDF()
            if weight:
                self.PDF = hl.Hist(self.signal_PDF.bins, self.signal_PDF.values*(1+self.alpha_flux+self.alpha_xsec+self.alpha_efficiency) + self.background_PDF.values * (1+self.alpha_background))
            else:
                self.PDF = hl.Hist(self.signal_PDF.bins, self.signal_PDF.values + self.background_PDF.values)

    def _set_flux_uncertainty(self, sigma):
        self.sigma_flux = sigma
    
    def _set_cross_section_uncertainty(self, sigma):
        self.sigma_xsec = sigma
    
    def _set_efficiency_uncertainty(self, sigma):
        self.sigma_efficiency = sigma
    
    def _set_background_rate_uncertainty(self, sigma):
        self.sigma_background = sigma

    def _get_flux_uncertainty(self):
        return self.sigma_flux
    
    def _get_cross_section_uncertainty(self):
        return self.sigma_xsec
    
    def _get_efficiency_uncertainty(self):
        return self.sigma_efficiency
    
    def _get_background_rate_uncertainty(self):
        return self.sigma_background
        

    def _add_other_experiment_filename(self, filename):
        self.other_exp_filename.append( filename )

    def _load_other_experiment_sensitivity(self):
        for file in self.other_exp_filename:
            try:
                sens = np.loadtxt(file)
                self.other_exp_sensitivity.append( sens )
            except Exception as e:
                print(f'ERROR occurs load_other_experiment_sensitivity -> {e}')
        
    def _set_fitting_filename(self, filename):
        self.fitting_filename = filename
        
    def _load_fitting_results(self):
        try:
            self.fitting_results = pd.read_csv(self.fitting_filename)
            self.fitting_loadflag = True
        except Exception as e:
            print(f'ERROR occurs load_fitting_results -> {e}')
          
    def _load_allPDFs(self):
        self._load_asimov_dataset()
        self._load_signal_PDF()
        self._load_background_PDF()
        self._get_PDF()
    
    def _set_channel(self, cha):
        self.channel = cha

    def _set_fit_dimension(self, d):
        self.fit_dimension = d

    def _set_energy_threshold(self, E):
        self.Ethr = E
    
    def _set_dmsquare_fit(self, dm_square):
        self.dm_square_fit = dm_square
        
    def _set_sin2theta_square_fit(self, sin2theta_square):
        self.sin2theta_square_fit = sin2theta_square
                    
    def _set_dmsquare_data(self, dm_square):
        self.dm_square_data = dm_square
        
    def _set_sin2theta_square_data(self, sin2theta_square):
        self.sin2theta_square_data = sin2theta_square

    def _set_fixed_parameters(self, params):
        self.fixed_parameters = params

    def _set_fit_channel(self, cha):
        self.channel = cha

    def energy_cut(self):
        if self.channel == 'ES':
            E_idx = self.asimov_dataset.index(self.Ethr, axis=1)
            self.Ethr_index = E_idx
        
    def chi_square(self, alpha_flux, alpha_xsec, alpha_efficiency, alpha_background):

        if self.channel == 'ES':
            if self.fit_dimension == 2:
                measured = self.asimov_dataset.values[:, self.Ethr_index:-1]
                predicted_signal0 = self.signal_PDF.values[:, self.Ethr_index:-1]
            elif self.fit_dimension == 1:
                measured = self.asimov_dataset_projected.values[self.Ethr_index:-1]
                predicted_signal0 = self.signal_PDF_projected.values[self.Ethr_index:-1]
        elif self.channel == 'CC':
            measured = self.asimov_dataset.values
            predicted_signal0 = self.signal_PDF.values
        if self.channel == "ES":
            if self.fit_dimension == 2:
                predicted_background0 = self.background_PDF.values[:, self.Ethr_index:-1]
                predicted = predicted_signal0 * (1+alpha_flux+alpha_xsec+alpha_efficiency) + predicted_background0*(1+alpha_background)
            elif self.fit_dimension == 1:
                predicted_background0 = self.background_PDF_projected.values[self.Ethr_index:-1]
                predicted = predicted_signal0 * (1+alpha_flux+alpha_xsec+alpha_efficiency) + predicted_background0*(1+alpha_background)
        elif self.channel == 'CC':
            predicted = predicted_signal0 * (1+alpha_flux+alpha_xsec+alpha_efficiency)

        dchi2 = 0.
        stat_shape_err2 = measured
        dchi2 = np.where(stat_shape_err2>0, (measured-predicted)**2/stat_shape_err2, 0)
        dchi2 = np.sum(dchi2)
        ### pull terms
        dchi2 += alpha_flux**2/self.sigma_flux**2
        dchi2 += alpha_xsec**2/self.sigma_xsec**2
        dchi2 += alpha_efficiency**2/self.sigma_efficiency**2
        dchi2 += alpha_background**2/self.sigma_background**2
        return dchi2


    def minimize_chi_square(self, alpha_flux0=0, alpha_xsec0=0, alpha_efficiency0=0, alpha_background0=0):
        is_fit_valid = False
        n_fit, n_fit_max0 = 0, 10
        
        while (not is_fit_valid) and (n_fit < n_fit_max0):
            # initial values of fitting parameters
            initials = { \
                        'alpha_flux' : np.random.normal(alpha_flux0, 0.005), \
                        'alpha_xsec' : np.random.normal(alpha_xsec0, 0.005), \
                        'alpha_efficiency' : np.random.normal(alpha_efficiency0, 0.005), \
                        'alpha_background' : np.random.normal(alpha_background0, 0.005) \
                        }
        
            for k in self.fixed_parameters:
                initials[k]     = 0.0

            m = Minuit(self.chi_square, **initials)
            for k in self.all_parameters:
                if k in self.fixed_parameters:
                    m.fixed[k]      = True
                else:
                    m.limits[k] = (-10.0, 10.0)
        
            
            m.migrad()
            is_fit_valid = m.valid
            n_fit += 1

        if not is_fit_valid:
            print(f'Fitting fails after {n_fit} trys, will suspend it for now.')
    
        self.alpha_flux         = m.values['alpha_flux']
        self.alpha_xsec         = m.values['alpha_xsec']
        self.alpha_efficiency   = m.values['alpha_efficiency']
        self.alpha_background   = m.values['alpha_background']

        self.minimizer          = m
        self.fval               = m.fval
        self.fitted_values      = m.values
        self.fitted_errors      = m.errors
    

    def delta_chisquare_decomposition(self):
        if self.channel == 'ES':
            if self.fit_dimension == 2:
                measured = self.asimov_dataset.values[:, self.Ethr_index:-1]
                predicted_signal0 = self.signal_PDF.values[:, self.Ethr_index:-1]
            elif self.fit_dimension == 1:
                measured = self.asimov_dataset_projected.values[self.Ethr_index:-1]
                predicted_signal0 = self.signal_PDF_projected.values[self.Ethr_index:-1]
        elif self.channel == 'CC':
            measured = self.asimov_dataset.values
            predicted_signal0 = self.signal_PDF.values
        if self.channel == "ES":
            if self.fit_dimension == 2:
                predicted_background0 = self.background_PDF.values[:, self.Ethr_index:-1]
                predicted = predicted_signal0 * (1+self.alpha_flux+self.alpha_xsec+self.alpha_efficiency) + predicted_background0*(1+self.alpha_background)
            elif self.fit_dimension == 1:
                predicted_background0 = self.background_PDF_projected.values[self.Ethr_index:-1]
                predicted = predicted_signal0 * (1+self.alpha_flux+self.alpha_xsec+self.alpha_efficiency) + predicted_background0*(1+self.alpha_background)
        elif self.channel == 'CC':
            predicted = predicted_signal0 * (1+self.alpha_flux+self.alpha_xsec+self.alpha_efficiency)

        stat_shape_err2 = measured
        if self.channel == 'ES' and self.fit_dimension == 2:
            dchi2 = np.zeros(measured.shape)
            for i in range(measured.shape[0]):
                for j in range(measured.shape[1]):
                    if stat_shape_err2[i, j] != 0:
                        dchi2[i, j] = (measured[i,j]-predicted[i,j])**2 / stat_shape_err2[i,j]
        else:
            dchi2 = np.zeros(len(measured))
            for i in range(len(measured)):
                if stat_shape_err2[i] != 0:
                    dchi2[i] = (measured[i]-predicted[i])**2 / stat_shape_err2[i]
            
        penalty_flux = self.alpha_flux**2/self.sigma_flux**2
        penalty_xsec = self.alpha_xsec**2/self.sigma_xsec**2
        penalty_efficiency = self.alpha_efficiency**2/self.sigma_efficiency**2
        penalty_background = self.alpha_background**2/self.sigma_background**2
        return dchi2, penalty_flux, penalty_xsec, penalty_efficiency, penalty_background
        
    

    def parse_oscillation_parameters(self, dm_square, sin2theta_square):
        label = r'$\Delta m^2=$' + f'{dm_square:.5f}' + r' eV$^2, \sin^2(2\theta) = $' + f'{sin2theta_square:.5f}'
        return label
        

    def plot_fitting_histograms(self, data=True, pdf=True, fitted_pdf=False):
        fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(8, 14))
        if data:
            ## Baseline
            asimov_dataset_projection0 = self.asimov_dataset.project(axes=[0])
            h_data_asimov = hl.Hist(asimov_dataset_projection0.bins, \
                                    asimov_dataset_projection0.values, \
                                    errors = np.sqrt(asimov_dataset_projection0.values))
            h_data_asimov_linestyle = hl.LineStyle(line=True, errorbar=True, lw=2)
            hl.plot1d(ax0, h_data_asimov, style=h_data_asimov_linestyle, label=f'Asimov data {self.parse_oscillation_parameters(self.dm_square_data, self.sin2theta_square_data)}')

            ## Energy
            asimov_dataset_projection1 = self.asimov_dataset.project(axes=[1])
            h_data_asimov = hl.Hist(asimov_dataset_projection1.bins, \
                                    asimov_dataset_projection1.values, \
                                    errors = np.sqrt(asimov_dataset_projection1.values))
            h_data_asimov_linestyle = hl.LineStyle(line=True, errorbar=True, lw=2)
            hl.plot1d(ax1, h_data_asimov, style=h_data_asimov_linestyle, label=f'Asimov data {self.parse_oscillation_parameters(self.dm_square_data, self.sin2theta_square_data)}')
            
        if pdf:
            self._get_PDF(weight=False)
            pdf_projection0 = self.PDF.project(axes=[0])
            h_pdf = hl.Hist(pdf_projection0.bins, \
                                    pdf_projection0.values, \
                                    errors = np.sqrt(pdf_projection0.values))
            h_pdf_linestyle = hl.LineStyle(line=True, errorbar=True, lw=2, capsize=8)
            hl.plot1d(ax0, h_pdf, style=h_pdf_linestyle, label=f'Fit PDF (before) {self.parse_oscillation_parameters(self.dm_square_fit, self.sin2theta_square_fit)}')


            pdf_projection1 = self.PDF.project(axes=[1])
            h_pdf = hl.Hist(pdf_projection1.bins, \
                                    pdf_projection1.values, \
                                    errors = np.sqrt(pdf_projection1.values))
            h_pdf_linestyle = hl.LineStyle(line=True, errorbar=True, lw=2, capsize=8)
            hl.plot1d(ax1, h_pdf, style=h_pdf_linestyle, label=f'Fit PDF (before) {self.parse_oscillation_parameters(self.dm_square_fit, self.sin2theta_square_fit)}')

        if fitted_pdf:
            #self.minimize_chi_square()
            self._get_PDF(weight=True)
            fitted_pdf_projection0 = self.PDF.project(axes=[0])
            h_fitted_pdf = hl.Hist(fitted_pdf_projection0.bins, \
                                    fitted_pdf_projection0.values, \
                                    errors = np.sqrt(fitted_pdf_projection0.values))
            h_fitted_pdf_linestyle = hl.LineStyle(line=True, errorbar=True, lw=2, capsize=8)
            hl.plot1d(ax0, h_fitted_pdf, style=h_fitted_pdf_linestyle, label=f'Fit PDF (after) {self.parse_oscillation_parameters(self.dm_square_fit, self.sin2theta_square_fit)}')

            fitted_pdf_projection1 = self.PDF.project(axes=[1])
            h_fitted_pdf = hl.Hist(fitted_pdf_projection1.bins, \
                                    fitted_pdf_projection1.values, \
                                    errors = np.sqrt(fitted_pdf_projection1.values))
            h_fitted_pdf_linestyle = hl.LineStyle(line=True, errorbar=True, lw=2, capsize=8)
            hl.plot1d(ax1, h_fitted_pdf, style=h_fitted_pdf_linestyle, label=f'Fit PDF (after) {self.parse_oscillation_parameters(self.dm_square_fit, self.sin2theta_square_fit)}')
        
        ax0.set_xlabel('Baseline [m]', fontsize=14)
        ax0.set_ylabel('Signal count', fontsize=14)
        ax0.legend(fontsize=13, bbox_to_anchor=(1.0, 1.3))
        ax0.tick_params(labelsize=13)
        ax1.set_xlabel('Electron energy [MeV]', fontsize=14)
        ax1.set_ylabel('Signal count', fontsize=14)
        ax1.tick_params(labelsize=13)
        plt.tight_layout()
        plt.show()
        return fig
            
        
        