import os
import numpy as np
#from iminuit import cost, Minuit
import histlite as hl
import sys
import h5py as h5
import matplotlib.pyplot as plt

from iminuit import cost, Minuit

from MC_generator import MC_generator


class fitter:
    def __init__(self, source, detector) -> None:
        self.source = source
        self.detector = detector

        self.MC_gen = MC_generator(self.source, self.detector)

        self.scale_count_flag = False
        self.MC_gen.scale_counts()

        self.dataset = []
        self.PDF  = None
        self.data_filename = ''
        self.load_data_flag = False
        self.pdf_filename = ''
        self.load_pdf_flag = False

        self.other_exp_files = []    
        self.other_exp_sens = []    

        self.data_dm2 = 0.0
        self.data_sin2 = 0.0
        self.fit_dm2 = 0.0
        self.fit_sin2 = 0.0

        self.bin_width = 0.03

    ######################################################################## 
    ## Setters:
    def _set_data_dm2(self, val):
        self.data_dm2 = val
    
    def _set_data_sin2(self, val):
        self.data_sin2 = val
        
    def _set_fit_dm2(self, val):
        self.fit_dm2 = val
        
    def _set_fit_sin2(self, val):
        self.fit_sin2 = val

    def _set_data_filename(self, name):
        self.data_filename = name
        
    def _set_pdf_filename(self, name):
        self.pdf_filename = name

    def _set_bin_width(self, val):
        self.bin_width = val

    ######################################################################## 
    ## Gettters:
    def _get_data_dm2(self):
        return self.data_dm2
    
    def _get_data_sin2(self):
        return self.data_sin2
    
    def _get_fit_dm2(self):
        return self.fit_dm2
    
    def _get_fit_sin2(self):
        return self.fit_sin2

    ######################################################################## 
    # file loaders  
    def load_data_file(self, evno=[]):
        '''
        Load toy MC data from data_filename.
        '''
        if not os.path.exists(self.data_filename):
            print(f'{self.data_filename} does not exist :(')
            sys.exit(-1)

        with h5.File(self.data_filename, 'r') as f:

            if len(evno) == []:
                dataset_names = f.keys()
                for name in dataset_names:
                    self.dataset.append( f[name][:]) 
            
            else:
                for i in evno:
                    name = f'event{i}'
                    if name not in f.keys():
                        continue
                    self.dataset.append( f[name][:])

        self.load_data_flag = True 

    def load_other_experiments_sensitivity(self):
        for file in self.other_exp_files:
            if not os.path.exists(file):
                print(f"{file} does not exist !")
                continue
            arr = np.loadtxt(file)
            self.other_exp_sens.append( arr )
        
    ######################################################################## 
    def expected_signal_count_onebin(self, bl, R):
        h_fit = self.MC_gen.generate_asimov_dataset(data_dm2=self.fit_dm2, data_sin2=self.fit_sin2, coarse_step_bl=self.bin_width)
        f = lambda x: np.interp(x, h_fit.centers[0], h_fit.values)
        return f(bl) * R
    
    ######################################################################## 
    ### calculate dchi2 for Asimov dataset: total 3 ways for now.
    
    def calculate_rateonly_dchi2_asimov(self):
        ## There is actually no only fitting, just calculating rate-only delta_chi2 for (sin2, dm2) pairs.
        # Fitting histogram
        h_fit = self.MC_gen.generate_asimov_dataset(data_dm2=self.fit_dm2, data_sin2=self.fit_sin2, coarse_step_bl=self.bin_width)
        fit_nevt = np.sum( h_fit.values )
        
        # Data asimov histogram
        h_data = self.MC_gen.generate_asimov_dataset(data_dm2=self.data_dm2, data_sin2=self.data_sin2, coarse_step_bl=self.bin_width)
        data_nevt = np.sum( h_data.values )

        # For now, only consider Poisson statistics fluctuation.
        sigma_data = np.sqrt( data_nevt )
        
        chi2 = (fit_nevt - data_nevt)**2 / sigma_data**2

        return chi2
        

    def calculate_shape_only_dchi2_asimov_fixedRate(self):
        ## There is actually no only fitting, just calculating rate-only delta_chi2 for (sin2, dm2) pairs.
        # Fitting histogram
        h_fit = self.MC_gen.generate_asimov_dataset(data_dm2=self.fit_dm2, data_sin2=self.fit_sin2 , coarse_step_bl=self.bin_width)
        fit_spec = h_fit.values
        
        # Data asimov histogram
        h_data = self.MC_gen.generate_asimov_dataset(data_dm2=self.data_dm2, data_sin2=self.data_sin2, coarse_step_bl=self.bin_width )
        data_spec = h_data.values

        # For now, only consider Poisson statistics fluctuation.
        sigma_data = np.sqrt( data_spec )
        
        mask_id = np.where(sigma_data != 0)
        chi2 = np.sum( (fit_spec[mask_id] - data_spec[mask_id])**2 / sigma_data[mask_id]**2 ) 

        return chi2

    def calculate_shape_only_dchi2_asimov_scanRate(self, coarse_scan_step=0.1, coarse_Nscan=21, fine_scan_step=0.001, fine_Nscan=21, draw=False):

        def calculate_scan_range(center, scan_step, scan_number):
            x = np.zeros(scan_number)
            N_oneSide = int((scan_number - 1) / 2)
            low = center - N_oneSide * scan_step
            for i in range(scan_number):
                Ri = low + i * scan_step
                x[i] = Ri
            return x
            
        def scanning(x, data, fit, err):
            y = np.zeros(len(x))
            for i in range(len(x)):
                R = x[i]
                dchi2 = np.sum( (data - R*fit)**2/err**2 )
                y[i] = dchi2 
            return y
                
        h_fit = self.MC_gen.generate_asimov_dataset(data_dm2=self.fit_dm2, data_sin2=self.fit_sin2, coarse_step_bl=self.bin_width)
        fit_spec = h_fit.values
        
        # Data asimov histogram
        h_data = self.MC_gen.generate_asimov_dataset(data_dm2=self.data_dm2, data_sin2=self.data_sin2, coarse_step_bl=self.bin_width)
        data_spec = h_data.values

        sigma_data = np.sqrt( data_spec )
        
        mask_id = np.where( sigma_data != 0 )

        
        ## Coarse scanning
        find_coarse_best = False
        find_coarse_time = 0
        start = 1.0
        while (not find_coarse_best) and (find_coarse_time < 10):
            coarse_R = calculate_scan_range(start, coarse_scan_step, coarse_Nscan)
            coarse_dchi2 = scanning(coarse_R, data_spec[mask_id], fit_spec[mask_id], sigma_data[mask_id])
            min_idx = np.argmin(coarse_dchi2) 
            if min_idx == 0:
                #print('xxxxx Scanning not good -> minimum at the left-edge of the coarse scanning range !')
                start = coarse_R[0]
            elif min_idx == len(coarse_R) - 1:
                #print('xxxxx Scanning not good -> minimum at the right-edge of the coarse scanning range !')
                start = coarse_R[-1]
            else:
                find_coarse_best = True
            find_coarse_time += 1
                
        if (not find_coarse_best):
            print(f"xxx Something wrong with coarse scanning -> could not find minimum after {find_coarse_time} times scanning.")

        # Fine scanning
        find_fine_best = False
        find_fine_time = 0
        start = coarse_R[min_idx]
        while (not find_fine_best) and (find_fine_time < 10):
            fine_R = calculate_scan_range(start, fine_scan_step, fine_Nscan)
            fine_dchi2 = scanning(fine_R, data_spec[mask_id], fit_spec[mask_id], sigma_data[mask_id])
            min_idx = np.argmin(fine_dchi2) 
            if min_idx == 0:
                #print('xxxxx Scanning not good -> minimum at the left-edge of the fine scanning range !')
                start = fine_R[0]
            elif min_idx == len(fine_R) - 1:
                #print('xxxxx Scanning not good -> minimum at the right-edge of the fine scanning range !')
                start = fine_R[-1]
            else:
                find_fine_best = True
            find_fine_time += 1

        if (not find_fine_best):
            print(f"xxx Something wrong with fine scanning -> could not find minimum after {find_fine_time} times scanning.")
        
        fine_min_R = fine_R[min_idx]
            
        if draw:
            bins = h_fit.bins[0]
            values = h_fit.values * fine_min_R
            h_fit_scaled = hl.Hist(bins, values)
            fig = self.draw_fits([h_fit_scaled, h_fit, h_data], ['scaled fit', 'fit', 'data'])
            return coarse_R, coarse_dchi2, fine_R, fine_dchi2, fig
        
        return coarse_R, coarse_dchi2, fine_R, fine_dchi2



    def calculate_shape_only_dchi2_asimov_fitRate(self, draw=False):
        fit_valid = False
        N_fit, N_fit_max = 0, 5
        while (not fit_valid) and (N_fit < N_fit_max):
            h_fit  = self.MC_gen.generate_asimov_dataset(data_dm2=self.fit_dm2, data_sin2=self.fit_sin2,  coarse_step_bl=self.bin_width)
            h_data = self.MC_gen.generate_asimov_dataset(data_dm2=self.data_dm2, data_sin2=self.data_sin2, coarse_step_bl=self.bin_width)

            xe = h_data.bins[0]
            xc = (xe[1:]+xe[:-1])/2.
            vals = h_data.values
            errs = np.sqrt( vals )
            mask_id = np.where(vals != 0 )
            xc = xc[mask_id]
            vals = vals[mask_id]
            errs = errs[mask_id]
            c = cost.LeastSquares(xc, vals, errs, self.expected_signal_count_onebin )
            m = Minuit(c, R=1.0,)
            m.limits['R'] = (0, 2)

            m.migrad()
            fit_valid = m.valid
            N_fit += 1
        
        if not fit_valid:
            print(f'Fitting failed after {N_fit} time trys.')

        if draw:
            fit_R = m.values[0]
            bins = h_fit.bins[0]
            values = h_fit.values * fit_R
            h_fit_scaled = hl.Hist(bins, values)
            fig = self.draw_fits([h_fit_scaled, h_fit, h_data], ['scaled fit', 'fit', 'data'])
            return m.values[0], m.errors[0], m, fig
        
        return m.values[0], m.errors[0], m
        

        
    ######################################################################## 
    ### Fit toy MC fluctuated dataset
    
    def calculate_rateonly_dchi2_MCdata(self):
        ## There is actually no only fitting, just calculating rate-only delta_chi2 for (sin2, dm2) pairs.
        # Fitting histogram
        h_fit = self.MC_gen.generate_asimov_dataset(data_dm2=self.fit_dm2, data_sin2=self.fit_sin2, coarse_step_bl=self.bin_width)
        fit_nevt = np.sum( h_fit.values )
        
        if not self.load_data_flag:
            self.load_data_file()
        if len( self.dataset ) == 0:
            print("There is no dataset loaded yet.")
            sys.exit(-1)

        for dat in self.dataset:
            data_nevt = len(dat)
            sigma_nevt = np.sqrt( data_nevt )
            
        

    ######################################################################## 
    ## Draw fits
    def draw_fits(self, histo, labels):
        fig, ax = plt.subplots(figsize=(8, 6))
        for h, lb in zip(histo, labels):
            hl.plot1d(ax, h, label=lb)
        ax.set_xlabel('Baseline [m]', fontsize=14)
        ax.set_ylabel('Count per bin', fontsize=14)
        ax.legend(fontsize=13)
        ax.tick_params(labelsize=14)
        plt.tight_layout()
        return fig