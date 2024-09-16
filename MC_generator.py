import numpy as np
from scipy.integrate import dblquad, quad, tplquad, nquad
from scipy.interpolate import interp2d, LinearNDInterpolator
from histlite import hist, Hist
import h5py as h5
import os
from scipy.stats import norm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from IPython import get_ipython
def isnotebook():
        try:
                shell = get_ipython().__class__.__name__
                if shell == 'ZMQInteractiveShell':
                        return True   # Jupyter notebook or qtconsole
                elif shell == 'TerminalInteractiveShell':
                        return False  # Terminal running IPython
                else:
                        return False  # Other type (?)
        except NameError:
                return False      # Probably standard Python interpreter
if(isnotebook()):
        from tqdm.notebook import tqdm
else:
        from tqdm import tqdm

from signal_spectrum import signal_calculation
from detector import detector
from neutrino_source import neutrino_source
from oscillation import *
from xsection import xsection

from unit_conversion import *

class MC_generator:
    def __init__(self, source, det, dm2=1.0, sin2theta_square=0.1, int_type='nue', seed=42):
        self.source = source
        self.det = det
        self.xsec = xsection()
        self.interaction =  int_type # 'nue' or CC, only two options now
        
        #### calculate baseline ranges:
        self.baseline_min = np.abs(self.det.position[2] - self.source.position[2]) - self.det.height/2.
        self.baseline_max = np.sqrt( (np.abs(self.det.position[2]-self.source.position[2]) + self.det.height/2.)**2 + self.det.radius**2 )
        
        # neutrino mixing parameters:
        self.dm2 = dm2
        self.sin2theta_square = sin2theta_square

        # expected un-oscillated neutrino number:
        #self.source.scale_counts(self.det.height. self.det.volume, self.det.baseline, self.det.run_time )
        self.n_events_noosc = self.scale_counts()
        
        self.h_maps = None
        
        self.scale_noosc = 11.978 # hard-coded now,

        np.random.seed(seed)

        # If I do the full-radius-range integral as in the following codes, I get a normalization factor as,
        self.integral_norm_factor = 1.0
        
        
        ### Load 
        self.osc_event_rate_file = '/p/lustre1/yu47/Sterile_Neutrino/sensitivity/data/event_rate_bl10cm_Ev750keV.h5'
        self.osc_event_rate_func = None


    def scale_counts(self):
        # Scale the number of neutrinos generated by the source from the LZ paper (JHEP11(2014) 042)
        # Approximately take the average arrived neutrinos scaled by the square of distances from the source to the detector center (L+H/2.)**2
        # The LZ detector is hard-coded here.
        #H0 = 1.38 # unit: m, LZ detector height
        #R0 = H0 / 2. # unit: m, LZ detector radius
        #V0 = np.pi * R0**2 * H0 # unit: m^3, LZ detector volume
        #L0 = 1 # unit: m, below the fiducial volume, along the central axis of the cylinder. Distance from the source to the bottom plane of the cylinder.
        
        N0 = 12518 # 100 days run
        T0 = 100 # unit: days, 100 days run
        activity0 = 5e6 # unit: Ci, 1 Ci = 3.7e10 dacays per second
        #geometry_factor_LZ =  0.37128 # hard-coded, I did convolve_oscillation_event_rate for LZ detector without oscillation

        spheric_geometry_factor = 494.620
        spheric_geometry_factor_LZ = 123.382

        geometry_factor_ratio = 4.16
        
        #H = self.det.height
        #V = self.det.volume
        #L = self.det.baseline - self.det.height/2.
        #L = np.sqrt((self.source.position[0]-self.det.position[0])**2 + (self.source.position[1]-self.det.position[1])**2 + (self.source.position[2]-self.det.position[2])**2 ) - self.det.height/2.
        T = self.det.run_time
        
        Enu = self.source.energies[0]
        #geometry_factor = self.convolve_oscillation_event_rate(Enu)[0]
        br  = self.source.ratios[0]
        xsec_nue = self.xsec.total_xsec_nuescatter(Enu)
        xsec_cc = self.xsec.total_xsec_CC(Enu)

        totxsec = 0.
        if self.interaction == 'nue':
            totxsec = xsec_nue
        elif self.interaction == 'CC':
            totxsec = xsec_cc
        else:
            print("There is no such interaction in this package -> must in {'CC', 'nue'}.")

        if False:
            print(f'Pre-scaling event count: {N0:.3f}')
            print(f"Scaling exposure time: {T:.3f}, {T0:.3f}")
            print(f'Geometry factor: {spheric_geometry_factor:.3f}, {spheric_geometry_factor_LZ:.3f}')
            print(f"Activity: {self.source.activity}, {activity0}")
            print(f"Cross section: {totxsec:.3e}, {xsec_nue:.3e}")

        self.n_events_noosc = N0 * (T/T0) * geometry_factor_ratio * (self.source.activity/activity0) * br * (totxsec/xsec_nue)
        #self.n_events_noosc = N0 * (T/T0) * (spheric_geometry_factor/spheric_geometry_factor_LZ) * (self.source.activity/activity0) * br * (totxsec/xsec_nue)
        #self.n_events_noosc = N0 * (T/T0) * (geometry_factor/geometry_factor_LZ) * (self.source.activity/activity0) * br * (totxsec/xsec_nue)
        #self.n_events_noosc = (N0 / V0 / T0 * (L0+ H0/2.)**2 ) * V * T / (L+H/2.)**2 * (self.source.activity / activity0) * br * totxsec / xsec_nue
        return self.n_events_noosc
    
    def draw_experiment_layout(self, events=[], elev=0, azim=90, roll=0):
        # Create meshgrid for the source
        theta = theta = np.linspace(0, 2 * np.pi, 50)
        z0 = np.linspace(self.source.position[2] - self.source.height/2., self.source.position[2] + self.source.height/2., 50)
        z1 = np.linspace(self.det.position[2] - self.det.height/2., self.det.position[2] + self.det.height/2., 50)
        theta_grid0, z_grid0 = np.meshgrid(theta, z0)    
        theta_grid1, z_grid1 = np.meshgrid(theta, z1)    
        x_grid0 = self.source.diameter/2. * np.cos(theta_grid0)
        y_grid0 = self.source.diameter/2. * np.sin(theta_grid0)
        x_grid1 = self.det.radius * np.cos(theta_grid1)
        y_grid1 = self.det.radius * np.sin(theta_grid1)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(x_grid0, y_grid0, z_grid0, color='b', alpha=0.2, rstride=5, cstride=5) 
        ax.plot_surface(x_grid1, y_grid1, z_grid1, color='r', alpha=0.2, rstride=5, cstride=5) 

        for evt in events:
            xs, ys, zs = evt[0::3], evt[1::3], evt[2::3]
            ax.plot3D(xs, ys, zs, color='gray')
        
        ax.view_init(elev=elev, azim=azim, roll=roll)

        plt.show()
        return fig
    
    ## getter methods:

    def _get_source(self):
        return self.gen
    
    def _get_detector(self):
        return self.det
    
    def _get_nevents(self):
        return self.n_events_noosc
    
    def _get_baseline_range(self):
        return self.baseline_min, self.baseline_max

    def _get_oscillation_length(self):
        return oscillation_length(self.dm2, self.source.energies[0])

    ####################################
    ### setter methods:
    def _set_dm2(self, dm2):
        self.dm2 = dm2
        
    def _set_sin2theta_square(self, sin2):
        self.sin2theta_square = sin2

    def _set_event_rate_filename(self, filename):
        self.osc_event_rate_file = filename

    ####################################
    ### Generating events uniformly in the detector, which could be incorrect because I did not consider the flux scaling in unit volumes with different distances.

    def generate_in_cubic(self, x0, y0, z0, H0, L0, N0):
        x = np.random.uniform(-L0/2+x0, L0/2+x0, size=N0)
        y = np.random.uniform(-L0/2+y0, L0/2+y0, size=N0)
        z = np.random.uniform(-H0/2+z0, H0/2+z0, size=N0)
        return x, y, z
    
    def Is_in_detector(self, x, y, z):
        if ( np.sqrt((x-self.det.position[0])**2 + (y-self.det.position[1])**2) < self.det.radius) and (-self.det.height/2. <= z-self.det.position[2] < self.det.height/2.) :
            return True
        else:
            return False

    def Is_in_source(self, x, y, z):
        if ( np.sqrt((x-self.source.position[0])**2 + (y-self.source.position[1])**2) < self.source.diameter/2.) and (-self.source.height/2. <= z-self.source.position[2] < self.source.height/2.) :
            return True
        else:
            return False

    def smear_position(self, pos, res):
        x_smeared = np.random.normal(loc=pos[:,0], scale=res, size=len(pos))
        y_smeared = np.random.normal(loc=pos[:,1], scale=res, size=len(pos))
        z_smeared = np.random.normal(loc=pos[:,2], scale=res, size=len(pos))
        pos_smeared = np.array([x_smeared, y_smeared, z_smeared]).T
        bl = np.sqrt( (x_smeared - self.source.position[0])**2 + (y_smeared - self.source.position[1])**2 + (z_smeared - self.source.position[2])**2 )
        return pos_smeared, bl


    # Do numerical integral to as: int_V (Phi) dV \propto int_V (P(V)/L^2) dV (where P is the survival probability, L is the distance) 
    # Considering flux scaling at each unit volume, and also the survival probability 
    def flux_scaling(self, bl):
        # scaled to 1m baseline
        return 1. / bl**2


    def expected_flux_UnitSourceCylindricalVolume(self, r, z, theta, E, dm2, sin2theta_square, xd, yd, zd):
        #(xd, yd, zd) is the point coordinate in the detector we are calculating now.
        xs, ys, zs = r*np.cos(theta), r*np.sin(theta), z # (xs, ys, zs) is the unit volume coordinate in the cylindrical source.
        bl = np.sqrt((xs-xd)**2 + (ys-yd)**2 + (zs-zd)**2 )
        P = electron_neutrino_survival_probability(dm2, sin2theta_square, E, bl)
        flux = self.flux_scaling(bl)
        return P * flux # There is no volume here.
    
    
    def expected_rate_UnitSourceCylinderVolume_UnitDetectorSphericalVolume(self, E, dm2, sin2, rs, thetas, zs, rd, thetad, phid):
        xs, ys = rs * np.cos(thetas), rs * np.sin(thetas)
        xd, yd, zd = rd * np.sin(thetad) * np.cos(phid), rd * np.sin(thetad) * np.sin(phid), rd * np.cos(thetad)
        l = np.sqrt((xs-xd)**2 +(ys-yd)**2 + (zs-zd)**2)
        P = electron_neutrino_survival_probability(dm2, sin2, E, l)
        F = self.flux_scaling(l)
        vol_source = rs
        vol_det = rd**2 * np.sin(thetad)
        return P * F * vol_source * vol_det
        
    
        
    def expected_rate_UnitDetectorSphericVolume_withSourceGeometry_noNumericalIntegral(self, rs, phi, theta, E, dm2, sin2):
        xd, yd, zd = rs*np.sin(theta)*np.cos(phi), rs*np.sin(theta)*np.sin(phi), rs*np.cos(theta) + self.source.position[2]
        num_samples = int(1e6)
        z_samples = np.random.uniform(self.source.position[2]-self.source.height/2., self.source.position[2]+self.source.height/2., num_samples)
        theta_samples = np.random.uniform(0, 2 * np.pi, num_samples)
        r_samples = (self.source.diameter/2.) * np.sqrt(np.random.uniform(0, 1, num_samples))
        value_samples = self.expected_flux_UnitSourceCylindricalVolume(r_samples, z_samples, theta_samples, E, dm2, sin2, xd, yd, zd) * r_samples
        source_cylinder_volume = np.pi * (self.source.diameter/2.)**2 * self.source.height
        integral_estimate = np.mean(value_samples) * source_cylinder_volume
        return integral_estimate

  
    def expected_rate_UnitDetectorSphericVolume_withSourceGeometry(self, rs, phi, theta, E, dm2, sin2theta_square):
        ## In shperical coordinate, the origin is set as the center of the source so that the radius = baseline
        xd, yd, zd = rs*np.sin(theta)*np.cos(phi), rs*np.sin(theta)*np.sin(phi), rs*np.cos(theta) + self.source.position[2]
        f = lambda r, z, theta: self.expected_flux_UnitSourceCylindricalVolume(r, z, theta, E, dm2, sin2theta_square, xd, yd, zd) * r
        res, err = tplquad(f, 0, 2*np.pi, self.source.position[2]-self.source.height/2., self.source.position[2]+self.source.height/2., 0, self.source.diameter/2.)
        return res, err


    def expected_rate_UnitDetectorSphericVolume_noSourceGeometry(self, rs, E, dm2, sin2theta_square):
        flux = self.flux_scaling(rs)
        P = electron_neutrino_survival_probability(dm2, sin2theta_square, E, rs)
        return P * flux


    def expected_rate_DetectorSphericalAngularIntegral_withSourceGeometry(self, rs, E, dm2, sin2theta_square, mc=False):
        theta_low, theta_high = self.spheric_theta_limits(rs)
        if mc:
            f = lambda phi, theta: self.expected_rate_UnitDetectorSphericVolume_withSourceGeometry_noNumericalIntegral(rs, phi, theta, E, dm2, sin2theta_square) * np.sin(theta) * rs**2
        else:
            f = lambda phi, theta : self.expected_rate_UnitDetectorSphericVolume_withSourceGeometry(rs, phi, theta, E, dm2, sin2theta_square)[0] * np.sin(theta) * rs**2
        res, err = dblquad(f, theta_low, theta_high, 0, 2*np.pi)
        return res, err

    def expected_rate_DetectorSphericalVolumeIntegral_withSourceGeometry(self, E, dm2, sin2theta_square, mc=False):
        if mc:
            f = lambda phi, theta, r: self.expected_rate_UnitDetectorSphericVolume_withSourceGeometry_noNumericalIntegral(r, phi, theta, E, dm2, sin2theta_square) * np.sin(theta) * r**2
        else:
            f = lambda phi, theta, r : self.expected_rate_UnitDetectorSphericVolume_withSourceGeometry(r, phi, theta, E, dm2, sin2theta_square)[0] * np.sin(theta) * r**2
        res, err = tplquad(f, self.baseline_min, self.baseline_max, lambda r: self.spheric_theta_limits(r)[0], lambda r: self.spheric_theta_limits(r)[1], 0, 2*np.pi)

        return res, err
    
    def expected_rate_DetectorSphericalAngularIntegral_noSourceGeometry(self, rs, E, dm2, sin2):
        f = lambda theta : self.expected_rate_UnitDetectorSphericVolume_noSourceGeometry(rs, E, dm2, sin2) * rs**2 * np.sin(theta)
        theta_low, theta_high = self.spheric_theta_limits(rs)
        res, err = quad(f, theta_low, theta_high)
        return res, err

        
    def expected_rate_DetectorSphericalVolumeIntegral_noSourceGeometry(self, E, dm2, sin2):
        f = lambda theta, r : self.expected_rate_UnitDetectorSphericVolume_noSourceGeometry(r, E, dm2, sin2) * r**2 * np.sin(theta)
        res, err = dblquad(f, self.baseline_min, self.baseline_max, lambda r: self.spheric_theta_limits(r)[0], lambda r: self.spheric_theta_limits(r)[1])
        return res, err

    
    ### The high-dimensional integral is very computation-heavy so that we should -> we need to do a tplquad for the source volume and then a dbl/tplquad for the detector.
    def expected_rate_MonteCarloIntegral(self, E, dm2, sin2):
        dist = np.abs( self.det.position[2] - self.source.position[2] )
        #f = lambda zs, thetas, rs, phi, theta, r: self.expected_flux_UnitSourceCylindricalVolume(rs, zs, thetas, E, dm2, sin2, r*np.sin(theta)*np.cos(phi), r*np.sin(theta)*np.sin(phi), r*np.cos(theta)-(dist)) * electron_neutrino_survival_probability(dm2, sin2, E, ) * rs * r**2 * np.sin(theta)
        pass
        
    def expected_rate_DetectorSphericalAngularIntegral_withSourceGeometry_1(self, r, E, dm2, sin2):
        f = lambda thetad, zs, thetas, rs: self.expected_rate_UnitSourceCylinderVolume_UnitDetectorSphericalVolume(E, dm2, sin2, r, thetas, zs, rs, thetad, 0.)
        zmin, zmax = self.source.position[2] - self.source.height/2., self.source.position[2] + self.source.height/2.
        res, err = nquad(f, [[self.spheric_theta_limits(r)[0],  self.spheric_theta_limits(r)[1]], [zmin, zmax], [0, 2*np.pi], [0, self.source.diameter/2.] ])
        res, err = nquad()
        return res, err
        
        
        
        
        
  
    ##################################################
    ##### Previous codes
    ### Integral inside the detector
    def expected_rate_UnitCylindricalVolume(self, r, z, E, dm2, sin2theta_square):
        bl = np.sqrt(r**2 + (z - self.source.position[2])**2) 
        P = electron_neutrino_survival_probability(dm2, sin2theta_square, E, bl)
        flux = self.flux_scaling(bl)
        return  P * flux

    def expected_rate_cylindrical_integral(self, E, data_dm2, data_sin2):
        f = lambda r, z, E: self.expected_rate_UnitCylindricalVolume(r, z, E, data_dm2, data_sin2) * r * 2*np.pi
        res, err = dblquad(f, -self.det.height/2., self.det.height/2., lambda z: 0., lambda z: self.det.radius, args=(E, ), )
        return res, err

    def expected_rate_UnitSphericVolume(self, rs, E, dm2, sin2theta_square):
        P = electron_neutrino_survival_probability(dm2, sin2theta_square, E, rs)
        flux = self.flux_scaling(rs)
        return P * flux

    def expected_rate_spheric_integral(self, E, data_dm2, data_sin2):
        f = lambda theta, r, E: self.expected_rate_UnitSphericVolume(r, E, data_dm2, data_sin2) * r**2 * np.sin(theta) * 2*np.pi
        res, err = dblquad(f, self.baseline_min, self.baseline_max, lambda r: self.spheric_theta_limits(r)[0], lambda r: self.spheric_theta_limits(r)[1], args=(E, ))
        return res, err
        
        
    def convolve_resolution(self, h0, rebin_width=0.03):
        bin_edges = h0.bins[0]
        bin_conts = h0.values
        bin_cents = h0.centers[0]
        smeared_bin_cents = np.zeros( len(bin_cents) )
        lowbin, higbin, binwidth = bin_edges[0], bin_edges[-1], bin_edges[1]-bin_edges[0]
        nsigma = 5 # only calculate 5-sigma region smearing
        for ibin, ibincenter in enumerate(bin_cents):
            lowi  = max([lowbin, ibincenter-nsigma*self.det.spatial_resolution])
            lowid = int((lowi-lowbin)/binwidth)
            highi = min([higbin, ibincenter+nsigma*self.det.spatial_resolution])
            highid = int((highi-lowbin)/binwidth)
            for ibin_smear in np.arange(lowid, highid):
                smeared_bin_cents[ibin_smear] += bin_conts[ibin] * norm.pdf(bin_cents[ibin_smear], loc=ibincenter, scale=self.det.spatial_resolution) * binwidth

        smeared_h0 = Hist(bin_edges, smeared_bin_cents)

        baseline_edges_new = np.arange(bin_edges[0], bin_edges[-1], rebin_width)
        baseline_edges_new = np.append(baseline_edges_new, bin_edges[-1])
        smeared_h = smeared_h0.rebin(0, baseline_edges_new)
        return smeared_h
        


    def generate_dataset(self, data_dm2=0.0, data_sin2=0.0, poisson=False, smear=False):
        # The current acceptance-rejection sampling method seems not efficient :( Around 5 sec per event
        n_events_generated = 0
        gen_pos = []
        gen_bl  = []
        h = self.generate_asimov_dataset(data_dm2, data_sin2)
        total_expected_event = np.sum(h.values)
        if poisson:
            total_expected_event = np.random.poisson( total_expected_event )
        while n_events_generated < total_expected_event:
            x0, y0, z0 = self.generate_in_cubic(*self.det.position, self.det.height, 2*self.det.radius, 1)
            x0, y0, z0 = x0[0], y0[0], z0[0]
            if not self.Is_in_detector(x0, y0, z0):
                continue
            bl = np.sqrt( x0**2 + y0**2 + (z0-self.source.position[2])**2 )
            flux_max = self.flux_scaling(self.baseline_min)
            flux = self.flux_scaling(bl)
            flux_prob = flux / flux_max
            surprob = electron_neutrino_survival_probability( data_dm2, data_sin2, self.source.energies[0], bl)
            p0 = np.random.uniform()
            if p0 < flux_prob * surprob:
                n_events_generated += 1
                gen_pos.append( [x0, y0, z0] )
                gen_bl.append( bl )
            else:
                continue

        gen_pos = np.array(gen_pos)
        gen_bl = np.array(gen_bl)
        if smear:
            gen_pos, gen_bl = self.smear_position(gen_pos, self.det.spatial_resolution)

        return gen_pos, gen_bl


    def generate_asimov_dataset(self, data_dm2=0.0, data_sin2=0.0, fine_step_bl=0.001, coarse_step_bl=0.03, source_geom=False, smear=False):
        #total_expected_event = self.interp_oscillated_event_rate(data_dm2, data_sin2)
        #print(f"Total {total_expected_event} events for oscillation parameters ({data_dm2:.4f}, {data_sin2:.4f}).")
        baseline_edges = np.arange(self.baseline_min, self.baseline_max+fine_step_bl, fine_step_bl)
        baseline_cents = (baseline_edges[1:] + baseline_edges[:-1]) / 2.
        counts = []
        for bl in tqdm( baseline_cents ):
            #tmp_shell_counts = self.theta_integral(bl, data_dm2, data_sin2)[0]
            if source_geom:
                tmp_shell_counts = self.expected_rate_DetectorSphericalAngularIntegral_withSourceGeometry(bl, self.source.energies[0], data_dm2, data_sin2, mc=False)[0]
            else:
                tmp_shell_counts = self.expected_rate_DetectorSphericalAngularIntegral_noSourceGeometry(bl, self.source.energies[0], data_dm2, data_sin2)[0]
            counts.append(tmp_shell_counts)
        #scale_factor = total_expected_event / np.sum(counts)
        #print(total_expected_event, np.sum(counts), scale_factor)
        counts = np.array(counts)
        #counts = counts * scale_factor
        counts = counts * self.scale_noosc
        tmp_hist = Hist(baseline_edges, counts)

        if smear:
            # will smear the positions only for now to consider the detector spatial resolution
            smeared_h = self.convolve_resolution(tmp_hist, rebin_width=coarse_step_bl)
            return smeared_h

        # Rebin the original histogram into coarse bins
        baseline_edges_new = np.arange(baseline_edges[0], baseline_edges[-1], coarse_step_bl)
        baseline_edges_new = np.append(baseline_edges_new, baseline_edges[-1])
        new_hist = tmp_hist.rebin(0, baseline_edges_new)
        return new_hist
        
            
            
        

            
    def spheric_theta_limits(self, rs):
        z_min, z_max = -self.det.height/2., self.det.height/2.
        rs_min, rs_max = -self.source.position[2]-self.det.height/2., np.sqrt((-self.source.position[2]+self.det.height/2.)**2 + self.det.radius**2)
        dist = z_min - self.source.position[2]
        theta_max = np.arctan(self.det.radius / (dist))

        if False:
            print('******************** General limits on rs and theta **********************')
            print(f'theta_c in [0, {theta_max/np.pi*180:.4f}]')
            print(f'r_c in [{rs_min:.3f}, {rs_max:.3f}] ')
            print('**************************************************************************')
    
    
        if rs > rs_max or rs < rs_min:
            print(f'Radius is out of range [{rs_min:.3f}, {rs_max:.3f}].')
            return 0, 0

        theta_uplimt0 = 100.
        if self.det.radius < rs:
            #print('Cylindrical radisu condition should be satisfied !')
            theta_uplimt0 = np.arcsin(self.det.radius/rs)
    
        theta_uplimt1 = 100.
        up_cond = (-self.det.height/2.-self.source.position[2]) / rs
        if 0<= up_cond <= 1:
            #print('Lower cylindrical height condition should be satisfied !')
            theta_uplimt1 = np.arccos(up_cond)
        theta_lowlimt1 = -100
        low_cond = (self.det.height/2.-self.source.position[2]) / rs
        if 0<= low_cond <= 1:
            #print('Upper cylindrical height condition should be satisfied !')
            theta_lowlimt1 = np.arccos(low_cond)

        if False:
            print('******************** Strict limits on rs and theta **********************')
            print(f'Upper limits of theta: {theta_uplimt0:.4f}, {theta_uplimt1:.4f}.')
            print(f"Lower limits of theta: {theta_lowlimt1:.4f}.") 
    
    
        uplimit = min([theta_max, theta_uplimt0, theta_uplimt1])
        lowlimit = max([0, theta_lowlimt1])
        return lowlimit, uplimit 

        
        
    def theta_integral(self, rs, dm2, sin2theta_square):
        # rs is actually the baseline
        flux = self.flux_scaling(rs)
        surprob = electron_neutrino_survival_probability(dm2, sin2theta_square, self.source.energies[0], rs)
        low, up = self.spheric_theta_limits(rs)
        f = lambda theta, r: flux * surprob * r**2 * np.sin(theta)
        res, err = quad(f, low, up, args=(rs,))
        return res, err

    def radius_integral_fullrange(self, dm2, sin2theta_square):
        f = lambda r: self.theta_integral(r, dm2=dm2, sin2theta_square=sin2theta_square)[0]
        res, err = quad(f, self.baseline_min, self.baseline_max)
        self.integral_norm_factor = res
        return res
        
    def radius_integral(self, rmin, rmax, dm2, sin2theta_square):
        '''
        In my codes, I did theta integral first, then I do radius integral.
        '''
        if self.integral_norm_factor == 0:
            self.radius_integral_fullrange(dm2, sin2theta_square)
        f = lambda r: self.theta_integral(r, dm2=dm2, sin2theta_square=sin2theta_square)[0]
        res, err = quad(f, rmin, rmax, )
        return res/self.integral_norm_factor, err/self.integral_norm_factor
        
    
    def expected_rate_and_integral(self, r, dm2, sin2theta_square):
        if self.integral_norm_factor == 0:
            self.radius_integral_fullrange(dm2, sin2theta_square)
        density     = self.theta_integral(r, dm2, sin2theta_square)[0] / self.integral_norm_factor
        integral    = self.radius_integral(self.baseline_min, r, dm2, sin2theta_square)[0]
        return density, integral

        
    def monte_carlo_sampling(self, n_event=1e5, return_pos=False):
        # First, randomly sample points within the source and detector volume separately.
        x_source, y_source, z_source = [], [], []
        delta_n = int(n_event - len(x_source))
        while delta_n > 0:
            x_source0, y_source0, z_source0 = self.generate_in_cubic(*self.source.position, self.source.height, self.source.diameter, int(delta_n*1.5))
            for x, y, z in zip(x_source0, y_source0, z_source0):
                if self.Is_in_source(x, y, z):
                    x_source.append(x)
                    y_source.append(y)
                    z_source.append(z)
            delta_n = int(n_event - len(x_source))
            
        x_source = np.array(x_source)[0:n_event]        
        y_source = np.array(y_source)[0:n_event]        
        z_source = np.array(z_source)[0:n_event]        
        
        
        x_detector, y_detector, z_detector = [], [], []
        delta_n = int(n_event - len(x_detector))
        while delta_n > 0:
            x_detector0, y_detector0, z_detector0 = self.generate_in_cubic(*self.det.position, self.det.height, self.det.radius*2, int(n_event*1.5))
            for x, y, z in zip(x_detector0, y_detector0, z_detector0):
                if self.Is_in_detector(x, y, z):
                    x_detector.append(x)
                    y_detector.append(y)
                    z_detector.append(z)
            delta_n = int(n_event - len(x_detector))
        x_detector = np.array(x_detector)[0:n_event]        
        y_detector = np.array(y_detector)[0:n_event]        
        z_detector = np.array(z_detector)[0:n_event]        
        
        # match pairs: calculate baseline and weight.
        pairs = []
        x3, y3, z3 = self.source.position[0], self.source.position[1], self.source.position[2]
        for i in range(n_event):
            x1, y1, z1 = x_source[i], y_source[i], z_source[i]
            x2, y2, z2 = x_detector[i], y_detector[i], z_detector[i]
            bl0 = np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
            bl = np.sqrt((x2-x3)**2 + (y2-y3)**2 + (z2-z3)**2)
            weight = self.flux_scaling(bl0) * electron_neutrino_survival_probability(self.dm2, self.sin2theta_square, self.source.energies[0], bl0)
            # I have made a mistake here that I should not use the real baseline here, as we do not know the precise position in the source.
            pairs.append([bl0, bl, weight])

        pairs = np.array(pairs) 
        #return pairs 
        if return_pos:
            return pairs, x_source, y_source, z_source, x_detector, y_detector, z_detector
        else:
            return pairs
        
    
    



# LZ exp
source_LZ = neutrino_source('Cr51', 1e5, [0.75], [1.0])
det_LZ = detector('LZ')
det_LZ.update_geometry(1.38, 1.38)
det_LZ.position = (0, 0, 0)
dist = 1.0
source_LZ.position = ( 0, 0, -det_LZ.height/2.-dist)
det_LZ.run_time = 100

dm2, sin2theta_square = 0.0, 0.0
gen_LZ = MC_generator(source_LZ, det_LZ,  dm2=dm2, sin2theta_square=sin2theta_square, int_type='nue');  