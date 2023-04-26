import numpy as np
import pandas as pd
from numpy import genfromtxt
import matplotlib.pyplot as plt

out_data_file = open("plots/veto_toy_mc_2/data.txt", "w")

#START OF INPUTS *******************************************************************

run_time = 365*10 # [days] data taking exposure

n_sites = 4
site_arr = ["SNOLAB", "LNGS", "WIPP", "SURF"]

#trials = [400, 100, 20, 100]  # No. of simulated exposures to run per site. Must be >= 1
trials = [1, 1, 1, 1]  # No. of simulated exposures to run per site

SNOLAB_flux = 54 # [muons/day] mean muon rate through the OD
LNGS_flux = SNOLAB_flux * 3.41E-4/3.31E-6 # [muons/day] mean muon rate through the OD. Scaled by ratio of Borexino/SNOLAB flux measurements
WIPP_flux = SNOLAB_flux * 4.07E-3/3.31E-6 # [muons/day] mean muon rate through the OD. Scaled by ratio of EXO-200/SNOLAB flux measurements
SURF_flux = SNOLAB_flux * 5.31E-5/3.31E-6 # [muons/day] mean muon rate through the OD. Scaled by ratio of MAJORANA/SNOLAB flux measurements
muon_flux_arr = [SNOLAB_flux, LNGS_flux, WIPP_flux, SURF_flux]

muon_tag_eff = 0.99 # probability of tagging a muon passing through the OD

simulated_muons_snolab = 7.9E6 # 2023_02_10: 2.485E6; 2023_02_17: 7.5E6; 2023_02_24: 7.9E6; 
simulated_muons_lngs = 1.007E8 # 2023_02_10: 52.5E6; 2023_02_17: 6.05E7; 2023_02_24: 1.007E8; 
simulated_muons_wipp = 1.007E8 # 2023_02_10: 52.5E6; 2023_02_17: 6.05E7; 2023_02_24: 1.007E8; - same as LNGS
simulated_muons_surf = 7.9E6 # 2023_02_10: 52.5E6; 2023_02_17: 6.05E7; 2023_02_24: 7.9E6; - previously same as LNGS. now same as SNOLAB
simulated_muons = [simulated_muons_snolab, simulated_muons_lngs, simulated_muons_wipp, simulated_muons_surf]
mult_df_snolab = pd.read_pickle("inputs/capture_multiplicity_snolab_2023_02_24.pkl")# probability distribution of number of captures per muon at SNOLAB. Zeros excluded
mult_df_lngs = pd.read_pickle("inputs/capture_multiplicity_lngs_2023_02_24.pkl")# probability distribution of number of captures per muon at LNGS. Zeros excluded
mult_df_wipp = pd.read_pickle("inputs/capture_multiplicity_lngs_2023_02_24.pkl")# probability distribution of number of captures per muon at WIPP. Zeros excluded. Assume same as LNGS
mult_df_surf = pd.read_pickle("inputs/capture_multiplicity_snolab_2023_02_24.pkl")# probability distribution of number of captures per muon at SURF. Zeros excluded. Assume same as SNOLAB
mult_df_arr = [mult_df_snolab, mult_df_lngs, mult_df_wipp, mult_df_surf]

mult_arr_snolab = mult_df_snolab.to_numpy(dtype="int64")
mult_arr_lngs = mult_df_lngs.to_numpy(dtype="int64")
mult_arr_wipp = mult_df_wipp.to_numpy(dtype="int64")
mult_arr_surf = mult_df_surf.to_numpy(dtype="int64")
mult_arr = [mult_arr_snolab, mult_arr_lngs, mult_arr_wipp,  mult_arr_surf]

muon_nonzero_cap_frac_snolab = len(mult_arr_snolab)/simulated_muons_snolab #fraction of muons with non-zero captures (others were removed from file)
muon_nonzero_cap_frac_lngs = len(mult_arr_lngs)/simulated_muons_lngs #fraction of muons with non-zero captures (others were removed from file)
muon_nonzero_cap_frac_wipp = len(mult_arr_wipp)/simulated_muons_wipp #fraction of muons with non-zero captures (others were removed from file)
muon_nonzero_cap_frac_surf = len(mult_arr_surf)/simulated_muons_surf #fraction of muons with non-zero captures (others were removed from file)
muon_nonzero_cap_frac_arr = [muon_nonzero_cap_frac_snolab, muon_nonzero_cap_frac_lngs, muon_nonzero_cap_frac_wipp, muon_nonzero_cap_frac_surf]
print("Fraction of muons with non-zero captures: ", site_arr, muon_nonzero_cap_frac_arr, file=out_data_file)


#Load capture tagging efficiency curves
bin_centers = np.load("inputs/cap_eff_bin_cent.npy")
cap_tag_file = "inputs/cap_eff_cuml_[Iso][EneVar][Reg][Ene]_2023_02_24.npy" 
cap_tag_eff = np.load(cap_tag_file)

#veto window multiplitcity thresh - open veto window if x neutron captures detected, even if no muon detected
veto_mult_thresh = 3

#Scan different energy cuts for neutron captures
n_cap_tag_emin_scan = [0, 400, 1400, 1800, 2600, 400, 1400] #[2500, 3000, 3500] # [keV] Must be multiples of 40 keV
n_cap_tag_emax_scan = [4600, 4600, 4600, 4600, 4600, 9000, 9000] #[4250, 4250, 4250] # [keV] Must be multiples of 40 keV
n_scans = len(n_cap_tag_emin_scan)

cap_tag_vol_idx = 0 # 0: All LXe, 1: Active LXe
cap_tag_vol = ["All LXe", "Active LXe"]

cap_tag_enevar_idx = 0
cap_tag_enevar = ["MCTruth 1%", "LightOnly"]


#END OF INPUTS *******************************************************************

print(file=out_data_file)
mean_muon_arr = [0] * n_sites #initialize to zero
muon_xe137_cap_prod_arr = [0] * n_sites
for h in range(n_sites):
    mean_muon_arr[h] = muon_flux_arr[h]*run_time
    print("Mean number of muons at " + site_arr[h] + ": ", mean_muon_arr[h], file=out_data_file)

    plt.figure()
    max_mult = mult_df_arr[h].max().max()
    plt.hist(mult_df_arr[h][["Xe137Count"]].to_numpy(dtype="int64"), bins=np.arange(max_mult+1), alpha=0.8, histtype='step', label="Xe137")
    plt.hist(mult_df_arr[h][["Cu64Count"]].to_numpy(dtype="int64"), bins=np.arange(max_mult+1), alpha=0.8, histtype='step',  label="Cu64")
    plt.hist(mult_df_arr[h][["Cu66Count"]].to_numpy(dtype="int64"), bins=np.arange(max_mult+1), alpha=0.8, histtype='step', label="Cu66")
    plt.hist(mult_df_arr[h][["F20Count"]].to_numpy(dtype="int64"), bins=np.arange(max_mult+1), alpha=0.8, histtype='step', label="F20")
    plt.hist(mult_df_arr[h][["H2Count"]].to_numpy(dtype="int64"), bins=np.arange(max_mult+1), alpha=0.8, histtype='step', label="H2")
    plt.xlabel("Captures/Muon given atleast one capture, at " + site_arr[h])
    plt.yscale("log")
    plt.legend()
    plt.savefig("plots/veto_toy_mc_2/captures_per_muon_"+site_arr[h]+".png")
    print("Mean [137Xe, 64Cu, 66Cu, 20F, 2H] captures per muon, given atleast one capture, at " + site_arr[h] + ": ", mult_df_arr[h][["Xe137Count"]].mean()[0], mult_df_arr[h][["Cu64Count"]].mean()[0],  mult_df_arr[h][["Cu66Count"]].mean()[0], mult_df_arr[h][["F20Count"]].mean()[0], mult_df_arr[h][["H2Count"]].mean()[0], file=out_data_file)
    print("Fraction of muons producing [137Xe, 64Cu, 66Cu, 20F, 2H] at " + site_arr[h] + ": ", len(mult_df_arr[h][["Xe137Count"]].to_numpy(dtype="int64").nonzero()[0])/simulated_muons[h], len(mult_df_arr[h][["Cu64Count"]].to_numpy(dtype="int64").nonzero()[0])/simulated_muons[h],  len(mult_df_arr[h][["Cu66Count"]].to_numpy(dtype="int64").nonzero()[0])/simulated_muons[h], len(mult_df_arr[h][["F20Count"]].to_numpy(dtype="int64").nonzero()[0])/simulated_muons[h], len(mult_df_arr[h][["H2Count"]].to_numpy(dtype="int64").nonzero()[0])/simulated_muons[h], file=out_data_file)
    print(file=out_data_file)

print("Muon tagging efficiency: ", muon_tag_eff, file=out_data_file)

print("Capture tagging efficiency file: ", cap_tag_file, file=out_data_file)
xe137_cap_tag_eff_scan = []
cu64_cap_tag_eff_scan = []
cu66_cap_tag_eff_scan = []
f20_cap_tag_eff_scan = []
h2_cap_tag_eff_scan = []
for scan in range (n_scans):
    ene_min_idx = np.where(bin_centers==n_cap_tag_emin_scan[scan]+20)[0][0]
    ene_max_idx = np.where(bin_centers==n_cap_tag_emax_scan[scan]+20)[0][0]
    xe137_cap_tag_eff_scan.append(cap_tag_eff[0][cap_tag_enevar_idx][cap_tag_vol_idx][ene_min_idx] - cap_tag_eff[0][cap_tag_enevar_idx][cap_tag_vol_idx][ene_max_idx])
    cu64_cap_tag_eff_scan.append(cap_tag_eff[1][cap_tag_enevar_idx][cap_tag_vol_idx][ene_min_idx] - cap_tag_eff[1][cap_tag_enevar_idx][cap_tag_vol_idx][ene_max_idx])
    cu66_cap_tag_eff_scan.append(cap_tag_eff[2][cap_tag_enevar_idx][cap_tag_vol_idx][ene_min_idx] - cap_tag_eff[2][cap_tag_enevar_idx][cap_tag_vol_idx][ene_max_idx])
    f20_cap_tag_eff_scan.append(cap_tag_eff[3][cap_tag_enevar_idx][cap_tag_vol_idx][ene_min_idx] - cap_tag_eff[3][cap_tag_enevar_idx][cap_tag_vol_idx][ene_max_idx])
    h2_cap_tag_eff_scan.append(cap_tag_eff[4][cap_tag_enevar_idx][cap_tag_vol_idx][ene_min_idx] - cap_tag_eff[4][cap_tag_enevar_idx][cap_tag_vol_idx][ene_max_idx])


print("Capture multiplicity threshold even if muon is missed (x events in 10 ms):", veto_mult_thresh, file=out_data_file)

veto_time = 25/(60*24) #veto time after muon [days]
xe137_hl = 3.82/(60*24) # 137Xe half-life [days]
veto_time_eff = 1 - np.power(0.5, veto_time/xe137_hl) #probability of 137Xe decay falling within veto time
print("Veto Time Efficiency:", veto_time_eff, file=out_data_file)

print(file=out_data_file)
print("Run Time [days]: ", run_time, file=out_data_file)
print(file=out_data_file)

for s in range(n_sites): #loop over locations
    print("*******************************", file=out_data_file)
    print("Site: " + site_arr[s], file=out_data_file)
    print("Trials: ", trials[s], file=out_data_file)
    print(file=out_data_file)
    print("Site: " + site_arr[s])
    after_veto_xe137_mean_scan = np.array([])
    after_veto_xe137_std_scan = np.array([])
    deadtime_mean_scan = np.array([])
    deadtime_std_scan = np.array([])

    for scan in range(n_scans):#loop over possible neutron capture tagging efficiencies

        print("Neutron capture window [keV]: [", n_cap_tag_emin_scan[scan], ", ",n_cap_tag_emax_scan[scan], "]")
        xe137_arr=np.array([])
        after_veto_xe137_arr=np.array([])
        deadtime_arr=np.array([])

        xe137_cap_tag_eff = xe137_cap_tag_eff_scan[scan]
        cu64_cap_tag_eff = cu64_cap_tag_eff_scan[scan]
        cu66_cap_tag_eff = cu66_cap_tag_eff_scan[scan]
        f20_cap_tag_eff = f20_cap_tag_eff_scan[scan]
        h2_cap_tag_eff = h2_cap_tag_eff_scan[scan]
        print("Neutron capture energy variable: ", cap_tag_enevar[cap_tag_enevar_idx], file=out_data_file)
        print("Neutron capture volume: ", cap_tag_vol[cap_tag_vol_idx], file=out_data_file)
        print("Neutron capture window [keV]: [", n_cap_tag_emin_scan[scan], ", ",n_cap_tag_emax_scan[scan], "]",  file=out_data_file)
        print("137Xe Tagging efficiency: ", xe137_cap_tag_eff_scan[scan], file=out_data_file)
        print("64Cu Tagging efficiency: ", cu64_cap_tag_eff_scan[scan], file=out_data_file)
        print("66Cu Tagging efficiency: ", cu66_cap_tag_eff_scan[scan], file=out_data_file)
        print("20F Tagging efficiency: ", f20_cap_tag_eff_scan[scan], file=out_data_file)
        print("2H Tagging efficiency: ", h2_cap_tag_eff_scan[scan], file=out_data_file)
        
        for t in range(trials[s]):   #loop over each exposure trial 
            if(t%10 == 0):
                print(f"On trial {t} of {trials[s]}")
                
            xe137 = 0
            cu64 = 0
            cu66 = 0
            f20 = 0
            h2 = 0
            after_veto_xe137 = 0
            deadtime = 0

            num_muons = np.random.poisson(mean_muon_arr[s]*muon_nonzero_cap_frac_arr[s]) #number of muons in exposure that have at least one capture
            
            for mu in range(num_muons): #loop over each muon (with non-zero captures) in exposure
                #if(mu%1000000 == 0):
                #    print(f"On muon {mu} of trial {t} of {trials} for efficiency {xe137_cap_tag_eff}")

                muon_id = np.random.choice(len(mult_arr[s]))
                xe137_mult = mult_arr[s][muon_id][0]
                cu64_mult = mult_arr[s][muon_id][1]
                cu66_mult = mult_arr[s][muon_id][2]
                f20_mult = mult_arr[s][muon_id][3]
                h2_mult = mult_arr[s][muon_id][4]

                xe137 += xe137_mult
                cu64 += cu64_mult
                cu66 += cu66_mult
                f20 += f20_mult
                h2 += h2_mult
                
                total_det_cap = np.random.binomial(xe137_mult, xe137_cap_tag_eff) + np.random.binomial(cu64_mult, cu64_cap_tag_eff) + np.random.binomial(cu66_mult, cu66_cap_tag_eff) + np.random.binomial(f20_mult, f20_cap_tag_eff) + np.random.binomial(h2_mult, h2_cap_tag_eff)
                
                if ((np.random.rand() < muon_tag_eff and total_det_cap > 0) or total_det_cap >= veto_mult_thresh): #veto window is opened
                    deadtime += veto_time
                    for c in range(xe137_mult):
                        if (np.random.rand() > veto_time_eff):
                            after_veto_xe137+=1 #it lived long enough to survive the veto
                else: #no veto window opened, all captures missed
                    after_veto_xe137 += xe137_mult
                            
            xe137_arr=np.append(xe137_arr, xe137)
            after_veto_xe137_arr=np.append(after_veto_xe137_arr, after_veto_xe137)
            deadtime_arr = np.append(deadtime_arr, deadtime)


        print("Total Xe137 [atoms] Mean: ", xe137_arr.mean(), " STD: ", xe137_arr.std(), file=out_data_file)
        print("Unvetoed Xe137 [atoms] Mean: ", after_veto_xe137_arr.mean(), " STD: ", after_veto_xe137_arr.std(), file=out_data_file)
        print("Deadtime [days] Mean: ", deadtime_arr.mean(), " STD: ", deadtime_arr.std(), file=out_data_file)
        after_veto_xe137_mean_scan = np.append(after_veto_xe137_mean_scan, after_veto_xe137_arr.mean())
        after_veto_xe137_std_scan = np.append(after_veto_xe137_std_scan, after_veto_xe137_arr.std())
        deadtime_mean_scan = np.append(deadtime_mean_scan, deadtime_arr.mean())
        deadtime_std_scan = np.append(deadtime_std_scan, deadtime_arr.std())
        
        plt.figure()
        counts, bins = np.histogram(after_veto_xe137_arr, bins=20)
        plt.stairs(counts, bins)
        label = "Xe137 atoms in full LXe after [" + str(n_cap_tag_emin_scan[scan])+ ", "+str(n_cap_tag_emax_scan[scan])+ "] keV veto at "+str(site_arr[s])+" after 10 yr"
        plt.xlabel(label)
        print("Unvetoed 137Xe trials/bin", np.array2string(counts, separator=", "), file=out_data_file)
        print("Unvetoed 137Xe bins", np.array2string(bins, separator=", "), file=out_data_file)
        plt.savefig("plots/veto_toy_mc_2/xe137_untagged_captag_"+"[" + str(n_cap_tag_emin_scan[scan]) + ", "+str(n_cap_tag_emax_scan[scan])+ "]_"+str(site_arr[s])+".png")
        plt.cla()
        plt.clf()
        
        fig2, ax2 = plt.subplots()
        plt.hist(deadtime_arr, bins=20)
        plt.xlabel("Minimum deadtime [days] after"+" [" + str(n_cap_tag_emin_scan[scan])+ ", "+str(n_cap_tag_emax_scan[scan])+ "] keV veto at "+str(site_arr[s])+" after 10 yr")
        plt.savefig("plots/veto_toy_mc_2/deadtime_captag_"+"[" + str(n_cap_tag_emin_scan[scan]) + ", "+str(n_cap_tag_emax_scan[scan])+ "]_"+str(site_arr[s])+".png")
        plt.cla()
        plt.clf()
        
        print(file=out_data_file)
        
    fig3, ax3 = plt.subplots()
    plt.errorbar(xe137_cap_tag_eff_scan, after_veto_xe137_mean_scan, yerr=after_veto_xe137_std_scan, marker="o")
    plt.ylabel("Mean Xe137 atoms after veto at "+str(site_arr[s])+" (10 yr)")
    plt.xlabel("137Xe Neutron Capture Tagging Efficiency")
    plt.savefig("plots/veto_toy_mc_2/xe137_untagged_vs_captag_"+str(site_arr[s])+".png")
    plt.cla()
    plt.clf()
        
    fig4, ax4 = plt.subplots()
    plt.errorbar(xe137_cap_tag_eff_scan, deadtime_mean_scan, yerr=deadtime_std_scan, marker="o")
    plt.ylabel("Minimum deadtime at "+str(site_arr[s])+" (10 yr) [days]")
    plt.xlabel("137Xe Neutron Capture Tagging Efficiency")
    plt.savefig("plots/veto_toy_mc_2/deadtime_vs_captag_"+str(site_arr[s])+".png")
    plt.cla()
    plt.clf()
    plt.close("all")
        
    print(file=out_data_file)
    print("*******************************", file=out_data_file)
    print(file=out_data_file)
#plt.show()

out_data_file.close()   
