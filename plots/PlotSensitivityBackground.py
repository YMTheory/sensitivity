
import ROOT

import cPickle as pickle
import glob, os
import matplotlib.pyplot as plt
import numpy as np

import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'

nEXOlib = '/data/data033/exo/software/nEXO_Sensitivity/quick/v5/lib/libnEXOSensitivity.so'
ROOT.gSystem.Load(nEXOlib)

def ShowPickleResults( inpkl ):

    results = pickle.load(open(inpkl,'rb'))
    print results

def GetBkgdLimit(filename,branch,prob,max=10):

    chain = ROOT.TChain('tree')
    chain.Add(filename)
    chain.SetEstimate(chain.GetEntries()+1)
    branch = branch + '>>hist(%d,%d,%d)'%(10000,max/100.,max)
    chain.Draw(branch,'','goff')
    print 'using', chain.GetSelectedRows()
    hist = ROOT.gDirectory.Get('hist')
    hist.Scale(1./hist.Integral())
    probs = ROOT.TGraph()
    for b in range(1,hist.GetNbinsX()+1):
        probs.SetPoint(probs.GetN(),hist.Integral(1,b),hist.GetBinCenter(b))
    
    return probs.Eval(prob)

def ReadSensFiles( infiles, outpkl, scales = ['1'], livetime = 5.0, mass = 3740, efficiency = 0.82, sens_unit = 1e27 ):

    results = {}
    results['livetime'] = livetime
    results['efficiency'] = efficiency
    results['mass'] = mass
    results['sens_unit'] = sens_unit
    results['scales'] = scales
    
    results['sens'] = {}
    results['median_fwhm'] = {'1t':{},'3t':{}}
    results['ul_90_fwhm'] = {'1t':{},'3t':{}}

    for scale in scales:
        scale_files = infiles % (scale)

        tc = ROOT.TChain( "tree" )
        tc.Add( scale_files ) 

        print "Loaded files for scale", scale, "containing %d toy MCs" % tc.GetEntries()
    
        ## First make plots of the sensitivity distributions (minos and profile)
        tc.Draw("num_signal+num_signal_eHi","stat_sig == 0 && covQual_sig == 3","goff")
        nSel = tc.GetSelectedRows()
        ul_minos, sens_minos = np.zeros(nSel), np.zeros(nSel)
        for n in range( nSel ):
            ul_minos[n] = tc.GetV1()[n]/efficiency
            sens_minos[n] = ROOT.nEXOUtils.GetHalfLife(ul_minos[n],livetime,mass)/sens_unit

        results['sens'][scale] = ROOT.TMath.Median(nSel,np.asarray(sens_minos))
        
        max = 10*float(scale)
        results['median_fwhm']['1t'][scale] = GetBkgdLimit(scale_files,'bkg_fwhm_1t',0.5,max)
        results['median_fwhm']['3t'][scale] = GetBkgdLimit(scale_files,'bkg_fwhm_3t',0.5,max)
        results['ul_90_fwhm']['1t'][scale] = GetBkgdLimit(scale_files,'bkg_fwhm_1t',0.9,max)
        results['ul_90_fwhm']['3t'][scale] = GetBkgdLimit(scale_files,'bkg_fwhm_3t',0.9,max)

        print 'bkgd scale', scale, 'sens', results['sens'][scale], 'median', results['median_fwhm']['3t'][scale]

    pickle.dump(results, open( outpkl, 'wb'))

def MakePlot( inpkl, inpkl_bb2n, outname = 'out_sensbkgds_plot' ):

    axis_label_fontsize = 15
    legend_fontsize = 14
    labels_fontsize = 14

    min_x, max_x = 0.00001, 35
    min_y, max_y = 0, 72.5

    results = pickle.load(open(inpkl,'rb'))
    livetime = results['livetime']
    sens = results['sens']
    sens_unit = results['sens_unit']
    scales = results['scales']        

    bb2n_results = pickle.load(open(inpkl_bb2n,'rb'))
    bb2n_livetime = bb2n_results['livetime']
    bb2n_sens = bb2n_results['sens']
    bb2n_sens_unit = bb2n_results['sens_unit']
    bb2n_scales = bb2n_results['scales']        

    fv_scale = 1.035
    n = len(scales)
    x_meds, x_uls, y_sens = np.zeros(n), np.zeros(n), np.zeros(n)
    x_meds_1t, x_uls_1t = np.zeros(n), np.zeros(n)
    sens_max = 0
    for i, scale in enumerate(scales):
        x_meds[i], x_uls[i], y_sens[i] = results['median_fwhm']['3t'][scale], results['ul_90_fwhm']['3t'][scale], results['sens'][scale]
        x_meds_1t[i], x_uls_1t[i] = results['median_fwhm']['1t'][scale], results['ul_90_fwhm']['1t'][scale]
        if y_sens[i] > sens_max:
            sens_max = y_sens[i]

    bb2n_n = len(bb2n_scales)
    bb2n_x_meds, bb2n_x_uls, bb2n_y_sens = np.zeros(n), np.zeros(n), np.zeros(n)
    bb2n_sens_max = 0
    for i, scale in enumerate(bb2n_scales):
        bb2n_x_meds[i], bb2n_x_uls[i], bb2n_y_sens[i] = bb2n_results['median_fwhm']['3t'][scale], bb2n_results['ul_90_fwhm']['3t'][scale], bb2n_results['sens'][scale]
        if bb2n_y_sens[i] > bb2n_sens_max:
            bb2n_sens_max = bb2n_y_sens[i]

    fig, ax = plt.subplots()
    plt.plot( x_meds , y_sens, 'o-', markerfacecolor='b', markeredgecolor='b', label="Median All Bkgs.")
    plt.plot( x_uls , y_sens, 'o:', markerfacecolor='none', markeredgecolor='b', color='b', label="90% C.L. All Bkgs.")
    plt.plot( bb2n_x_meds , bb2n_y_sens, 'o-', markerfacecolor='r',color='r', markeredgecolor='r', label=r"Median $2\nu\beta\beta$ Excl.")
    plt.plot( bb2n_x_uls , bb2n_y_sens, 'o:', markerfacecolor='none', markeredgecolor='r', color='r', label=r"90% C.L. $2\nu\beta\beta$ Excl.")
    #plt.plot( x_meds_1t , y_sens, 'o-', markerfacecolor='r', markeredgecolor='r', color='r', label="Median FWHM-1t")
    #plt.plot( x_uls_1t , y_sens, 'o:', markerfacecolor='none', markeredgecolor='r', color='r', label="90% C.L. FWHM-1t")

    plt.plot( [min_x,max_x], [bb2n_sens_max,bb2n_sens_max], 'r--', linewidth=1.5, alpha=0.75)
    plt.plot( [min_x,max_x], [sens_max,sens_max], 'b--', linewidth=1.5, alpha=0.75)

    #plt.gca().set_xticks(range(0.01,10))
    plt.plot( [0.03,0.03], [min_y, max_y], 'k--', linewidth=1.5, alpha=0.75)
    plt.plot( [0.41,0.41], [min_y, max_y], 'k--', linewidth=1.5, alpha=0.75)
    plt.plot( [0.99,0.99], [min_y, max_y], 'k--', linewidth=1.5, alpha=0.75)
    plt.plot( [2.1,2.1], [min_y, max_y], 'k--', linewidth=1.5, alpha=0.75)

    plt.plot( [min_x,max_x], [34.2, 34.2], 'g--', linewidth=2, label="Feld.-Cous. Null Bkgd.", alpha=0.75)

    plt.title("nEXO Sensitivity (90% C.L.) in 10 Years")
    plt.xlabel("Backgrounds in FWHM-3t [cts/yr]", fontsize=axis_label_fontsize)
    plt.ylabel(r"$^{136}$Xe  $0\nu\beta\beta$  T$_{1/2}$ [10$^{%d}$yr]"%(np.log10(sens_unit)), fontsize=axis_label_fontsize)
    plt.xlim([min_x, max_x])
    plt.ylim([min_y, max_y])

    step = 5.0
    plt.gca().set_yticks(np.arange(min_y,max_y,step))#range(min_y,max_y+1,2))

    plt.xscale('log')
    #plt.yscale('log')
    plt.legend(loc="lower left",numpoints=1,prop={'size':legend_fontsize})
    fig.set_size_inches(8,5)

    plt.minorticks_on
    ax.grid(True, which='both')
    plt.savefig( outname+'.pdf', bbox_inches='tight' )
    plt.savefig( outname+'.png', bbox_inches='tight' )

    plt.show()

if __name__ == "__main__":

    # ReadSensFiles( "../quick/v5/results/done/fits_db_v73_2016-09-09_0nu_bkgs_scale%s_rdm_10.0_years_0.0_counts_*.root","sens_vs_bkgds_v73_10yr.pkl", livetime = 10.0, scales = ['0.00001','0.0001','0.0002','0.0005','0.001','0.002','0.005','0.010','0.02','0.05','0.100','0.2','0.5','1.000','2.0','5.0','10.000'])
    # ReadSensFiles( "../quick/v5/results/done/fits_db_v73_2016-09-09_0nu_scale%s_not2nu_allbkgs_rdm_10.0_years_0.0_counts_*.root","sens_vs_not2nu_bkgds_v73_10yr.pkl", livetime = 10.0, scales = ['0.00001','0.0001','0.0002','0.0005','0.001','0.002','0.005','0.01','0.02','0.05','0.1','0.2','0.5','1','2','5','10'])
    
    # ShowPickleResults("sens_vs_bkgds_v73_10y.pkl")
    
    MakePlot("sens_vs_bkgds_v73_10yr.pkl","sens_vs_not2nu_bkgds_v73_10yr.pkl","plot_nexo_sens_vs_allbkgds_dbv73_10yr_3t")
