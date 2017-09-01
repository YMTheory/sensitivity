
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

def GetBkgdLimit(filename,branch,prob):

    chain = ROOT.TChain('tree')
    chain.Add(filename)
    chain.SetEstimate(chain.GetEntries()+1)
    branch = branch + '>>hist(10000,0,10)'
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
        
        results['median_fwhm']['1t'][scale] = GetBkgdLimit(scale_files,'bkg_fwhm_1t',0.5)
        results['median_fwhm']['3t'][scale] = GetBkgdLimit(scale_files,'bkg_fwhm_3t',0.5)
        results['ul_90_fwhm']['1t'][scale] = GetBkgdLimit(scale_files,'bkg_fwhm_1t',0.9)
        results['ul_90_fwhm']['3t'][scale] = GetBkgdLimit(scale_files,'bkg_fwhm_3t',0.9)

        print 'bkgd scale', scale, 'sens', results['sens'][scale], 'median', results['median_fwhm']['3t'][scale]

    pickle.dump(results, open( outpkl, 'wb'))

def MakePlot( inpkl, outname = 'out_sensbkgds_plot.pdf' ):

    axis_label_fontsize = 15
    legend_fontsize = 14
    labels_fontsize = 14

    min_y, max_y = 0, 30

    results = pickle.load(open(inpkl,'rb'))
    livetime = results['livetime']
    sens = results['sens']
    sens_unit = results['sens_unit']
    scales = results['scales']        

    fv_scale = 1.035
    n = len(scales)
    x_meds, x_uls, y_sens = np.zeros(n), np.zeros(n), np.zeros(n)
    #x_meds_1t, x_uls_1t = np.zeros(n), np.zeros(n)
    for i, scale in enumerate(scales):
        x_meds[i], x_uls[i], y_sens[i] = results['median_fwhm']['3t'][scale], results['ul_90_fwhm']['3t'][scale], results['sens'][scale]
        #x_meds_1t[i], x_uls_1t[i] = results['median_fwhm']['1t'][scale], results['ul_90_fwhm']['1t'][scale]

    fig, ax = plt.subplots()
    plt.plot( x_meds , y_sens, 'o-', markerfacecolor='b', markeredgecolor='b', label="Median")
    plt.plot( x_uls , y_sens, 'o-', markerfacecolor='none', markeredgecolor='r', color='r', label="90% C.L.")
    #plt.plot( x_meds_1t , y_sens, 'ro', linewidth=2, label="Median")

    #plt.gca().set_xticks(range(0.01,10))
    plt.plot( [0.03,0.03], [min_y, max_y], 'g--', linewidth=1)
    plt.plot( [0.6,0.6], [min_y, max_y], 'g--', linewidth=1)
    plt.plot( [0.99,0.99], [min_y, max_y], 'g--', linewidth=1)
    plt.plot( [2.1,2.1], [min_y, max_y], 'g--', linewidth=1)

    plt.xlabel("Background FWHM-3t [cts/yr]", fontsize=axis_label_fontsize)
    plt.ylabel(r"$^{136}$Xe  $0\nu\beta\beta$  T$_{1/2}$ [10$^{%d}$yr]"%(np.log10(sens_unit)), fontsize=axis_label_fontsize)

    plt.xlim([0.01,10])
    plt.ylim([min_y, max_y])
    plt.xscale('log')
    plt.legend(loc="upper right",numpoints=1,prop={'size':legend_fontsize})
    fig.set_size_inches(6,5)

    plt.minorticks_on
    ax.grid(True, which='both')
    plt.savefig( outname )

    plt.show()
