
import ROOT

import cPickle as pickle
import glob, os
import matplotlib.pyplot as plt
import numpy as np
import copy

import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
#matplotlib.rcParams['axes.unicode_minus'] = True

nEXOlib = '../lib/libnEXOSensitivity.so'
ROOT.gSystem.Load(nEXOlib)

def GetVolvsSD():
    gr = ROOT.TGraph()
    gr.SetPoint(0,90,3)
    gr.SetPoint(1,122,2.5)
    gr.SetPoint(2,159,2)
    gr.SetPoint(3,202,1.5)
    gr.SetPoint(4,256,1)
    gr.SetPoint(5,333,0.5)
    for i in range(5):
        gr.Fit('pol3','Q')
    pol = gr.GetFunction('pol3')
    gr.SetMarkerStyle(20)
    c = ROOT.TCanvas()
    gr.Draw('APL')
    #raw_input('continue')
    print pol
    return copy.copy(pol)

def ShowPickleResults( inpkl ):

    results = pickle.load(open(inpkl,'rb'))
    print results

def iterrange( x, nsig=4.0 ):
    ## iterate throwing out outliers a few times and recomputing sigma to make a 
    ## guess at the range
    xcurr = 1.0*x  # make a copy of input array
    for i in range( 3 ):
        cmu, cstd = np.mean(xcurr), np.std(xcurr)
        gpts = np.logical_and( x > cmu-nsig*cstd, x < cmu + nsig*cstd )
        xcurr = x[gpts]
    return [np.max([0,cmu-nsig*cstd]), cmu + nsig*cstd]

def FindCILimits( s, nbins=None, crange=[], CL=0.68 ):

    ## estimate desired bins if not specified
    if( not nbins):
        nbins = np.max([len(s)/10., 10.])
    
    ## iterate range to throw out outliers, if not specified
    if( not crange ):
        crange = iterrange( s )

    hh, be = np.histogram( s, bins=nbins, range=crange )
    bcent = be[:-1] + np.diff(be)/2.0

    ## find median and +/- sigma band from CDF rather than vector directly
    xx = be[:-1]+np.diff(be)/2.
    cdf = np.cumsum(hh)
    cdf = 1.0*cdf/cdf[-1]
    ## alpha is in percent, so divide by an extra 100
    mpos = np.interp([0.5-CL/2., 0.5, 0.5+CL/2.], cdf, xx)

    return mpos

def ReadSensFiles( infiles, outpkl, sdList = [10], livetime = 10.0, efficiency = 0.82, sens_unit = 1e27 ):

    volf = GetVolvsSD()

    volList = {}
    for sd in sdList:
        if sd == 10:
            volList[sd] = 3.86
        else:
            volList[sd] = volf.Eval(sd)

    results = {}
    results['livetime'] = livetime
    results['efficiency'] = efficiency
    results['sens_unit'] = sens_unit
    results['sd_list'] = sdList
    results['fv_list'] = volList
    
    results['sens_cls'] = {}

    for sd in sdList:
        sd_files = infiles % ((640-sd)//10,sd)

        tc = ROOT.TChain( "tree" )
        tc.Add( sd_files ) 

        print "Loaded files for standoff-distance", sd, ", volume", volList[sd], "containing %d toy MCs" % tc.GetEntries()
    
        ## First make plots of the sensitivity distributions (minos and profile)
        tc.Draw("num_signal+num_signal_eHi","stat_sig == 0 && covQual_sig == 3","goff")
        nSel = tc.GetSelectedRows()
        ul_minos, sens_minos = np.zeros(nSel), np.zeros(nSel)
        for n in range( nSel ):
            ul_minos[n] = tc.GetV1()[n]/efficiency
            sens_minos[n] = ROOT.nEXOUtils.GetHalfLife(ul_minos[n],livetime,volList[sd]*1000)/sens_unit
            sens_minos[n] /= 1.034 # finer bins correction

        results['sens_cls'][sd] = FindCILimits(sens_minos)

        print 'volume', volList[sd], 'sens', results['sens_cls'][sd]
        
    pickle.dump(results, open( outpkl, 'wb'))

def EvalCounting( inpkl, outpkl, livetime = 10.0, efficiency = 0.82 ):

    # if already has pkl file from bkgd vs mass plots

    roiList = ['fwhm','1sigma','2sigma','3sigma']
    volList = {'fv':3.86,'3t':3,'2t':2,'1t':1,'3p5t':3.5,'2p5t':2.5,'1p5t':1.5,'0p5t':0.5}
    effList = {'fwhm':0.761,'1sigma':0.683,'2sigma':0.954,'3sigma':0.997}
    ssfrac = 0.94
    clList = sorted([0.5])

    fc = ROOT.TFeldmanCousins(0.9)
    rdm = ROOT.TRandom3()

    fc.SetMuMin(0)
    fc.SetMuMax(300) #(10*mean,10)[mean < 1])
    fc.SetMuStep(0.001)

    results = pickle.load(open(inpkl,'rb'))
    for roi in roiList:
        for vol in volList:
            for cl in clList:
                bkgyr = results[roi][vol][cl]
                mean = livetime*bkgyr # linear approximation better than 0.1%
                
                #fc.SetMuStep((fc.GetMuMax()-fc.GetMuMin())/10000.)

                #draws = []
                #for draw in range(1000):
                #    draws.append(fc.CalculateUpperLimit(rdm.Poisson(mean),mean))
                #signal = ROOT.TMath.Median(len(draws),np.asarray(draws))

                signal = 0
                obs = 0
                diff = 1
                while diff > 0.01 or obs < mean:
                    prev = signal
                    fcul = fc.CalculateUpperLimit(obs,mean)
                    pmo = ROOT.TMath.PoissonI(obs,mean)
                    signal += pmo * fcul
                    diff = signal - prev
                    obs += 1
                    #print mean, obs, pmo, fcul, signal

                results[roi][vol][cl] = ROOT.nEXOUtils.GetHalfLife(signal,livetime,volList[vol]*1000)
                results[roi][vol][cl] *= (efficiency*effList[roi]*ssfrac)
                print roi, vol, cl, mean, signal, results[roi][vol][cl]
               # raw_input('continue')
    
    pickle.dump(results, open( outpkl, 'wb'))


def MakePlot( inpkl, outname, countpkl = None ):

    axis_label_fontsize = 15
    legend_fontsize = 14
    labels_fontsize = 14

    min_x, max_x = 0.25, 4.25
    min_y, max_y = 1.1, 9.9

    results = pickle.load(open(inpkl,'rb'))
    fvList = results['fv_list']
    sens = results['sens_cls']
    sens_unit = results['sens_unit']
    badSd = [40,120]

    xList = [fvList[sd] for sd in sorted(fvList) if not sd in badSd]
    yList = [sens[sd][1] for sd in sorted(sens) if not sd in badSd]

    if countpkl:
        #roiList = ['fwhm','1sigma','2sigma'] #,'3sigma']
        #volList = {'fv':3.74,'3t':3,'2t':2,'1t':1,'3p5t':3.5,'2p5t':2.5,'1p5t':1.5,'0p5t':0.5}
        #roiList = {'fwhm':'FWHM (2428 - 2488 keV)','1sigma':'1$\sigma$ (2433 - 2483 keV)','2sigma':'2$\sigma$ (2408 - 2507 keV)','3sigma':'3$\sigma$ (2384 - 2532 keV)'}
        # roiList = {'fwhm':'FWHM','1sigma':'$\pm$ 1$\sigma$','2sigma':'$\pm$ 2$\sigma$'}#,'3sigma':'\pm 3$\sigma$'}
        roiList = {'fwhm':'FWHM',}
        volList = {0.5:'0p5t',1:'1t',1.5:'1p5t',2:'2t',2.5:'2p5t',3:'3t',3.86:'fv'}#,3.5:'3p5t'
        markerList = {'fwhm':'bD-','1sigma':'ro-','2sigma':'g^-','3sigma':'m*-'}
        orderList = ['fwhm', ] #'3sigma'

        count_results = pickle.load(open(countpkl,'rb'))
        count_xList = [vol[0] for vol in sorted(volList.items())]
        count_yList = {}
        for roi in roiList:
            count_yList[roi] = [count_results[roi][vol[1]][0.5]/sens_unit for vol in sorted(volList.items())]   
    # fix fv
    xList[0] = 3.86

    fig, ax = plt.subplots()
    plt.plot( xList, yList, 'ks-', label='Full nEXO 2D Analysis')
    # if countpkl:
    #     for roi in orderList:
    #         plt.plot( count_xList, count_yList[roi], markerList[roi], label='Count. Exp. Feld.'+u'\u2212'+'Cous. ('+roiList[roi]+')')

    plt.plot( count_xList, count_yList['fwhm'], markerList['fwhm'], label='Counting Experiment (FWHM)')

    #plt.title("nEXO Sensitivity (90% C.L.) in 10 Years")
    plt.xlabel("Fiducial Volume Mass [tonne]", fontsize=axis_label_fontsize)
    plt.ylabel(r"$^{136}$Xe  $0\nu\beta\beta$  T$_{1/2}$ [$10^{%d}$ yr]"%(np.log10(sens_unit)), fontsize=axis_label_fontsize)
    plt.xlim([min_x, max_x])
    plt.ylim([min_y, max_y])
    
    #plt.yscale('log')

    plt.legend(loc="lower right",prop={'size':legend_fontsize}) #,frameon=False
    ax.get_legend().get_title().set_fontsize(legend_fontsize)

    ax.grid(True, which='both')
    fig.set_size_inches(7.5,5.5)

    plt.savefig("pdf/" + outname+'.pdf', bbox_inches='tight' )
    plt.savefig("png/" + outname+'.png', bbox_inches='tight' )
    # plt.show()


if __name__ == "__main__":

    # ReadSensFiles("../quick/v5/results/done/fits_db_v73_2016-09-09_0nu_bkgs_rdm_bin270_800_3500_%d_%d_640_10.0_years_0.0_counts_*.root","sens_vs_vol_new.pkl",sdList=[10,40,90,120,160,200,260,330])
    # EvalCounting( "bkgd_vs_mass.pkl", "count_sens_fvs_new.pkl", livetime = 10.0, efficiency = 0.82 )
    MakePlot("sens_vs_vol.pkl","plot_sens_vs_fv_v1", countpkl = "count_sens_fvs.pkl")
