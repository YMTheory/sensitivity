
import ROOT

import cPickle as pickle
import glob, os
import matplotlib.pyplot as plt
import numpy as np
import array

import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'

#nEXOlib = '/data/data033/exo/software/nEXO_Sensitivity/quick/v5/lib/libnEXOSensitivity.so'
#ROOT.gSystem.Load(nEXOlib)

def ShowPickleResults( inpkl ):

    results = pickle.load(open(inpkl,'rb'))
    print results

def GetBkgdLimit(filename,branch,prob,max=10):

    chain = ROOT.TChain('tree')
    chain.Add(filename)
    chain.SetEstimate(chain.GetEntries()+1)
    branch = branch + '>>hist(%d,%d,%d)'%(100000,max/1000.,max)
    chain.Draw(branch,'','goff')
    print 'using', chain.GetSelectedRows()
    hist = ROOT.gDirectory.Get('hist')
    hist.Scale(1./hist.Integral())
    probs = ROOT.TGraph()
    for b in range(1,hist.GetNbinsX()+1):
        probs.SetPoint(probs.GetN(),hist.Integral(1,b),hist.GetBinCenter(b))
    
    return probs.Eval(prob)

def GetHistPercent(nVals,valArray,oneArray,probs,nbins=1000):

    h = ROOT.TH1D('h','',nbins,ROOT.TMath.MinElement(nVals,valArray),ROOT.TMath.MaxElement(nVals,valArray))
    h.FillN(nVals,valArray,oneArray)
   
    nq = len(probs)
    xq = np.asarray(probs)        
    yq = np.asarray(np.zeros(nq))

    h.GetQuantiles(nq,yq,xq)

    return yq
    

def ReadFiles( infiles, outpkl ):

    results = {}

    roiList = ['fwhm','1sigma','2sigma','3sigma']
    volList = ['fv','3t','2t','1t','3p5t','2p5t','1p5t','0p5t']
    clList = sorted([0.5,0.9])

    idx, idxDict = 1, {}    
    drawString = '1'
    for roi in roiList:
        idxDict[roi] = {}
        for vol in volList:
            idxDict[roi][vol] = idx
            idx += 1
            drawString += ':bkg_%s_%s'%(roi,vol)
    #drawString = drawString[:-1]

    chain = ROOT.TChain('tree')
    chain.Add(infiles)
    chain.SetEstimate(chain.GetEntries()+1)
    chain.Draw(drawString,'','para goff')

    for roi in roiList:
        results[roi] = {}
        for vol in volList:
            results[roi][vol] = {}
            results[roi][vol]['median'] = ROOT.TMath.Median(chain.GetSelectedRows(),chain.GetVal(idxDict[roi][vol]))
            print roi, vol, 'median', results[roi][vol]['median']
            percList = GetHistPercent(chain.GetSelectedRows(),chain.GetVal(idxDict[roi][vol]),chain.GetVal(0),clList)
            for cl,perc in zip(clList,percList):
                results[roi][vol][cl] = perc
                print roi, vol, cl, results[roi][vol][cl]

    pickle.dump(results, open( outpkl, 'wb'))

def MakePlot( inpkl, outname, logscale=False ):

    axis_label_fontsize = 15
    legend_fontsize = 14
    labels_fontsize = 14

    min_x, max_x = 0, 4
    min_y, max_y = 0, 3 
    if logscale:
        min_y, max_y = 0.02, 5 # 0, 3 #11 #8.5

    roiList = {'fwhm':'FWHM (2428 - 2488 keV)','1sigma':'$\pm$1$\sigma$ (2433 - 2483 keV)','2sigma':'$\pm$2$\sigma$ (2408 - 2507 keV)','3sigma':'$\pm$3$\sigma$ (2384 - 2532 keV)'}
    volList = {0.5:'0p5t',1:'1t',1.5:'1p5t',2:'2t',2.5:'2p5t',3:'3t',3.8:'fv'}#,3.5:'3p5t'}#,3.8:'fv'}
    clList = sorted([0.5])
    orderList = ['1sigma','fwhm','2sigma']#,'3sigma']

    #markerList = {'fwhm':'gD-','1sigma':'ro-','2sigma':'bs-','3sigma':'m^-'}
    markerList = {'fwhm':'bo-','1sigma':'gD-','2sigma':'rs-','3sigma':'m^-'}

    results = pickle.load(open(inpkl,'rb'))

    xList = [vol[0] for vol in sorted(volList.items())]
    print xList
    yList = {}
    for roi in roiList:
        yList[roi] = {}
        for cl in clList:
            yList[roi][cl] = [results[roi][vol[1]][cl]/vol[0] for vol in sorted(volList.items())]

    fig, ax = plt.subplots()
    for roi in orderList:
        plt.plot( xList, yList[roi][0.5], markerList[roi], label=roiList[roi])

    #plt.title("nEXO Sensitivity (90% C.L.) in 10 Years")
    plt.xlabel("Liquid Xe Mass [tonne]", fontsize=axis_label_fontsize)
    plt.ylabel("Bgd. Index [cts/(ROI$\cdot$tonne$\cdot$yr)]", fontsize=axis_label_fontsize)
    plt.xlim([min_x, max_x])
    plt.ylim([min_y, max_y])
    
    if logscale:
        plt.yscale('log')

    plt.legend(loc=("upper left","lower right")[logscale],title="Region of Interest",prop={'size':legend_fontsize}) #,frameon=False
    ax.get_legend().get_title().set_fontsize(legend_fontsize)

    step_y = 0.25
    if not logscale:
        plt.gca().set_yticks(np.arange(min_y, max_y+step_y,step_y))
    ax.grid(True, which='both')
    fig.set_size_inches(7.5,5.5)

    plt.savefig( outname+'.pdf', bbox_inches='tight' )
    plt.savefig( outname+'.png', bbox_inches='tight' )
    plt.show()



if __name__ == "__main__":

    # ReadFiles("../quick/v5/results/done/fits_db_v73_2016-09-09_0nu_scale1_not2nu_allbkgs_rdm_10.0_years_0.0_counts_*.root","bkgd_vs_mass.pkl")
    # ReadFiles("../quick/v5/results/done/fits_db_v73b_2016-09-28_0nu_plots_rdm_10.0_years_0.0_counts_*.root","bkgd_vs_mass_2nu.pkl")

    MakePlot("bkgd_vs_mass_2nu.pkl","plot_bkgd_vs_mass_v5_lin",logscale=False)
