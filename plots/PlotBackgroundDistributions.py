
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

def fill_between(ax, x, y1, y2=0, **kwargs):
    """ version of fill_between that adds an empty patch
        to make it easier to add to the legend """
    ax = plt.gca()
    ax.fill_between(x, y1, y2, **kwargs)
    #p = ax.Rectangle((1e-20, 1e-20), 1e-20, 1e-20, **kwargs)
    #ax.add_patch(p)
    #return p

def iterrange( x, nsig=4.0 ):
    ## iterate throwing out outliers a few times and recomputing sigma to make a 
    ## guess at the range
    xcurr = 1.0*x  # make a copy of input array
    for i in range( 3 ):
        cmu, cstd = np.mean(xcurr), np.std(xcurr)
        gpts = np.logical_and( x > cmu-nsig*cstd, x < cmu + nsig*cstd )
        xcurr = x[gpts]
    return [np.max([0,cmu-nsig*cstd]), cmu + nsig*cstd]

def ShowPickleResults( inpkl ):

    results = pickle.load(open(inpkl,'rb'))
    print results

def ReadFiles( infiles, outpkl):

    results = {}
    
    chain = ROOT.TChain('tree')
    chain.Add( infiles )

    bkgdList = ['bkg_fwhm_1t','bkg_fwhm_3t','bkg_fwhm_fv']
    
    bkgdString = ''
    for bkgd in bkgdList:
        bkgdString += '%s:'%(bkgd)
    bkgdString = bkgdString[:-1]

    chain.Draw(bkgdString,'','goff')
    nSel = chain.GetSelectedRows()
    for b,bkgd in enumerate(bkgdList):
        results[bkgd] = np.zeros(nSel)
        for i in range(nSel):
            results[bkgd][i] = chain.GetVal(b)[i]
            
    pickle.dump(results, open( outpkl, 'wb'))    
    
def PlotDistribution( s, ax, nbins = None, crange = [] ):

    legend_fontsize = 14

    ## estimate desired bins if not specified
    if( not nbins):
        nbins = np.max([len(s)/200., 10.])
    
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
    mpos = np.interp([0, 0.5, 0.9], cdf, xx)

    gpts = np.logical_and( bcent > mpos[0], bcent < mpos[2] )
    xx1,xx2,yy = be[gpts], be[1:][gpts], hh[gpts]
    xx = np.ndarray.flatten(np.vstack((xx1,xx2)).T)
    yy = np.ndarray.flatten(np.vstack((yy,yy)).T)

    #fig = plt.figure()
    #fill_between(ax, xx, 0, yy, color='c', edgecolor='None', alpha=0.5, label="90% Conf. Level")
    ax.fill_between(xx, 0, yy, color='c', edgecolor='None', alpha=0.5, label="90% Interval")
    ax.step(be[:-1], hh, where='post', color='k', linewidth=1.5)
    med_val = hh[ np.argmin( np.abs( bcent - mpos[1]) )]
#    plt.plot( [mpos[1], mpos[1]], [0, med_val], 'r', linewidth=1.5, label="Median = %.1f " % mpos[1])
    label_unit = ("Median = %.2f" % mpos[1],"Median = %.1f" % mpos[1])[mpos[1] > 1]
    ax.plot( [mpos[1],mpos[1]], [0, med_val], 'r', linewidth=1.5, label=label_unit)

    ax.legend(loc="upper center",prop={'size':legend_fontsize})

    #plt.xlabel(xlab, fontsize=axis_label_fontsize)
    #plt.ylabel(ylab, fontsize=axis_label_fontsize)
    #if( plot_title ):
    #    plt.title(plot_title)
    
    #return fig


def MakePlots(inpkl, outname):

    results = pickle.load(open(inpkl,'rb'))

    axis_label_fontsize = 15
    fig, (ax_fv, ax_3t, ax_1t) = plt.subplots(1, 3, sharey=True, squeeze=True)

    PlotDistribution(results['bkg_fwhm_fv'],ax_fv)
    PlotDistribution(results['bkg_fwhm_3t'],ax_3t)
    PlotDistribution(results['bkg_fwhm_1t'],ax_1t)
    
    fig.subplots_adjust(wspace = 0, hspace = 0)

    fig.set_size_inches(12,4.5)

    plt.ylim([0,950])


    ax_fv.set_ylabel("Number of Toy Fits", fontsize=axis_label_fontsize)
    ax_fv.set_xlabel("Bkgd. FWHM-FV [cts/yr]", fontsize=axis_label_fontsize)
    ax_fv.locator_params(axis='x',nbins=7)
    ax_fv.get_xticklabels()[0].set_visible(False)
    ax_fv.get_xticklabels()[-1].set_visible(False)
    ax_3t.set_xlabel("Bkgd. FWHM-3t [cts/yr]", fontsize=axis_label_fontsize)
    ax_3t.locator_params(axis='x',nbins=9)
    ax_3t.get_xticklabels()[0].set_visible(False)
    ax_3t.get_xticklabels()[-1].set_visible(False)
    ax_1t.set_xlabel("Bkgd. FWHM-1t [cts/yr]", fontsize=axis_label_fontsize)
    ax_1t.locator_params(axis='x',nbins=7)
    ax_1t.get_xticklabels()[0].set_visible(False)
    ax_1t.get_xticklabels()[-1].set_visible(False)

    plt.savefig( outname + '.pdf', bbox_inches='tight' )
    plt.savefig( outname + '.png', bbox_inches='tight' )

    plt.show()



if __name__ == "__main__":

    # ReadFiles("/data/data033/exo/software/nEXO_Sensitivity/quick/v5/results/done/fits_db_v73_2016-09-09_0nu_rdm_5.0_years_0.0_counts_*.root","v73_bkg_dist.pkl")
    # ShowPickleResults("v73_bkg_dist.pkl")
    MakePlots("v73_bkg_dist.pkl", "v73_bkg_dist_v1")
