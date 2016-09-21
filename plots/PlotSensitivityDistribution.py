
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

def fill_between(x, y1, y2=0, **kwargs):
    """ version of fill_between that adds an empty patch
        to make it easier to add to the legend """
    ax = plt.gca()
    plt.fill_between(x, y1, y2, **kwargs)
    p = plt.Rectangle((1e-20, 1e-20), 1e-20, 1e-20, **kwargs)
    #ax.add_patch(p)
    return p

def iterrange( x, nsig=4.0 ):
    ## iterate throwing out outliers a few times and recomputing sigma to make a 
    ## guess at the range
    xcurr = 1.0*x  # make a copy of input array
    for i in range( 3 ):
        cmu, cstd = np.mean(xcurr), np.std(xcurr)
        gpts = np.logical_and( x > cmu-nsig*cstd, x < cmu + nsig*cstd )
        xcurr = x[gpts]
    return [np.max([0,cmu-nsig*cstd]), cmu + nsig*cstd]

def PlotSensitivity( s, median, unit = 1, nbins=None, crange=[], CL=0.68, xlab="90% CL upper limit", ylab="Number of Toy Fits", plot_title="" ):

    axis_label_fontsize = 15
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
    mpos = np.interp([0.5-CL/2., 0.5, 0.5+CL/2.], cdf, xx)

    gpts = np.logical_and( bcent > mpos[0], bcent < mpos[2] )
    xx1,xx2,yy = be[gpts], be[1:][gpts], hh[gpts]
    xx = np.ndarray.flatten(np.vstack((xx1,xx2)).T)
    yy = np.ndarray.flatten(np.vstack((yy,yy)).T)

    fig = plt.figure()
    fill_between(xx, 0, yy, color='c', edgecolor='None', alpha=0.5, label="68% Confidence Interval")
    plt.step(be[:-1], hh, where='post', color='k', linewidth=1.5)
    med_val = hh[ np.argmin( np.abs( bcent - mpos[1]) )]
#    plt.plot( [mpos[1], mpos[1]], [0, med_val], 'r', linewidth=1.5, label="Median = %.1f " % mpos[1])
    label_unit = "Median = %.1f" % median
    if unit != 1:
        label_unit = r"%s $\times$ 10$^{%d}$"%(label_unit, np.log10(unit))
    plt.plot( [median, median], [0, med_val], 'r', linewidth=1.5, label=label_unit)

    plt.legend(loc="upper right",prop={'size':legend_fontsize})
    plt.xlabel(xlab, fontsize=axis_label_fontsize)
    plt.ylabel(ylab, fontsize=axis_label_fontsize)
    if( plot_title ):
        plt.title(plot_title)
    
    return fig


def MakePlots(inpkl, outdir):

    results = pickle.load(open(inpkl,'rb'))
    ul_minos = results['ul_minos']
    sens_minos = results['sens_minos']
    livetime = results['livetime']
    ul_minos_median = results['ul_minos_median']
    sens_minos_median = results['sens_minos_median']
    sens_unit = results['sens_unit']

    sens_fig_minos = PlotSensitivity( ul_minos, ul_minos_median, plot_title="Upper Limit Distribution in %.1f Years"%(livetime), xlab=r"$0\nu\beta\beta$ Upper Limit (90% CL) [cts]" )
    plt.savefig( os.path.join( outdir, "UpperLimitDistribution.png" ) )
    
    sens_fig_prof = PlotSensitivity( sens_minos, sens_minos_median, unit=sens_unit, plot_title="Sensitivity Distribution in %.1f Years"%(livetime), xlab=r"$^{136}$Xe  $0\nu\beta\beta$  T$_{1/2}$ [10$^{27}$ yr]" )
    plt.savefig( os.path.join( outdir, "SensitivityDistribution.png" ) )

    plt.show()

def ReadFiles( livetime, mass, efficiency, infiles, outpkl, sens_unit = 1e27 ):
    
    tc = ROOT.TChain( "tree" )
    tc.Add( infiles ) #os.path.join( inpath, "*.root") )

    print "Loaded files containing %d toy MCs" % tc.GetEntries()
    
    ## First make plots of the sensitivity distributions (minos and profile)
    tc.Draw("num_signal+num_signal_eHi","stat_sig == 0 && covQual_sig == 3","goff")
    nSel = tc.GetSelectedRows()
    ul_minos, sens_minos = np.zeros(nSel), np.zeros(nSel)
    for n in range( nSel ):
        ul_minos[n] = tc.GetV1()[n]/efficiency
        sens_minos[n] = ROOT.nEXOUtils.GetHalfLife(ul_minos[n],livetime,mass)/sens_unit
    ul_minos_median = ROOT.TMath.Median(nSel,np.asarray(ul_minos))
    sens_minos_median = ROOT.nEXOUtils.GetHalfLife(ul_minos_median,livetime,mass)/sens_unit
    print ul_minos_median, 'cts', sens_minos_median, 'x 10^', np.log10(sens_unit), 'yr'

    results = {}
    results['livetime'] = livetime
    results['efficiency'] = efficiency
    results['mass'] = mass
    results['ul_minos'] = ul_minos
    results['sens_minos'] = sens_minos
    results['ul_minos_median'] = ul_minos_median
    results['sens_minos_median'] = sens_minos_median
    results['sens_unit'] = sens_unit
        
    pickle.dump(results, open( outpkl, 'wb'))

if __name__ == "__main__":
    from optparse import OptionParser

    usage = '''%prog -l livetime [yr] -m mass [kg] -e efficiency -i infiles -o outdir
            '''

    parser = OptionParser(usage)
    parser.add_option("-l","--livetime", type=float, nargs=1, default=5.0)
    parser.add_option("-m","--mass", type=float, nargs=1, default=3740)
    parser.add_option("-e","--efficiency", type=float, nargs=1, default=1.00)

    parser.add_option("-i","--infiles", default=None)
    parser.add_option("-o","--outdir", default=".")

    options,args = parser.parse_args()
    if options.infiles == None:
        print '*** ERROR: Must specify input file path ***'
        parser.print_help()
        sys.exit(1)

    MakePlots( options.livetime, options.mass, options.efficiency, options.infiles, options.outdir)
