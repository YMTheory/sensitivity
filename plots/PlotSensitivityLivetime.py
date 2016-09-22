
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

def fill_between(x, y1, y2=0, **kwargs):
    """Plot filled region between `y1` and `y2`.

    This function works exactly the same as matplotlib's fill_between, except
    that it also plots a proxy artist (specifically, a rectangle of 0 size)
    so that it can be added it appears on a legend.
    """
    ax = plt.gca()
    plt.fill_between(x, y1, y2, **kwargs)
    #p = plt.Rectangle((1e-20, 1e-20), 1e-20, 1e-20, **kwargs)
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
    

def PowerFit( data, idx = None ):
    # get x,y from data [idx = lo, med, hi]
    # fit y = p0 * x^p1
    # return p0 and p1

    xs, ys = [], []
    n = len(data)
    for livetime in sorted(data):
        xs.append(livetime)
        if not idx is None:
            ys.append(data[livetime][idx])
        else:
            ys.append(data[livetime])
    print n, xs, ys

    f = ROOT.TF1('f','[0]*pow(x,[1])',xs[0],xs[-1])
    f.SetParameters(ys[0],1)
    
    g = ROOT.TGraph(len(data), np.asarray(xs), np.asarray(ys))
    for i in range(10):
        g.Fit(f,'QW')
    g.Fit(f)

    g.Draw('AP*')
    if not idx is None:
        raw_input('continue idx = %d'%(idx))
    else:
        raw_input('continue')

    return f.GetParameter(0), f.GetParameter(1)


def ReadSensFiles( infiles, outpkl, livetimes = [5.0], mass = 3740, efficiency = 0.82, sens_unit = 1e27, create_disc = False ):

    results = {}
    results['livetimes'] = livetimes
    results['efficiency'] = efficiency
    results['mass'] = mass
    results['sens_unit'] = sens_unit
    
    results['ul_minos_cls'] = {}
    results['sens_minos_cls'] = {}

    for livetime in livetimes:
        livetime_files = infiles % (livetime)

        tc = ROOT.TChain( "tree" )
        tc.Add( livetime_files ) 

        print "Loaded files for livetime", livetime, "containing %d toy MCs" % tc.GetEntries()
    
        ## First make plots of the sensitivity distributions (minos and profile)
        tc.Draw("num_signal+num_signal_eHi","stat_sig == 0 && covQual_sig == 3","goff")
        nSel = tc.GetSelectedRows()
        ul_minos, sens_minos = np.zeros(nSel), np.zeros(nSel)
        for n in range( nSel ):
            ul_minos[n] = tc.GetV1()[n]/efficiency
            sens_minos[n] = ROOT.nEXOUtils.GetHalfLife(ul_minos[n],livetime,mass)/sens_unit

        results['ul_minos_cls'][livetime] = FindCILimits(ul_minos)
        results['sens_minos_cls'][livetime] = FindCILimits(sens_minos)

        print 'livetime', livetime, 'ul', results['ul_minos_cls'][livetime] , 'sens', results['sens_minos_cls'][livetime]

    results['ul_minos_lo_pow'] = PowerFit(results['sens_minos_cls'],0)
    results['ul_minos_med_pow'] = PowerFit(results['sens_minos_cls'],1)
    results['ul_minos_hi_pow'] = PowerFit(results['sens_minos_cls'],2)
        
    if create_disc:
        results['disc_minos'] = {}
        ba_sens_to_disc_livetimes = [0.1,0.5,1.0,2,3,4,5]
        ba_sens_to_disc_ratios = [.95,0.70,0.60,.50,.50,.50,.50]
        ba_sens_to_disc_conversion = ROOT.TGraph(len(ba_sens_to_disc_ratios),np.asarray(ba_sens_to_disc_livetimes),np.asarray(ba_sens_to_disc_ratios))
        for livetime in results['sens_minos_cls']:
            results['disc_minos'][livetime] = ba_sens_to_disc_conversion.Eval(livetime)*results['sens_minos_cls'][livetime][1]
        results['disc_minos_pow'] = PowerFit(results['disc_minos'])        

    pickle.dump(results, open( outpkl, 'wb'))

def ReadDiscFiles( infiles, outpkl, livetimes = [5.0], signals = [0.0], mass = 3740, sens_unit = 1e27, prob = 0.5, sigma = 3 ):

    results = {}
    results['livetimes'] = livetimes
    results['signals'] = signals
    results['mass'] = mass
    results['sens_unit'] = sens_unit
    results['probability'] = prob
    results['significance'] = sigma    
    
    results['disc_minos'] = {}

    for livetime in sorted(livetimes):
        results['disc_minos'][livetime] = ROOT.nEXOUtils.GetDiscHalfLife(prob,sigma,len(signals),np.asarray(signals),livetime,infiles,mass) 
        results['disc_minos'][livetime] /= sens_unit
        results['disc_minos'][livetime] *= 2.1 # factor of 2 is to scale old result with new ones, must disappear!!!
        print 'livetime', livetime, 'discovery', results['disc_minos'][livetime]

    results['disc_minos_pow'] = PowerFit(results['disc_minos'])
        
    pickle.dump(results, open( outpkl, 'wb'))


def MakePlot( inpkl, outname = 'out_sens_plot', discpkl = None, imprpkl = None, labels = False, exo200 = True):

    axis_label_fontsize = 15
    legend_fontsize = 14
    labels_fontsize = 14

    results = pickle.load(open(inpkl,'rb'))
    livetimes = results['livetimes']
    #sens_minos_cls = results['sens_minos_cls']
    sens_unit = results['sens_unit']
    median_pow = results['ul_minos_med_pow']
    lolim_pow = results['ul_minos_lo_pow']
    hilim_pow = results['ul_minos_hi_pow']

    if discpkl:
        disc_results = pickle.load(open(discpkl,'rb'))
        disc_pow = disc_results['disc_minos_pow']

    if imprpkl:
        impr_results = pickle.load(open(imprpkl,'rb'))
        impr_sens_unit = impr_results['sens_unit']
        impr_median_pow = impr_results['ul_minos_med_pow']
        

    fv_scale = 1.035
    time_scale = 0.0002
    xlist = np.arange(time_scale,livetimes[-1]+time_scale,time_scale)
    ylist, ylolist, yhilist, dlist, imprlist = [], [], [], [], []
    for x in xlist:
        ylist.append(median_pow[0]*(x**median_pow[1])*sens_unit*fv_scale)
        ylolist.append(lolim_pow[0]*(x**lolim_pow[1])*sens_unit*fv_scale)
        yhilist.append(hilim_pow[0]*(x**hilim_pow[1])*sens_unit*fv_scale)
        if discpkl:
            dlist.append(disc_pow[0]*(x**disc_pow[1])*sens_unit*fv_scale)
        if imprpkl:
            imprlist.append(impr_median_pow[0]*(x**impr_median_pow[1])*impr_sens_unit*fv_scale)

    xfom = 5.0 # yr
    yfom = median_pow[0]*(xfom**median_pow[1])*sens_unit*fv_scale

    y200 = [1.9e25/sens_unit]
    x200 = (np.log(y200[0]) - np.log(median_pow[0]))/median_pow[1]
    x200 = [np.exp(x200)]
    y200 = [y200[0]*sens_unit]
    print x200, y200
    #fig = plt.figure()
    fig, ax = plt.subplots()
    if imprpkl:
        plt.plot( xlist, imprlist, linestyle=':', color='orange', linewidth=2, label="nEXO Improved Sensitivty (90% C.L.)")
    plt.plot( xlist, ylist, 'r-', markersize=4, markerfacecolor='r', linewidth=2, label="nEXO Sensitivity (90% C.L.)")
    #fill_between( xlist, ylolist, yhilist, color='c', edgecolor='None', alpha=0.5, label="68% Confidence Interval" ) 
    if discpkl:
        plt.plot( xlist, dlist, 'b--', markersize=4, markerfacecolor='b', linewidth=2, label=r"nEXO Discovery 3$\sigma$, Prob. 50%")
    if exo200:
        plt.plot( x200, y200, 'o', markersize=10, markeredgewidth=2, markerfacecolor='none', label="EXO-200 Sensitivity (90% C.L.)", markeredgecolor='green')

    plt.gca().set_xticks(range(0,livetimes[-1]+1))

    plt.xlabel("Livetime [yr]", fontsize=axis_label_fontsize)
    #plt.ylabel("T$_{1/2}$ [10$^{%d}$ yr]"%(np.log10(sens_unit)), fontsize=16)
    plt.ylabel(r"$^{136}$Xe  $0\nu\beta\beta$  T$_{1/2}$ [yr]", fontsize=axis_label_fontsize)

    plt.xlim([-0.2, livetimes[-1]])
    #plt.ylim([(0.1*sens_unit,0.01*sens_unit)[exo200], 25*sens_unit])
    plt.ylim([(0.1*sens_unit,0.01*sens_unit)[exo200], 75*sens_unit])
    plt.yscale('log')
    plt.legend(loc="lower right",numpoints=1,prop={'size':legend_fontsize})
    fig.set_size_inches(7.5,5)

    if labels:
        if exo200:
            ax.annotate(r'1.9x10$^{25}$ yr', xy=(x200[0]+0.2,0.9*y200[0]), color='green', fontsize=labels_fontsize)#, xytext=(x200[0]+0.25,y200[0]*2.5), arrowprops=dict(arrowstyle='->'))
            ax.annotate('Nature 510, 229 (2014)', xy=(x200[0]+0.1,0.65*y200[0]), color='green', fontsize=9, style='italic')#, xytext=(x200[0]+0.25,y200[0]*2.5), arrowprops=dict(arrowstyle='->'))
        plt.plot([xfom],[yfom], 'rx', markersize=10, markeredgewidth=2)#, markerfacecolor='none' )
        ax.annotate('%.1fx10$^{%d}$ yr'%(yfom/(10**int(np.log10(yfom))),np.log10(yfom)), color='red',xy=(xfom-0.85,1.25*yfom), fontsize=labels_fontsize)#, xytext=(x200[0]+0.25,y200[0]*2.5), arrowprops=dict(arrowstyle='->'))

    plt.minorticks_on
    ax.grid(True, which='both')
    plt.savefig( outname+'.pdf', bbox_inches='tight' )
    plt.savefig( outname+'.png', bbox_inches='tight' )

    plt.show()



def MakePlot5p5( inpkl, bapkl, outname = 'out_sens_plot.pdf', discpkl = None, ba_discpkl = None):

    axis_label_fontsize = 15
    legend_fontsize = 14
    labels_fontsize = 14

    results = pickle.load(open(inpkl,'rb'))
    livetimes = results['livetimes']
    sens_unit = results['sens_unit']
    median_pow = results['ul_minos_med_pow']
    lolim_pow = results['ul_minos_lo_pow']
    hilim_pow = results['ul_minos_hi_pow']

    if discpkl and ba_discpkl:
        disc_results = pickle.load(open(discpkl,'rb'))
        disc_pow = disc_results['disc_minos_pow']
        ba_disc_results = pickle.load(open(ba_discpkl,'rb'))
        ba_disc_pow = ba_disc_results['disc_minos_pow']

    ba_results = pickle.load(open(bapkl,'rb'))
    ba_livetimes = ba_results['livetimes']
    ba_sens_unit = ba_results['sens_unit']
    ba_median_pow = ba_results['ul_minos_med_pow']
    ba_lolim_pow = ba_results['ul_minos_lo_pow']
    ba_hilim_pow = ba_results['ul_minos_hi_pow']

    fv_scale = 1.035
    time_scale = 0.0002
    max_time = 5.0
    xlist = np.arange(time_scale,max_time+time_scale,time_scale)
    ylist, ylolist, yhilist, dlist = [], [], [], []
    ba_list, ba_lolist, ba_hilist, ba_disc_list = [], [], [], []

    for x in xlist:
        ylist.append(median_pow[0]*(x**median_pow[1])*sens_unit*fv_scale)
        ylolist.append(lolim_pow[0]*(x**lolim_pow[1])*sens_unit*fv_scale)
        yhilist.append(hilim_pow[0]*(x**hilim_pow[1])*sens_unit*fv_scale)
        ba_list.append(ba_median_pow[0]*(x**ba_median_pow[1])*ba_sens_unit*fv_scale)
        ba_lolist.append(ba_lolim_pow[0]*(x**ba_lolim_pow[1])*ba_sens_unit*fv_scale)
        ba_hilist.append(ba_hilim_pow[0]*(x**ba_hilim_pow[1])*ba_sens_unit*fv_scale)         
        if discpkl and ba_discpkl:
            dlist.append(disc_pow[0]*(x**disc_pow[1])*sens_unit*fv_scale)
            ba_disc_list.append(ba_disc_pow[0]*(x**ba_disc_pow[1])*ba_sens_unit*fv_scale)
    
    ba_comb_list, ba_comb_disc_list = [], []
    for i,x in enumerate(xlist):
        ba_comb_list.append((ylist[-1]**2 + ba_list[i]**2)**0.5)
        if discpkl and ba_discpkl:
            ba_comb_disc_list.append((dlist[-1]**2 + ba_disc_list[i]**2)**0.5)

    y200 = [1.9e25/sens_unit]
    x200 = (np.log(y200[0]) - np.log(median_pow[0]))/median_pow[1]
    x200 = [np.exp(x200)]
    y200 = [y200[0]*sens_unit]
    print x200, y200

    ba_xlist, sens_cont_ylist, disc_cont_ylist = [], [], []
    for x in sorted(xlist):
        xx = x + max_time
        ba_xlist.append(xx)
        sens_cont_ylist.append(median_pow[0]*(xx**median_pow[1])*sens_unit*fv_scale)
        if discpkl and ba_discpkl:
            disc_cont_ylist.append(disc_pow[0]*(xx**disc_pow[1])*sens_unit*fv_scale)

    #fig = plt.figure()
    fig, ax = plt.subplots()
    plt.plot( ba_xlist, sens_cont_ylist, 'r-', linewidth=2, alpha=0.25)
    plt.plot( xlist, ylist, 'r-', markersize=4, markerfacecolor='r', linewidth=2, label="nEXO Sensitivity (90% C.L.)")
    #fill_between( xlist, ylolist, yhilist, color='c', edgecolor='None', alpha=0.5, label="68% Confidence Interval" ) 

    #plt.plot( ba_xlist, ba_list, 'r-', markersize=4, markerfacecolor='r', linewidth=2)
    plt.plot( ba_xlist, ba_comb_list, 'r-', markersize=4, markerfacecolor='r', linewidth=2)

    if discpkl and ba_discpkl:
        plt.plot( ba_xlist, disc_cont_ylist, 'b--', linewidth=2, alpha=0.25)
        plt.plot( xlist, dlist, 'b--', markersize=4, markerfacecolor='b', linewidth=2, label=r"nEXO Discovery 3$\sigma$, Prob. 50%")
        #plt.plot( ba_xlist, ba_disc_list, 'b--', markersize=4, markerfacecolor='b', linewidth=2)
        plt.plot( ba_xlist, ba_comb_disc_list, 'b--', markersize=4, markerfacecolor='b', linewidth=2)


    plt.plot( x200, y200, 'o', markersize=10, markeredgewidth=2, markerfacecolor='none', label="EXO-200 (Nature 510, 2014)", markeredgecolor='green')

    plt.plot( [max_time,max_time], [0.01*sens_unit, 50*sens_unit], 'k--', linewidth=1.5)

    plt.gca().set_xticks(range(0,livetimes[-1]+1))

    plt.xlabel("Livetime [yr]", fontsize=axis_label_fontsize)
    #plt.ylabel("T$_{1/2}$ [10$^{%d}$ yr]"%(np.log10(sens_unit)), fontsize=16)
    plt.ylabel(r"$^{136}$Xe  $0\nu\beta\beta$  T$_{1/2}$ [yr]", fontsize=axis_label_fontsize)

    plt.xlim([-0.2, 10.0])
    plt.ylim([0.01*sens_unit, 50*sens_unit])
    plt.yscale('log')
    plt.legend(loc="lower right",numpoints=1,prop={'size':legend_fontsize})
    fig.set_size_inches(7.5,5)

    ax.annotate('Start of Ba-tagging System', xy=(max_time,ba_comb_list[0]), color='k', fontsize=labels_fontsize, xytext=(0,ba_comb_list[-1]), arrowprops=dict(arrowstyle='->',color='r'))# xytext=(x200[0]+0.25,y200[0]*2.5), arrowprops=dict(arrowstyle='->'))
    ax.annotate('Start of Ba-tagging System', xy=(max_time,ba_comb_disc_list[0]), color='none', fontsize=labels_fontsize, xytext=(0,ba_comb_list[-1]), arrowprops=dict(arrowstyle='->',color='b'))# xytext=(x200[0]+0.25,y200[0]*2.5), arrowprops=dict(arrowstyle='->'))
    ax.annotate('1.9x10$^{25}$', xy=(x200[0]+0.2,0.9*y200[0]), color='green', fontsize=labels_fontsize)#, xytext=(x200[0]+0.25,y200[0]*2.5), arrowprops=dict(arrowstyle='->'))
    ax.annotate('Nature 510, 229 (2014)', xy=(x200[0]+0.1,0.65*y200[0]), color='green', fontsize=9, style='italic')#, xytext=(x200[0]+0.25,y200[0]*2.5), arrowprops=dict(arrowstyle='->'))

    plt.minorticks_on
    ax.grid(True, which='both')
    plt.savefig( outname, bbox_inches='tight' )

    plt.show()



def MakePlot3( senspkl, bapkl, imprpkl, outname = 'out_sens3_plot.pdf', labels = False, exo200 = True):

    axis_label_fontsize = 15
    legend_fontsize = 14
    labels_fontsize = 14

    results = pickle.load(open(senspkl,'rb'))
    livetimes = results['livetimes']
    sens_unit = results['sens_unit']
    median_pow = results['ul_minos_med_pow']

    ba_results = pickle.load(open(bapkl,'rb'))
    ba_sens_unit = ba_results['sens_unit']
    ba_median_pow = ba_results['ul_minos_med_pow']

    impr_results = pickle.load(open(imprpkl,'rb'))
    impr_sens_unit = impr_results['sens_unit']
    impr_median_pow = impr_results['ul_minos_med_pow']
        
    fv_scale = 1.035
    time_scale = 0.0002
    xlist = np.arange(time_scale,livetimes[-1]+time_scale,time_scale)
    ysens, yba, yimpr = [], [], []
    for x in xlist:
        ysens.append(median_pow[0]*(x**median_pow[1])*sens_unit*fv_scale)
        yba.append(ba_median_pow[0]*(x**ba_median_pow[1])*ba_sens_unit*fv_scale)
        yimpr.append(impr_median_pow[0]*(x**impr_median_pow[1])*impr_sens_unit*fv_scale)

    y200 = [1.9e25/sens_unit]
    x200 = (np.log(y200[0]) - np.log(median_pow[0]))/median_pow[1]
    x200 = [np.exp(x200)]
    y200 = [y200[0]*sens_unit]
    print x200, y200

    fig, ax = plt.subplots()

    plt.plot( xlist, ysens, 'r-', linewidth=2, label="Baseline Concept")
    plt.plot( xlist, yimpr, 'b--', linewidth=2, label="Considerable Improvements")
    line, = plt.plot( xlist, yba, 'g--', linewidth=2, label="Ba-tagging Scenario")
    dashes = [3,2]
    line.set_dashes(dashes)

    if exo200:
        plt.plot( x200, y200, 'o', markersize=10, markeredgewidth=2, markerfacecolor='none', label="EXO-200", markeredgecolor='green')

    plt.gca().set_xticks(range(0,livetimes[-1]+1))

    plt.xlabel("Livetime [yr]", fontsize=axis_label_fontsize)
    plt.ylabel(r"$^{136}$Xe  $0\nu\beta\beta$  T$_{1/2}$ [yr]", fontsize=axis_label_fontsize)

    plt.xlim([-0.2, livetimes[-1]])
    plt.ylim([(0.1*sens_unit,0.01*sens_unit)[exo200], 75*sens_unit])
    plt.yscale('log')
    plt.legend(loc="lower right",numpoints=1,prop={'size':legend_fontsize})
    fig.set_size_inches(7.5,5)

    plt.title('nEXO Sensitivity (90% C.L.)')

    if labels:
        if exo200:
            ax.annotate(r'1.9x10$^{25}$ yr', xy=(x200[0]+0.2,0.9*y200[0]), color='green', fontsize=labels_fontsize)#, xytext=(x200[0]+0.25,y200[0]*2.5), arrowprops=dict(arrowstyle='->'))
            ax.annotate('Nature 510, 229 (2014)', xy=(x200[0]+0.1,0.65*y200[0]), color='green', fontsize=9, style='italic')#, xytext=(x200[0]+0.25,y200[0]*2.5), arrowprops=dict(arrowstyle='->'))
        #plt.plot([xfom],[yfom], 'rx', markersize=10, markeredgewidth=2)#, markerfacecolor='none' )
        #ax.annotate('%.1fx10$^{%d}$ yr'%(yfom/(10**int(np.log10(yfom))),np.log10(yfom)), color='red',xy=(xfom-0.85,1.25*yfom), fontsize=labels_fontsize)#, xytext=(x200[0]+0.25,y200[0]*2.5), arrowprops=dict(arrowstyle='->'))

    plt.minorticks_on
    ax.grid(True, which='both')
    plt.savefig( outname+'.pdf', bbox_inches='tight' )
    plt.savefig( outname+'.png', bbox_inches='tight' )

    plt.show()

if __name__ == "__main__":

    
    # MakePlot( "batag_sens_time.pkl", "plot_nexo_sens_vs_time", discpkl = "batag_disc_time.pkl" , imprpkl = None, labels = True, exo200 = True)
    MakePlot3( "sens_time_dbv73.pkl","batag_sens_time.pkl","sens_optimistic_time.pkl",outname="nexo_sens_designs_v1", exo200=False)
