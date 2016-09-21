import ROOT

import cPickle as pickle
import glob, os
import matplotlib.pyplot as plt
import numpy as np

nEXOlib = '/data/data033/exo/software/nEXO_Sensitivity/quick/v5/lib/libnEXOSensitivity.so'
ROOT.gSystem.Load(nEXOlib)
random = ROOT.TRandom3(1)

def ShowPickleResults( inpkl ):

    results = pickle.load(open(inpkl,'rb'))
    print results


def EvalCounts(hitEfficiency, activity, time, halflife):
    if time*1./halflife > 0.01:  
        lhl = ROOT.TMath.Log2(halflife)
        time = lhl * (1 - np.exp(-time/lhl))
    return time * hitEfficiency * activity * 31556736 # seconds per year conversion


def EvalSignalCounts( halflife, livetime, eff = 0.82, mass = 3740, enrich = 0.9, avog = 6.022e23, molar = 0.136 ):
    
    atoms = mass * enrich * avog / molar    
    return atoms * livetime * eff * np.log(2.) / halflife

def EvalSignalActivity( halflife, enrich = 0.9, avog = 6.022e23, molar = 0.136e-3 ):

    return np.log(2.) * enrich * avog / molar / halflife /  31556736 # seconds per year conversion

def ReadTable( infile, outpkl, livetime = 5.0, bb0n_halflife = 6.2e27 ):

    groups = {}
    groups['bb0n'] = ['LXe_bb0n']
    groups['bb2n'] = ['LXe_bb2n']   
    groups['other'] = ['']   
    
    prebinx, prebiny, minbin = 20,10,1000 
    nbinsx, xmin, xmax = 140, 700, 3500
    nbinsy, ymin, ymax = 65, 0, 650
    histos = {}
    for group in groups:
        histos[group] = {}
    histos['other'] = {}

    specActivities = {}

    chain = ROOT.TChain("ExcelTableValues")
    table = ROOT.ExcelTableValues()
    chain.Add(infile)
    chain.SetBranchAddress("table",table)
    for i in range(chain.GetEntries()):
        chain.GetEntry(i)

        #if table.fIsotope == 'K-40':
        #    continue

        table.Print()
      
        for s in range(table.fSuffixes.size()):
            
            # make histogram
            groupFound = False
            for group in groups:
                if table.fPdf in groups[group]:
                    groupFound = True
                    if not s in histos[group]:
                        histos[group][s] = ROOT.TH2D('%s_%s'%(group,table.fSuffixes[s]),'',nbinsx,xmin,xmax,nbinsy,ymin,ymax)
                    curgroup = group
                    break
            if not groupFound:
                groups['other'].append(table.fPdf)
                if not s in histos['other']:
                    histos['other'][s] = ROOT.TH2D('%s_%s'%('other',table.fSuffixes[s]),'',nbinsx,xmin,xmax,nbinsy,ymin,ymax)
                curgroup = 'other'

            # get activity
            specActivity = table.fSpecActivCV
            activID = "%s_%s" % (table.fActivID,table.fIsotope)
            if activID in specActivities or table.fPdf == 'LXe_bb0n':
                if table.fPdf == 'LXe_bb0n':
                    print '#############', EvalSignalActivity(bb0n_halflife), 1./ ((bb0n_halflife / 2.165e21) * 0.0405362) * 60./46.
                    specActivities[activID] = EvalSignalActivity(bb0n_halflife)
                    #specActivities[activID] = 1./ ((bb0n_halflife / 2.165e21) * 0.0405362) * 60./46.  # use bb2n half-life to convert bb0n into specific activity
                specActivity = specActivities[activID]
            else:
                specActivity = []
                for iDraw in range(1000):
                    specActivityI = -1
                    nDraws = 1000
                    while specActivityI < 0 and nDraws > 0:
                        specActivityI = random.Gaus(table.fSpecActivCV,table.fSpecActivError)
                        nDraws -= 1
                    if specActivityI < 0 and nDraws == 0:
                        specActivityI = 1e-16
                    specActivity.append(specActivityI)
                specActivity = ROOT.TMath.Median(len(specActivity),np.asarray(specActivity))
                specActivities[activID] = specActivity
            # convert to counts (normalization)
            if table.fSpecActivError > 0:
                factor = table.fActivError*1./table.fSpecActivError    
            else:
                factor = table.fActivCV*1./table.fSpecActivCV
            activity = factor * specActivity;
            if table.fPdf == 'LXe_bb0n':
                counts = EvalCounts(table.fHitEffK[s]/table.fHitEffN[s],activity,livetime,bb0n_halflife)
            else:
                counts = EvalCounts(table.fHitEffK[s]/table.fHitEffN[s],activity,livetime,table.fHalflife)
            # get histogram
            hfile = ROOT.TFile.Open(table.fFileName,'read')
            h2din = hfile.Get('h_StandoffVsEnergy%s_Smear'%(table.fSuffixes[s]))
            if h2din.GetNbinsX() > minbin:
                h2din.Rebin2D(prebinx,prebiny)
            h2din.Scale(counts*1./h2din.Integral())
            h2dout = ROOT.TH2D('temp_%s_%s'%(table.fPdf,table.fSuffixes[s]),'',nbinsx,xmin,xmax,nbinsy,ymin,ymax)
            for bx in range(1,h2din.GetNbinsX()+1):
                centerx = h2din.GetXaxis().GetBinCenter(bx)
                for by in range(1,h2din.GetNbinsY()+1):
                    #print bx, by
                    centery = h2din.GetYaxis().GetBinCenter(by)
                    content = h2din.GetBinContent(bx,by)
                    h2dout.Fill(centerx,centery,content)
            #h2dout.Scale(counts*1./h2dout.Integral())
            histos[curgroup][s].Add(h2dout)

            if table.fPdf == 'LXe_bb0n':
                print counts, h2din.Integral(), h2dout.Integral()
            if table.fPdf == 'LXe_bb2n':
                print counts, h2din.Integral(), h2dout.Integral()

            hfile.Close()

    pickle.dump(histos, open( outpkl, 'wb'))

def ProjectVolume( name, mult, vol, histos ):

    sd = {'fv':0, '3t': 9, '2t': 16, '1t': 26}
    m = {'ss':0, 'ms':1}

    return histos[name][m[mult]].ProjectionX('%s_ene_%s'%(name,vol),sd[vol])
    

def ConvertHistogram( inhist ):
    
    cs, ws, nb = [], [], inhist.GetNbinsX()
    for i in range(1,nb+1):
        cs.append(inhist.GetBinCenter(i))
        ws.append(inhist.GetBinContent(i))
    return nb, cs, ws

def PlotVolume( vol, mult, ax, pdfs, nbins, hrange ):

    for component in pdfs:#reversed(sorted(pdfs.keys())):
        pdf = pdfs[component]
        pdf_vol = pdf[mult][vol]
        ax.hist(pdf_vol['cs'], weights=pdf_vol['ws'], edgecolor=pdf['color'], color=pdf['color'], bins=nbins, histtype='stepfilled', alpha=pdf['alpha'] , hatch=pdf['hatch'], range=hrange, label=pdf['legend'])


def MakePlot( inpkl , outname = 'energy_vols.pdf' ):

    axis_label_fontsize = 14
    legend_fontsize = 14
    labels_fontsize = 14

    histos = pickle.load(open(inpkl,'rb'))
    hrange = [800,3500]
    hnb = (hrange[-1]-hrange[0])/20

    vols = ['fv','3t','2t','1t']
    components = ['bb2n','bb0n','other']
    colors = ['lightgreen','red','darkgreen']
    alphas = [0.75,0.5,0.5]
    hatchs = ['','','']
    legends = [r'$2\nu\beta\beta$',r'$0\nu\beta\beta$','Other Bgd.']
    mults = ['ss','ms']

    pdfs = {}
    for component, color, alpha, hatch, legend in zip(components,colors, alphas, hatchs, legends):
        pdfs[component] = {}
        pdfs[component]['color'] = color
        pdfs[component]['alpha'] = alpha
        pdfs[component]['hatch'] = hatch
        pdfs[component]['legend'] = legend
        for mult in mults:
            pdfs[component][mult] = {}
            for vol in vols:
                pdfs[component][mult][vol] = {}
                pdfs[component][mult][vol]['nb'], pdfs[component][mult][vol]['cs'], pdfs[component][mult][vol]['ws'] = ConvertHistogram(ProjectVolume(component,mult,vol,histos))
       
    fig, ((ax_ss_fv, ax_ss_3t, ax_ss_2t, ax_ss_1t), (ax_ms_fv, ax_ms_3t, ax_ms_2t, ax_ms_1t)) = plt.subplots(2, 4, sharex=True, sharey=True)
    
    #fig.xlabel("Energy [keV]", fontsize=axis_label_fontsize)
    #fig.ylabel("Events / 20 keV [cts]", fontsize=axis_label_fontsize)

    PlotVolume('fv', 'ss', ax_ss_fv, pdfs, hnb, hrange)
    PlotVolume('3t', 'ss', ax_ss_3t, pdfs, hnb, hrange)
    PlotVolume('2t', 'ss', ax_ss_2t, pdfs, hnb, hrange)
    PlotVolume('1t', 'ss', ax_ss_1t, pdfs, hnb, hrange)
    PlotVolume('fv', 'ms', ax_ms_fv, pdfs, hnb, hrange)
    PlotVolume('3t', 'ms', ax_ms_3t, pdfs, hnb, hrange)
    PlotVolume('2t', 'ms', ax_ms_2t, pdfs, hnb, hrange)
    PlotVolume('1t', 'ms', ax_ms_1t, pdfs, hnb, hrange)

    fig.subplots_adjust(wspace = 0, hspace = 0)
    
    plt.xlim(hrange)
    plt.ylim([5e-2,7e6])
    plt.yscale('log')

    ax_ss_fv.set_ylabel("SS Events [cts]", fontsize=axis_label_fontsize)
    ax_ms_fv.set_ylabel("MS Events [cts]", fontsize=axis_label_fontsize)
    ax_ms_fv.set_xlabel("Energy [MeV]", fontsize=axis_label_fontsize)
    ax_ms_3t.set_xlabel("Energy [MeV]", fontsize=axis_label_fontsize)
    ax_ms_2t.set_xlabel("Energy [MeV]", fontsize=axis_label_fontsize)
    ax_ms_1t.set_xlabel("Energy [MeV]", fontsize=axis_label_fontsize)
    #plt.figtext(0.45,0.02,'Energy [MeV]',fontdict={'fontsize':axis_label_fontsize})

    #labels = [item.get_text() for item in ax_ms_fv.get_xticklabels()]
    #for i,label in enumerate(labels):
    #    print i, labels[i], label
    ax_ms_fv.set_xticklabels(['',1.0,'',2.0,'',3.0,''])
    #plt.gca().set_xticks([1500,2500])

    ax_ss_fv.set_title("Fid. Vol.", fontsize=axis_label_fontsize)
    ax_ss_3t.set_title("Inner 3 t", fontsize=axis_label_fontsize)
    ax_ss_2t.set_title("Inner 2 t", fontsize=axis_label_fontsize)
    ax_ss_1t.set_title("Inner 1 t", fontsize=axis_label_fontsize)

    leg = plt.legend(loc="upper left",bbox_to_anchor = (1,0.5),prop={'size':legend_fontsize})

    #plt.minorticks_on #off 
    #plt.grid(False, which='both')

    fig.set_size_inches(8,4.5)
    plt.savefig( outname, bbox_extra_artists=(leg,), bbox_inches='tight' )

    plt.show()

    #hh, be = np.histogram( cs, bins=nb-15, weights=ws, range=[1000,3500] )
    #plt.step(be[:-1], hh, where='post', color='k', linewidth=1.5)
