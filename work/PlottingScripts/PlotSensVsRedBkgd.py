
import ROOT

ROOT.gSystem.Load('../lib/libnEXOSensitivity.so')

resultsNameTemp = '../results/done/fits_hamamatsu_v68_2016-06-21_0nu_red%sx_fine_rdm_5.0_years_0.0_counts_*.root'

def GetBkgdCounts(filename,branch):

    chain = ROOT.TChain('tree')
    chain.Add(filename)
    chain.SetEstimate(chain.GetEntries()+1)
    chain.Draw(branch,'','goff')

    print 'using', chain.GetSelectedRows()
    
    return ROOT.TMath.Median(chain.GetSelectedRows(),chain.GetV1())

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

redList = ['0.25','0.5','0.75','1','2','4','6','8','12']

nameList = []
effList = {}
for red in reversed(redList):
    nameList.append(resultsNameTemp % (red))
    effList[nameList[-1]] = 0.82

#nameList.append('../results/done/fits_hamamatsu_v61_2016-02-24_Si_Cu_0nu_tpc_cryo_elec_fine_rdm_5.0_years_0.0_counts_*.root')
#effList[nameList[-1]] = 0.82
#nameList.append('../results/done/fits_hamamatsu_v61_0nu_eff_tpc_cryo_elec_4main_rdm_5.0_years_0.0_counts_*.root')
#effList[nameList[-1]] = 1
#nameList.append('../results/done/fits_hamamatsu_v61_0nu_eff_tpc_cryo_elec_4main_further_rdm_5.0_years_0.0_counts_*.root')
#effList[nameList[-1]] = 1

graph = ROOT.TGraph()
#graph.SetPoint(0,0,28)
graphul = ROOT.TGraph()
graph1t = ROOT.TGraph()
graph3t = ROOT.TGraph()
graphfv = ROOT.TGraph()

for name in nameList:
    print 'working on', name,
    
    n = graph.GetN()
    x = GetBkgdLimit(name,'bkg_fwhm_3t',0.5)#GetBkgdCounts(name)
    xul = GetBkgdLimit(name,'bkg_fwhm_3t',0.9)
    y = ROOT.nEXOUtils.GetSensHalfLife(name,5.0,3740,'signal',effList[name])
    #y1t = ROOT.nEXOUtils.EvalCountSensHalfLife(name,'bkg_fwhm_1t*5',5.0,1000)*effList[name]*0.762
    #y3t = ROOT.nEXOUtils.EvalCountSensHalfLife(name,'bkg_fwhm_3t*5',5.0,3000)*effList[name]*0.762
    #yfv = ROOT.nEXOUtils.EvalCountSensHalfLife(name,'bkg_fwhm_fv*5',5.0,3740)*effList[name]*0.762
    print x, y#, y1t, y3t, yfv
    graph.SetPoint(n,x,y/1e27)
    graphul.SetPoint(n,xul,y/1e27)
    #graph1t.SetPoint(n,x,y1t/1e27)
    #graph3t.SetPoint(n,x,y3t/1e27)
    #graphfv.SetPoint(n,x,yfv/1e27)

canvas = ROOT.TCanvas()
leg = ROOT.TLegend(0.1,0.7,0.48,0.9)
#graphul.SetMarkerStyle(24)
#graphul.SetMarkerColor(1)
#graphul.SetLineColor(1)
#graphul.GetYaxis().SetTitle('Sensitivity @ 90% CL (#times 10^{27} yrs)')
#graphul.GetXaxis().SetTitle('Background in FWHM-3t (cts/yr)')
#graphul.Draw('APL')
#leg.AddEntry(graphul,'90% UL','PL')
graph.SetMarkerStyle(20)
graph.SetMarkerColor(1)
graph.GetYaxis().SetTitle('Sensitivity @ 90% CL (#times 10^{27} yrs)')
graph.GetXaxis().SetTitle('Background in FWHM-3t (cts/yr)')
graph.Draw('APL')

pf = ROOT.TF1('power','[0]*pow(x,[1])',0,100)
pf.SetParameters(1,-0.5)
for i in range(10): graph.Fit(pf,'QW')
graph.Fit(pf,'W')

#leg.AddEntry(graph,'Median','PL')
leg.AddEntry(graph,'2D Fit Analysis','PL')
# graph1t.SetMarkerStyle(20)
# graph1t.SetMarkerColor(4)
# graph1t.SetLineColor(4)
# graph1t.Draw('PL')
# leg.AddEntry(graph1t,'Counting Expt. FWHM-1t','PL')
# graph3t.SetMarkerStyle(20)
# graph3t.SetMarkerColor(3)
# graph3t.SetLineColor(3)
# graph3t.Draw('PL')
# leg.AddEntry(graph3t,'Counting Expt. FWHM-3t','PL')
# graphfv.SetMarkerStyle(20)
# graphfv.SetMarkerColor(2)
# graphfv.SetLineColor(2)
# graphfv.Draw('PL')
#leg.AddEntry(graphfv,'Counting Expt. FWHM-FV','PL')

leg.Draw()
# leg = ROOT.TLegend(0.1,0.7,0.48,0.9)
# for r,red in enumerate(redList):
#     graph = graphs[group]
#     graph.GetYaxis().SetTitle('Sensitivity wrt 1% Resolution')
#     graph.GetXaxis().SetTitle('Resolution (%)')
#     graph.SetLineColor(g+1)
#     graph.SetMarkerColor(g+1)
#     graph.SetMarkerStyle(20)
#     if g == 0:
#         graph.Draw('APL')
#         graph.GetXaxis().SetRangeUser(0,50)
#     else:
#         graph.Draw('PL')
#     leg.AddEntry(graph,group,'PL')
# leg.Draw()

raw_input('end')
    
