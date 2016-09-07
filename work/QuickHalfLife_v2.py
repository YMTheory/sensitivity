
import ROOT

resultsFileName = 'allfits_hamamatsu_0nu_eff_tpc_rdm_%s_%0.2f_%0.1f_years_0.0_counts.root'

ROOT.gSystem.Load('../lib/libnEXOSensitivity.so')

nu = ROOT.nEXOUtil()

bkgs = [0.25,0.5,1.0,2.0,4.0] #[[1.0,2.0,2.5,3.0,4.0]
group = 'InternalTh232' #'VesselTh232' #'LXeXe137' #'Far' #'LXeXe137' #'VesselU238' #'LXeRn222' #'InternalU238' 'VesselU238', Far (cryo Th232), Xe137
all_years = [5.0] #[0.5,1.0,2.5,5.0,10.,12.5,15.0,20.0] #[5.0] #[0.5,1.0,2.5,5.0,10.]
graph = ROOT.TGraph()
invgr = ROOT.TGraph()
gr1t = ROOT.TGraph()
grFV = ROOT.TGraph()
for years in all_years:
    for bkg in bkgs:
        inFileName = '../results/full/%s' % (resultsFileName % (group,bkg,years))
        hl = nu.getHalfLife(inFileName,years,3740)
        bg3t = nu.getBkgdCounts(inFileName,2)
        bg1t = nu.getBkgdCounts(inFileName,3)
        bgFV = nu.getBkgdCounts(inFileName,1)
        graph.SetPoint(graph.GetN(),bg3t,hl)
        invgr.SetPoint(invgr.GetN(),hl,bg3t)
        gr1t.SetPoint(gr1t.GetN(),bg3t,bg1t/bg3t)
        grFV.SetPoint(grFV.GetN(),bg3t,bgFV/bg3t)


c = ROOT.TCanvas()
c.Divide(2,1)

c.cd(1)
evmod = ROOT.TF1('evmod','[0]*pow(x,[1])',0.2,4.1)
evmod.SetParameters(1e27,0.5)
for i in range(10):
    graph.Fit(evmod,'RQ')
graph.Fit(evmod,'R')
graph.GetXaxis().SetTitle('Background in FWHM-3t (counts)')
graph.GetYaxis().SetTitle('0nu Halflife (years)')
graph.Draw('Ap*')

c.cd(2)
spline = ROOT.TSpline3('spline',invgr)
invgr.GetXaxis().SetTitle('0nu Halflife (years)')
invgr.GetYaxis().SetTitle('Background in FWHM-3t (counts)')
invgr.Draw('AP*C')
#spline.Draw('')

fitx = evmod.GetX(1e28)
sply = spline.Eval(1e28)

print 'For 1e28 yrs =>', fitx, sply 

#c.cd(3)
grFV.GetXaxis().SetRangeUser(0,5)
grFV.Fit('pol0')
grFV.SetMarkerColor(1)
#grFV.Draw('Ap*')
print 'FV', grFV.GetFunction('pol0').Eval(fitx)*fitx, grFV.GetFunction('pol0').Eval(sply)*sply
gr1t.Fit('pol0')
gr1t.SetMarkerColor(4)
#gr1t.Draw('p*')
print '1t', gr1t.GetFunction('pol0').Eval(fitx)*fitx, gr1t.GetFunction('pol0').Eval(sply)*sply


raw_input('end')
