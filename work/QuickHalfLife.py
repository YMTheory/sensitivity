
import ROOT

resultsFileName = 'fits_hamamatsu_v68_2016-06-21_0nu_ssfrac20_fine_rdm_%0.1f_years_0.0_counts_*.root' # 'fits_hamamatsu_0nu_eff_tpc_klzfix_%0.1f_years_0.0_counts_*.root'

ROOT.gSystem.Load('../lib/libnEXOSensitivity.so')

#nu = ROOT.nEXOUtil()

all_years = [5.0] #[0.5,1.0,2.5,5.0,10.] #[5.0] #[0.5,1.0,2.5,5.0,10.] #[0.5,1.0,2.5,5.0,10.] #[1.4,2.9,4.2,8.7] #[0.5,1.0,2.5,5.0,10.,12.5,15.0,20.0] #[5.0] #[0.5,1.0,2.5,5.0,10.]
graph = ROOT.TGraph()
graph2 = ROOT.TGraph()
invgr = ROOT.TGraph()
gr1t = ROOT.TGraph()
grFV = ROOT.TGraph()
for years in all_years:
    inFileName = '../results/done/%s' % (resultsFileName % years) 

    hl = ROOT.nEXOUtils.GetSensHalfLife(inFileName,years,3740,"signal",0.82)
    #hl = ROOT.nEXOUtils.GetSensHalfLife(inFileName,1.45,344.5)#years,3740)

    #hl = ROOT.nEXOUtils.GetHalfLife(inFileName,years,3740)
    #bkg = ROOT.nEXOUtils.GetBackgroundCounts(inFileName,"bkg_fwhm_fv")
    #hlb = ROOT.nEXOUtils.GetHalfLife(inFileName,5,3740)
    #bg1t = ROOT.nEXOUtils.GetBackgroundCounts(inFileName,"bkg_fwhm_3t")
    #bgFV = ROOT.nEXOUtils.GetBackgroundCounts(inFileName,"bkg_fwhm_3t")
    print years, hl#, bkg
    #graph.SetPoint(graph.GetN(),years,hl)
    #graph2.SetPoint(graph2.GetN(),bkg*years,hlb)
    #invgr.SetPoint(invgr.GetN(),hlb,bkg*years)
    #gr1t.SetPoint(gr1t.GetN(),bkg,bg1t/bkg)
    #grFV.SetPoint(grFV.GetN(),bkg,bgFV/bkg)

raw_input('continue')

c = ROOT.TCanvas()
c.Divide(2,1)

c.cd(1)
evmod = ROOT.TF1('evmod','[0]*pow(x,[1])',0,21)
evmod.SetParameters(1e27,0.5)
for i in range(10):
    graph.Fit(evmod,'RQ')
graph.Fit(evmod,'R')
graph.GetXaxis().SetTitle('Detector Livetime (years)')
graph.GetYaxis().SetTitle('0nu Halflife (years)')
graph.Draw('Ap*')

c.cd(2)
evmod2 = ROOT.TF1('evmod2','[0]*pow(x,[1])',1e-5,100)
evmod2.SetParameters(1e27,-0.5)
for i in range(10):
    graph2.Fit(evmod2,'RQ')
graph2.Fit(evmod2,'R')
graph2.GetXaxis().SetTitle('Background in FWHM-3t (counts)')
graph2.GetYaxis().SetTitle('0nu Halflife (years)')
graph2.Draw('Ap*C')

spline = ROOT.TSpline3('spline',invgr)

fitx = evmod2.GetX(1e28)/5
sply = spline.Eval(1e28)/5

print 'For 1e28 yrs =>', fitx, sply


grFV.GetXaxis().SetRangeUser(0,5)
grFV.Fit('pol0','Q')
grFV.SetMarkerColor(1)
#grFV.Draw('Ap*')
print 'FV', grFV.GetFunction('pol0').Eval(fitx)*fitx, grFV.GetFunction('pol0').Eval(sply)*sply
gr1t.Fit('pol0','Q')
gr1t.SetMarkerColor(4)
#gr1t.Draw('p*')
print 'FV', gr1t.GetFunction('pol0').Eval(fitx)*fitx, gr1t.GetFunction('pol0').Eval(sply)*sply

raw_input('end')
