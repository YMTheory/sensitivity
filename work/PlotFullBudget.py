import ROOT

resultsFileName = 'allfits_hamamatsu_0nu_eff_tpc_rdm_%0.1f_years_0.0_counts.root'

ROOT.gSystem.Load('../lib/libnEXOSensitivity.so')

nu = ROOT.nEXOUtil()

all_years = [0.5,1.0,2.5,5.0,10.,12.5,15.0,20.0] 
tgt_year = 5.0

gr = ROOT.TGraph()

gr3t = ROOT.TGraph()
igr3t = ROOT.TGraph()
gr1t = ROOT.TGraph()
igr1t = ROOT.TGraph()
grfv = ROOT.TGraph()
igrfv = ROOT.TGraph()

for y,years in enumerate(all_years):
    inFileName = '../results/full/%s' % (resultsFileName % years) 

    hl = nu.getHalfLife(inFileName,years,3740)

    gr.SetPoint(y,years,hl)

    hlb = nu.getHalfLife(inFileName,tgt_year,3740)
    bg3t = nu.getBkgdCounts(inFileName,2)*years
    bg1t = nu.getBkgdCounts(inFileName,3)*years
    bgFV = nu.getBkgdCounts(inFileName,1)*years

    gr3t.SetPoint(y,bg3t,hlb)
    igr3t.SetPoint(y,hlb,bg3t)
    gr1t.SetPoint(y,bg1t,hlb)
    igr1t.SetPoint(y,hlb,bg1t)
    grfv.SetPoint(y,bgFV,hlb)
    igrfv.SetPoint(y,hlb,bgFV)
   

spl3t = ROOT.TSpline3('spl3t',igr3t)
spy3t = spl3t.Eval(1e28)/tgt_year
spl1t = ROOT.TSpline3('spl1t',igr1t)
spy1t = spl1t.Eval(1e28)/tgt_year
splfv = ROOT.TSpline3('splfv',igrfv)
spyfv = splfv.Eval(1e28)/tgt_year

fp3t = ROOT.TF1('fp3t','[0]*pow(x,[1])',0.01,100.)
fp3t.SetParameters(1e28,-0.5)
for i in range(10):
    gr3t.Fit(fp3t,'RQW')
fx3t = fp3t.GetX(1e28)/tgt_year

fp1t = ROOT.TF1('fp1t','[0]*pow(x,[1])',0.01,100.)
fp1t.SetParameters(1e28,-0.5)
for i in range(10):
    gr1t.Fit(fp1t,'RQW')
fx1t = fp1t.GetX(1e28)/tgt_year

fpfv = ROOT.TF1('fpfv','[0]*pow(x,[1])',0.01,100.)
fpfv.SetParameters(1e28,-0.5)
for i in range(10):
    grfv.Fit(fpfv,'RQW')
fxfv = fpfv.GetX(1e28)/tgt_year

print '%0.2f - %0.2f | %0.3f - %0.3f | %0.2f - %0.2f ' % (fx3t, spy3t, fx1t, spy1t, fxfv, spyfv)

c = ROOT.TCanvas()
c.Divide(2,1)

c.cd(1)
gr.GetXaxis().SetTitle('Detector Livetime (years)')
gr.GetYaxis().SetTitle('0nu Halflife Sensitivity (years)')
gr.SetMarkerStyle(20)
gr.SetMarkerSize(0.7)
gr.SetMarkerColor(1)
gr.SetLineColor(1)
gr.SetLineWidth(1)
fp = ROOT.TF1('fp','[0]*pow(x,[1])',0.1,21.)
fp.SetParameters(1e28,0.5)
for i in range(10):
    gr.Fit(fp,'QRW')
gr.GetFunction('fp').SetLineColor(1)
gr.GetFunction('fp').SetLineWidth(2)
gr.GetFunction('fp').SetLineStyle(7)
    
gr.Draw('APC')

c.cd(2)
leg = ROOT.TLegend(0.1,0.7,0.48,0.9)

gr3t.GetXaxis().SetTitle('Mean Background in FWHM-xt / Year (counts)')
gr3t.GetYaxis().SetTitle('Achieved 0nu Halflife Sensitivity in 5 years (years)')
gr3t.SetMarkerStyle(20)
gr3t.SetMarkerSize(0.7)
gr3t.SetMarkerColor(4)
gr3t.SetLineColor(4)
gr3t.SetLineWidth(1)
gr3t.GetFunction('fp3t').SetLineColor(4)
gr3t.GetFunction('fp3t').SetLineWidth(2)
gr3t.GetFunction('fp3t').SetLineStyle(7)
    
gr3t.Draw('APC')
leg.AddEntry(gr3t,'3 tonne','PL')

gr1t.SetMarkerStyle(20)
gr1t.SetMarkerSize(0.7)
gr1t.SetMarkerColor(3)
gr1t.SetLineColor(3)
gr1t.SetLineWidth(1)
gr1t.GetFunction('fp1t').SetLineColor(3)
gr1t.GetFunction('fp1t').SetLineWidth(2)
gr1t.GetFunction('fp1t').SetLineStyle(7)
    
gr1t.Draw('PC')
leg.AddEntry(gr1t,'1 tonne','PL')

grfv.SetMarkerStyle(20)
grfv.SetMarkerSize(0.7)
grfv.SetMarkerColor(1)
grfv.SetLineColor(1)
grfv.SetLineWidth(1)
grfv.GetFunction('fpfv').SetLineColor(1)
grfv.GetFunction('fpfv').SetLineWidth(2)
grfv.GetFunction('fpfv').SetLineStyle(7)
    
grfv.Draw('PC')
leg.AddEntry(grfv,'FV = 3.74 tonne','PL')

leg.Draw()

raw_input('end')
