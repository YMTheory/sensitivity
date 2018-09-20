import ROOT

resultsFileName = 'allfits_hamamatsu_0nu_eff_tpc_rdm_%s_%0.2f_%0.1f_years_0.0_counts.root'

ROOT.gSystem.Load('../lib/libnEXOSensitivity.so')

nu = ROOT.nEXOUtil()

bkgs = [0.25,0.5,1.0,2.0,4.0] #[[1.0,2.0,2.5,3.0,4.0]
groups = ['LXeXe137','LXeRn222','InternalU238','VesselU238','InternalTh232','VesselTh232','Far'] 
years = 5.0

gr3ts = {}
igr3ts = {}
gr1ts = {}
igr1ts = {}
grfvs = {}
igrfvs = {}

for group in groups:
    print 'Working on group', group

    gr3ts[group] = ROOT.TGraph()
    igr3ts[group] = ROOT.TGraph()
    gr1ts[group] = ROOT.TGraph()
    igr1ts[group] = ROOT.TGraph()
    grfvs[group] = ROOT.TGraph()
    igrfvs[group] = ROOT.TGraph()
    
    for b,bkg in enumerate(bkgs):
        inFileName = '../results/full/%s' % (resultsFileName % (group,bkg,years))

        hl = nu.getHalfLife(inFileName,years,3740)
        bg3t = nu.getBkgdCounts(inFileName,2)
        bg1t = nu.getBkgdCounts(inFileName,3)
        bgFV = nu.getBkgdCounts(inFileName,1)
        
        gr3ts[group].SetPoint(b,bg3t,hl)
        igr3ts[group].SetPoint(b,hl,bg3t)
        gr1ts[group].SetPoint(b,bg1t,hl)
        igr1ts[group].SetPoint(b,hl,bg1t)
        grfvs[group].SetPoint(b,bgFV,hl)
        igrfvs[group].SetPoint(b,hl,bgFV)

fp3ts = {}
fx3ts = {}
spl3ts = {}
spy3ts = {}

fp1ts = {}
fx1ts = {}
spl1ts = {}
spy1ts = {}

fpfvs = {}
fxfvs = {}
splfvs = {}
spyfvs = {}


print '***********************************************'
print 'GROUP | 3t | 1t | FV'
for group in groups:

    spl3ts[group] = ROOT.TSpline3('spl3t_%s'%(group),igr3ts[group])
    spy3ts[group] = spl3ts[group].Eval(1e28)

    fp3ts[group] = ROOT.TF1('fp3t_%s'%(group),'[0]*pow(x,[1])',0.2,4.1)
    fp3ts[group].SetParameters(1e28,-0.5)
    for i in range(10):
        gr3ts[group].Fit(fp3ts[group],'RQW')
    fx3ts[group] = fp3ts[group].GetX(1e28)

    spl1ts[group] = ROOT.TSpline3('spl1t_%s'%(group),igr1ts[group])
    spy1ts[group] = spl1ts[group].Eval(1e28)

    fp1ts[group] = ROOT.TF1('fp1t_%s'%(group),'[0]*pow(x,[1])',0.01,1)
    fp1ts[group].SetParameters(1e28,-0.5)
    for i in range(10):
        gr1ts[group].Fit(fp1ts[group],'RQW')
    fx1ts[group] = fp1ts[group].GetX(1e28)

    splfvs[group] = ROOT.TSpline3('splfv_%s'%(group),igrfvs[group])
    spyfvs[group] = splfvs[group].Eval(1e28)

    fpfvs[group] = ROOT.TF1('fpfv_%s'%(group),'[0]*pow(x,[1])',0.1,10)
    fpfvs[group].SetParameters(1e28,-0.5)
    for i in range(10):
        grfvs[group].Fit(fpfvs[group],'RQW')
    fxfvs[group] = fpfvs[group].GetX(1e28)


    print '%s | %0.2f - %0.2f | %0.3f - %0.3f | %0.2f - %0.2f ' % (group, fx3ts[group], spy3ts[group], fx1ts[group], spy1ts[group], fxfvs[group], spyfvs[group])
print '***********************************************'

c = ROOT.TCanvas()
c.Divide(2,1)

c.cd(1)
leg3t = ROOT.TLegend(0.1,0.7,0.48,0.9)
for g,group in enumerate(groups):
    graph = gr3ts[group]
    graph.GetXaxis().SetTitle('Mean Background in FWHM-3t / Year (counts)')
    graph.GetYaxis().SetTitle('Achieved 0nu Halflife Sensitivity in 5 Years (years)')
    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(0.7)
    graph.SetMarkerColor(g+1)
    graph.SetLineColor(g+1)
    graph.SetLineWidth(1)
    graph.GetFunction('fp3t_%s'%(group)).SetLineColor(g+1)
    graph.GetFunction('fp3t_%s'%(group)).SetLineWidth(1)
    graph.GetFunction('fp3t_%s'%(group)).SetLineStyle(7)
    
    graph.GetYaxis().SetRangeUser(1e27,2e28)
    graph.Draw(('PC','APC')[g==0])
    leg3t.AddEntry(graph,group,'PL')
leg3t.Draw()

c.cd(2)
leg1t = ROOT.TLegend(0.1,0.7,0.48,0.9)
for g,group in enumerate(groups):
    graph = gr1ts[group]
    graph.GetXaxis().SetTitle('Mean Background in FWHM-1t / Year (counts)')
    graph.GetYaxis().SetTitle('Achieved 0nu Halflife Sensitivity in 5 Years (years)')
    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(0.7)
    graph.SetMarkerColor(g+1)
    graph.SetLineColor(g+1)
    graph.SetLineWidth(1)
    graph.GetFunction('fp1t_%s'%(group)).SetLineColor(g+1)
    graph.GetFunction('fp1t_%s'%(group)).SetLineWidth(1)
    graph.GetFunction('fp1t_%s'%(group)).SetLineStyle(7)
    
    graph.GetXaxis().SetRangeUser(0.01,0.5)
    graph.GetYaxis().SetRangeUser(1e27,2e28)
    graph.Draw(('PC','APC')[g==0])
    leg1t.AddEntry(graph,group,'PL')
leg1t.Draw()

raw_input('end')
