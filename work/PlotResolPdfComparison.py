
import ROOT
import copy

resolList = [0.005,0.008,0.010,0.012,0.015,0.017,0.020]

canvas = ROOT.TCanvas()
leg = ROOT.TLegend(0.1,0.7,0.48,0.9)
hxs = {}

fin = ROOT.TFile.Open('../histos/individual_histos/nEXO_Histos_LXe_bb0n.root')
hin = fin.Get('h_StandoffVsEnergySS_Smear').Clone()
hx = hin.ProjectionX('%s_ref'%(hin.GetName()))
hx.Rebin(20)
hx.Scale(1/hx.Integral())
hx.SetLineColor(1)
hx.SetLineWidth(2)
hxs['ref'] = hx
canvas.cd()
hx.Draw('C')
leg.AddEntry(hx,'default','L')


for r,resol in enumerate(resolList):
    fin = ROOT.TFile.Open('../histos/individual_histos/nEXO_Histos_LXe_bb0n_resol%.3f.root'%(resol))
    hin = fin.Get('h_StandoffVsEnergySS_Smear').Clone()
    hx = hin.ProjectionX('%s_%.3f'%(hin.GetName(),resol))
    if resol < 0.016:
        hx.Rebin(20)
    hx.Scale(1/hx.Integral())
    hx.SetLineColor(r+2)
    hxs[resol] = hx
    canvas.cd()
    #draw = ('C','C same')[r > 0]
    
    hx.Draw('C same')#draw)
    leg.AddEntry(hx,'%.1f %%'%(resol*100.),'L')
    #fin.Close()


leg.Draw()
raw_input('end')
    
