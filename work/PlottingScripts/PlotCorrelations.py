
import ROOT
import copy

ROOT.gSystem.Load('../lib/libnEXOSensitivity.so')

#filename = '../v4/results/done/fits_hamamatsu_0nu_eff_tpc_rdm_5.0_years_0.0_counts_*.root'
filename = '../../v3/results/done/disc_hamamatsu_0nu_eff_tpc_rdm_5.0_years_10.0_counts_*.root'

chain = ROOT.TChain('tree')
chain.Add(filename)

#nfr = ROOT.nEXOFitResult()
#chain.SetBranchAddress('nEXOFitResult',nfr)
nfr = chain

reference = 'LXeBb0n'
comparisons = ['Far','FullTpcK40','FullTpcCo60','InternalTh232','InternalU238','LXeBb2n','LXeRn222','LXeXe137','VesselTh232','VesselU238']
legs = {'Far':' External Components','FullTpcK40':'TPC ^{40}K Backgrounds','FullTpcCo60':'TPC ^{60}Co Backgrounds','InternalTh232':'TPC ^{232}Th Backgrounds','InternalU238':'TPC ^{238}U Backgrounds','LXeBb2n':'#beta#beta2#nu','LXeRn222':'LXe Backgrounds','LXeXe137':'LXe Backgrounds','VesselTh232':'','VesselU238':''}

histos = {}
for compared in comparisons:
    histos[compared] = ROOT.TH1D('h_%s'%(compared),'',80,-1,1)

n = chain.GetEntries()
for i in range(n):
    chain.GetEntry(i)
    result = nfr.fitres_sig
    pars = result.floatParsFinal()

    if nfr.covQual_sig == 3 and nfr.stat_sig == 0:
        for compared in comparisons:
            histos[compared].Fill(result.correlation('num_%s'%(reference),'num_%s'%(compared)))

canvas = ROOT.TCanvas()
leg = ROOT.TLegend(0.1,0.7,0.48,0.9)
leg.AddEntry(0,'Mean Correlations for T^{0#nu#beta#beta}_{1/2} = 6#times10^{27} yrs','')
for c,compared in enumerate(comparisons):
    histo = histos[compared]
    
    histo.GetXaxis().SetTitle('Correlation between Component and 0#nu#beta#beta')
    histo.GetYaxis().SetTitle('Number of Toy Fits')

    histo.SetLineColor(c+1)

    draw_option = ''
    if c > 0:
        draw_option += 'same'

    canvas.cd()
    if abs(histo.GetMean()) > 0.01:
        histo.Draw(draw_option)
        leg.AddEntry(histo,'%s: %.2f'%(legs[compared],histo.GetMean()),'L')
    print compared, histo.GetMean()

canvas.cd()
leg.Draw()

raw_input('end')

        
