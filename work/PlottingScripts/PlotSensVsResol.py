
import ROOT

ROOT.gSystem.Load('../lib/libnEXOSensitivity.so')

# nexo_batag_bkgd = {0.005: {'1-sigma': 0.00012156702423737897, '3-sigma': 0.0033944465261427914, 'fwhm': 0.0001692518237632612, 'klz': 6.297484902155599, '2-sigma': 0.0007526356585967342}, 0.025: {'1-sigma': 1.3823112294776365, '3-sigma': 40.42842561928248, 'fwhm': 1.990805720575736, 'klz': 21.588934716815682, '2-sigma': 8.782721038129239}, 0.015: {'1-sigma': 0.07094488432994694, '3-sigma': 2.174251918499273, 'fwhm': 0.10469339584233239, 'klz': 10.598819211820723, '2-sigma': 0.4668243877624718}, 0.01: {'1-sigma': 0.0067560696852524416, '3-sigma': 0.20346954650063775, 'fwhm': 0.009434754255835287, 'klz': 7.8150010856290635, '2-sigma': 0.04540726057558686}, 0.047: {'1-sigma': 46.97428731061518, '3-sigma': 1285.7508966185196, 'fwhm': 67.73201160365716, 'klz': 96.1696785929671, '2-sigma': 295.85924841600354}, 0.029: {'1-sigma': 3.227244202804286, '3-sigma': 91.78470542957535, 'fwhm': 4.545140399015509, 'klz': 28.837731118066934, '2-sigma': 20.473612865109317}, 0.02: {'1-sigma': 0.3891544728685403, '3-sigma': 11.404611063743879, 'fwhm': 0.5427439321938436, 'klz': 15.018502268045767, '2-sigma': 2.5010184130460402}}
# fc = ROOT.TFeldmanCousins()
# fc.SetMuMax(500)

# InternalTh232,LXeBb2n,LXeXe137,LXeRn222,
groups = ['LXeBb2n','InternalU238','LXeRn222','Far','InternalTh232','LXeXe137','All']

resultsNameTemp = '../results/done/fits_hamamatsu_0nu_eff_tpc_rdm_groupon_[GROUP]_resol%.3f_5.0_years_0.0_counts_*.root'

resultsNames = {}
graphs = {}

for group in groups:
    resultsNames[group] = resultsNameTemp.replace('[GROUP]',group)
    graphs[group] = ROOT.TGraph()


resolList = [0.005,0.010,0.015] #,0.025,0.029,0.047] #[0.003,0.005,0.008,0.010,0.012,0.015,0.017,0.020,0.025,0.030,0.040]

for group in groups:
    yref = ROOT.nEXOUtils.GetSensHalfLife(resultsNames[group]%(0.010),5.0,3740)
    for resol in resolList:
        print 'working on', group, resol

        n = graphs[group].GetN()
        x = resol
        x *= 100

        y = ROOT.nEXOUtils.GetSensHalfLife(resultsNames[group]%(resol),5.0,3740)
        y /= yref #1e27
        graphs[group].SetPoint(n,x,y)

canvas = ROOT.TCanvas()
leg = ROOT.TLegend(0.1,0.7,0.48,0.9)
for g,group in enumerate(groups):
    graph = graphs[group]
    graph.GetYaxis().SetTitle('Sensitivity wrt 1% Resolution')
    graph.GetXaxis().SetTitle('Resolution (%)')
    graph.SetLineColor(g+1)
    graph.SetMarkerColor(g+1)
    graph.SetMarkerStyle(20)
    if g == 0:
        graph.Draw('APL')
        graph.GetXaxis().SetRangeUser(0,50)
    else:
        graph.Draw('PL')
    leg.AddEntry(graph,group,'PL')
leg.Draw()

raw_input('end')
    
