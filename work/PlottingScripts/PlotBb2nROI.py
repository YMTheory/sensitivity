
import ROOT

nexo = {0.005: {'1-sigma': 0.00012156702423737897, '3-sigma': 0.0033944465261427914, 'fwhm': 0.0001692518237632612, 'klz': 6.297484902155599, '2-sigma': 0.0007526356585967342}, 0.025: {'1-sigma': 1.3823112294776365, '3-sigma': 40.42842561928248, 'fwhm': 1.990805720575736, 'klz': 21.588934716815682, '2-sigma': 8.782721038129239}, 0.015: {'1-sigma': 0.07094488432994694, '3-sigma': 2.174251918499273, 'fwhm': 0.10469339584233239, 'klz': 10.598819211820723, '2-sigma': 0.4668243877624718}, 0.01: {'1-sigma': 0.0067560696852524416, '3-sigma': 0.20346954650063775, 'fwhm': 0.009434754255835287, 'klz': 7.8150010856290635, '2-sigma': 0.04540726057558686}, 0.047: {'1-sigma': 46.97428731061518, '3-sigma': 1285.7508966185196, 'fwhm': 67.73201160365716, 'klz': 96.1696785929671, '2-sigma': 295.85924841600354}, 0.029: {'1-sigma': 3.227244202804286, '3-sigma': 91.78470542957535, 'fwhm': 4.545140399015509, 'klz': 28.837731118066934, '2-sigma': 20.473612865109317}, 0.02: {'1-sigma': 0.3891544728685403, '3-sigma': 11.404611063743879, 'fwhm': 0.5427439321938436, 'klz': 15.018502268045767, '2-sigma': 2.5010184130460402}}
klzu = {0.005: {'1-sigma': 0.00011974399735947827, '3-sigma': 0.0033594812077865424, 'fwhm': 0.00016679324713209098, 'klz': 6.277691980182578, '2-sigma': 0.0007433006020418631}, 0.025: {'1-sigma': 1.2863696686399635, '3-sigma': 38.57113451546681, 'fwhm': 1.8575403469003504, 'klz': 20.505626336916766, '2-sigma': 8.286820582914515}, 0.015: {'1-sigma': 0.06785368563942029, '3-sigma': 2.1106198985184577, 'fwhm': 0.10029860512258892, 'klz': 10.35305068180718, '2-sigma': 0.4501921377859617}, 0.01: {'1-sigma': 0.006556537740834756, '3-sigma': 0.19938936392624984, 'fwhm': 0.009164839125332946, 'klz': 7.724708317088643, '2-sigma': 0.04431071400827591}, 0.047: {'1-sigma': 41.364399360027164, '3-sigma': 1186.7464309277852, 'fwhm': 59.93791936291382, 'klz': 85.51202678386471, '2-sigma': 267.35458876205666}, 0.029: {'1-sigma': 2.972276038664859, '3-sigma': 86.97410939495074, 'fwhm': 4.198150982731022, 'klz': 27.066324976875123, '2-sigma': 19.1597393540942}, 0.02: {'1-sigma': 0.36708573205396533, '3-sigma': 10.973135812079157, 'fwhm': 0.5129452611581655, 'klz': 14.47264017687168, '2-sigma': 2.3852696189496783}}
results = {}
results['nexo'] = nexo
results['klzu'] = klzu

graphs = {}
rois = ['fwhm','2-sigma']
funcs = ['nexo']#,'klzu']

legs = {}
legs['nexo'] = 'EXO-like'
legs['klzu'] = 'KLZ-like'
legs['2-sigma'] = '4 #sigma'
legs['fwhm'] = 'FWHM'

for func in funcs:
    graphs[func] = {}
    for roi in rois:
        graphs[func][roi] = ROOT.TGraph()

for resq in sorted(results['nexo']):
    x = resq*100.

    for func in funcs:
        for roi in rois:
            y = results[func][resq][roi]
            n = graphs[func][roi].GetN()
            graphs[func][roi].SetPoint(n,x,y)

canvas = ROOT.TCanvas()
leg = ROOT.TLegend(0.1,0.7,0.48,0.9)

for f, func in enumerate(graphs):
    for r, roi in enumerate(rois):
        graph = graphs[func][roi]
        graph.SetLineColor(f+1)
        graph.SetLineStyle(r+1)
        graph.SetLineWidth(2)
        graph.SetMarkerColor(f+1)
        graph.SetMarkerStyle(20+r)

        canvas.cd()
        draw_option = 'PL'
        if f+r == 0:
            draw_option += 'A'
            graph.GetXaxis().SetTitle('Resolution at the Q-value (%)')
            graph.GetXaxis().SetTitleSize(0.05)
            graph.GetXaxis().SetTitleOffset(0.85)
            graph.GetYaxis().SetTitle('#beta#beta2#nu Background (evt / t.yr)')
            graph.GetYaxis().SetTitleSize(0.05)
            graph.GetYaxis().SetTitleOffset(0.85)
        graph.Draw(draw_option)
        graph.GetYaxis().SetRangeUser(1e-4,1e3)

        #leg.AddEntry(graph,'ROI: %s, resolution: %s'%(legs[roi],legs[func]),'PL')
        leg.AddEntry(graph,'%s'%(legs[roi]),'PL')

leg.Draw()
raw_input('end')
