
import ROOT

nexo_batag_bkgd = {0.005: {'1-sigma': 0.00012156702423737897, '3-sigma': 0.0033944465261427914, 'fwhm': 0.0001692518237632612, 'klz': 6.297484902155599, '2-sigma': 0.0007526356585967342}, 0.025: {'1-sigma': 1.3823112294776365, '3-sigma': 40.42842561928248, 'fwhm': 1.990805720575736, 'klz': 21.588934716815682, '2-sigma': 8.782721038129239}, 0.015: {'1-sigma': 0.07094488432994694, '3-sigma': 2.174251918499273, 'fwhm': 0.10469339584233239, 'klz': 10.598819211820723, '2-sigma': 0.4668243877624718}, 0.01: {'1-sigma': 0.0067560696852524416, '3-sigma': 0.20346954650063775, 'fwhm': 0.009434754255835287, 'klz': 7.8150010856290635, '2-sigma': 0.04540726057558686}, 0.047: {'1-sigma': 46.97428731061518, '3-sigma': 1285.7508966185196, 'fwhm': 67.73201160365716, 'klz': 96.1696785929671, '2-sigma': 295.85924841600354}, 0.029: {'1-sigma': 3.227244202804286, '3-sigma': 91.78470542957535, 'fwhm': 4.545140399015509, 'klz': 28.837731118066934, '2-sigma': 20.473612865109317}, 0.02: {'1-sigma': 0.3891544728685403, '3-sigma': 11.404611063743879, 'fwhm': 0.5427439321938436, 'klz': 15.018502268045767, '2-sigma': 2.5010184130460402}}
fc = ROOT.TFeldmanCousins()
fc.SetMuMax(100)

resultsDir = '../results/done/'
resultsNameTemp = 'fits_hamamatsu_0nu_eff_tpc_batag_resol%.3f_5.0_years_0.0_counts_*.root'

ROOT.gSystem.Load('../lib/libnEXOSensitivity.so')

resolList = [0.005,0.010,0.015,0.025,0.029,0.047] #[0.003,0.005,0.008,0.010,0.012,0.015,0.017,0.020,0.025,0.030,0.040]
sensDict = {}

graph = ROOT.TGraph()
graphCount = ROOT.TGraph()

for resol in resolList:
    print 'working on resol', resol

    n = graph.GetN()
    x = resol
    x *= 100

    inFileName = '%s/%s'%(resultsDir,resultsNameTemp%(resol))
    sensDict[resol] = ROOT.nEXOUtils.GetBackgroundCounts(inFileName,"num_signal+num_signal_eHi") #ROOT.nEXOUtils.GetSensHalfLife(inFileName,5.0,3740)
    # sensDict[resol] /= (3740 * 5.)

    y = ROOT.nEXOUtils.GetHalfLife(sensDict[resol],5,3740) #1.4634,0.3445)
    #y = sensDict[resol] / 3.74 / 5.
    y /= 1e27
    graph.SetPoint(n,x,y)

    yct = nexo_batag_bkgd[resol]['fwhm']*3.74*5
    print yct
    fc.SetMuMax(yct*20)
    if yct < 50:
        avg = 0
        for i in range(int(yct - 5*yct**0.5),int(yct + 5*yct**0.5)):
            avg += fc.CalculateUpperLimit(i,yct) * ROOT.TMath.Poisson(i,yct)
        print sensDict[resol], avg
        yct = ROOT.nEXOUtils.GetHalfLife(avg,5,3740)
        yct /= 1e27
        graphCount.SetPoint(n,x,yct)
    
#graph.GetYaxis().SetTitle('Upper Limit Signal (cts/t.yr)')
graph.GetYaxis().SetTitle('Sensitivity (#times 10^{27} years)')
graph.GetXaxis().SetTitle('Resolution (%)')
graph.Draw('APL*')
graphCount.SetLineColor(2)
graphCount.SetMarkerColor(2)
graphCount.Draw('PL*')

raw_input('end')
    
