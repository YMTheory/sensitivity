
import ROOT

resultsDir = '../results/done/'
resultsNameTemp = 'fits_hamamatsu_0nu_eff_tpc_rdm_ssfrac%.1f_5.0_years_0.0_counts_*.root'

ROOT.gSystem.Load('../lib/libnEXOSensitivity.so')

resolList = [0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.3,1.5] #[0.5,0.6,0.7,1.0,1.3,1.5] #[0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5]
sensDict = {}

graph = ROOT.TGraph()
for resol in resolList:
    inFileName = '%s/%s'%(resultsDir,resultsNameTemp%(resol))
    sensDict[resol] = ROOT.nEXOUtils.GetSensHalfLife(inFileName,5.0,3740)
    
    x = resol
    y = sensDict[resol]
    n = graph.GetN()
    print x,y
    y /= 1e27
    graph.SetPoint(n,x,y)

graph.GetYaxis().SetTitle('Sensitivity (#times 10^{27} yrs)')
graph.GetXaxis().SetTitle('Improvement Fraction')
graph.Draw('APL*')
raw_input('end')
    
