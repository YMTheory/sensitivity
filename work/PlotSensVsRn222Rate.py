
import ROOT

resultsDir = '../results/done/'
resultsNameTemp = 'fits_hamamatsu_0nu_eff_tpc_rdm_%sx_rn222_5.0_years_0.0_counts_*.root'

ROOT.gSystem.Load('../lib/libnEXOSensitivity.so')

resolList = [0.17,0.33,1,3,9,12]
sensDict = {}

rn222Atoms = {0.17: 100, 0.33: 200, 1: 600, 3: 1800, 9: 5400, 12: 7200}

graph = ROOT.TGraph()
for resol in resolList:
    inFileName = '%s/%s'%(resultsDir,resultsNameTemp%(str(resol)))
    sensDict[resol] = ROOT.nEXOUtils.GetSensHalfLife(inFileName,5.0,3740)
    
    x = rn222Atoms[resol]
    y = sensDict[resol]
    y /= 1e27
    n = graph.GetN()
    print x,y
    graph.SetPoint(n,x,y)

graph.GetYaxis().SetTitle('Sensitivity (#times 10^{27} yrs)')
graph.GetXaxis().SetTitle('Rn222 (atoms)')
graph.Draw('APC*')
raw_input('end')
    
