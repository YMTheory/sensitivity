
import ROOT

resultsDir = '../results/done/'
resultsNameTemp = 'fits_hamamatsu_0nu_eff_tpc_rdm_thres%dkeV_5.0_years_0.0_counts_*.root'

ROOT.gSystem.Load('../lib/libnEXOSensitivity.so')

thresList = [700,800,900,1000]
sensDict = {}

graph = ROOT.TGraph()
for thres in thresList:
    inFileName = '%s/%s'%(resultsDir,resultsNameTemp%(thres))
    sensDict[thres] = ROOT.nEXOUtils.GetSensHalfLife(inFileName,5.0,3740)
    
    x = thres
    y = sensDict[thres]
    n = graph.GetN()
    print x,y
    y /= 1e27
    graph.SetPoint(n,x,y)

graph.GetYaxis().SetTitle('Sensitivity (#times 10^{27} yrs)')
graph.GetXaxis().SetTitle('Energy Threshold Cut (keV)')
graph.Draw('APL*')
raw_input('end')
    
