
import ROOT

resultsDir = '../results/done/'
resultsNameTemp = 'fits_hamamatsu_0nu_eff_tpc_rdm_sd%d_5.0_years_0.0_counts_*.root'

ROOT.gSystem.Load('../lib/libnEXOSensitivity.so')

resolList = [10,25,50,65]
sensDict = {}

graph = ROOT.TGraph()
for resol in resolList:
    inFileName = '%s/%s'%(resultsDir,resultsNameTemp%(resol))
    if resol == 65:
        inFileName = inFileName.replace('sd65','fv')
    sensDict[resol] = ROOT.nEXOUtils.GetSensHalfLife(inFileName,5.0,3740)
    
    x = resol
    y = sensDict[resol]
    y /= 1e27
    n = graph.GetN()
    print x,y
    graph.SetPoint(n,x,y)

graph.GetYaxis().SetTitle('Sensitivity (#times 10^{27} yrs)')
graph.GetXaxis().SetTitle('Number of standoff-distance bins')
graph.Draw('APL*')
raw_input('end')
    
