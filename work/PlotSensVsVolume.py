
import ROOT

resultsDir = '../results/done/'
resultsNameTemp = 'fits_hamamatsu_0nu_eff_tpc_rdm_%s_5.0_years_0.0_counts_*.root'

ROOT.gSystem.Load('../lib/libnEXOSensitivity.so')

resolList = ['fv','3t','1t']
volDict = {'fv':3740,'3t':3000,'1t':1000}
sensDict = {}

graph = ROOT.TGraph()
for r,resol in enumerate(resolList):
    print r, resol
    inFileName = '%s/%s'%(resultsDir,resultsNameTemp%(resol))
    sensDict[resol] = ROOT.nEXOUtils.GetSensHalfLife(inFileName,5.0,3740)#volDict[resol])
    
    x = r+0.5# resol
    y = sensDict[resol]
    n = graph.GetN()
    print x,y
    y /= 1e27
    graph.SetPoint(n,x,y)

xaxis = graph.GetHistogram().GetXaxis()
xaxis.Set(3,xaxis.GetBinLowEdge(1),xaxis.GetBinUpEdge(xaxis.GetNbins()))

for r,resol in enumerate(resolList):
    xaxis.SetBinLabel(r+1,resol)

graph.GetYaxis().SetTitle('Sensitivity (#times 10^{27} yrs)')
graph.GetXaxis().SetTitle('LXe Volume')
graph.Draw('APL*')
raw_input('end')
    
