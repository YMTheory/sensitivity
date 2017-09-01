
import ROOT
import copy

ROOT.gSystem.Load('../lib/libnEXOSensitivity.so')

def GetVolvsSD():
    gr = ROOT.TGraph()
    gr.SetPoint(0,90,3)
    gr.SetPoint(1,122,2.5)
    gr.SetPoint(2,159,2)
    gr.SetPoint(3,202,1.5)
    gr.SetPoint(4,256,1)
    gr.SetPoint(5,333,0.5)
    for i in range(5):
        gr.Fit('pol3','Q')
    pol = gr.GetFunction('pol3')
    gr.SetMarkerStyle(20)
    c = ROOT.TCanvas()
    gr.Draw('APL')
    #raw_input('continue')
    print pol
    return copy.copy(pol)

volf = GetVolvsSD()

def GetBkgdCounts(filename):

    chain = ROOT.TChain('tree')
    chain.Add(filename)
    chain.SetEstimate(chain.GetEntries()+1)
    chain.Draw('bkg_fwhm_3t','','goff')

    print 'using', chain.GetSelectedRows()
    
    return ROOT.TMath.Median(chain.GetSelectedRows(),chain.GetV1())


resultsNameTemp = '../results/done/fits_hamamatsu_v68_2016-06-21_0nu_sd%dm_fine_rdm_5.0_years_0.0_counts_*.root'

sdList = [10,100,130,160,190,220,250] #range(100,290,30) #+ [340]
#sdList = [10] + range(100,350,60)

nameList = []
effList = {}
volList = {}
for sd in sdList:
    nameList.append(resultsNameTemp % (sd))
    effList[nameList[-1]] = 0.82
    if sd == 10:
        volList[nameList[-1]] = 3.74
    else:
        volList[nameList[-1]] = volf.Eval(sd)

#nameList.append('../results/done/fits_hamamatsu_v61_2016-02-24_Si_Cu_0nu_tpc_cryo_elec_fine_rdm_5.0_years_0.0_counts_*.root')
#effList[nameList[-1]] = 0.82
#nameList.append('../results/done/fits_hamamatsu_v61_0nu_eff_tpc_cryo_elec_4main_rdm_5.0_years_0.0_counts_*.root')
#effList[nameList[-1]] = 1
#nameList.append('../results/done/fits_hamamatsu_v61_0nu_eff_tpc_cryo_elec_4main_further_rdm_5.0_years_0.0_counts_*.root')
#effList[nameList[-1]] = 1

graph = ROOT.TGraph()
#graph.SetPoint(0,0,28)

for sd, name in zip(sdList,nameList):
    print 'working on', name,
    
    n = graph.GetN()
    x = volList[name]
    y = ROOT.nEXOUtils.GetSensHalfLife(name,5.0,volList[name]*1000,'signal',effList[name])
    print x, y
    y /= 1e27
    graph.SetPoint(n,x,y)

# graph.SetPoint(graph.GetN(),3740,6.2)

#pf = ROOT.TF1('pf','[0]*(1+[1]*pow(x,[2]))',0.3,5)
#pf.SetParameters(1e27,1.5,0.3)
pf = ROOT.TF1('pf','[0]+[1]*x',1.5,5)
pf.SetParameters(1e27,1)
#for i in range(10):
#    graph.Fit(pf,'QR')

#print pf.Eval(3.74), pf.Eval(4.66), pf.Eval(4.78)

canvas = ROOT.TCanvas()
graph.SetMarkerStyle(20)
graph.SetMarkerColor(1)
graph.GetYaxis().SetTitle('Sensitivity @ 90% CL (#times 10^{27} years)')
graph.GetXaxis().SetTitle('Volume (tonne)')

graph.Draw('APL')

raw_input('end')
