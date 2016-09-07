
import ROOT

import SmearHisto as sh

rootFileInName = '../histos/individual_histos/nEXO_Histos_TPC_Co60.root'#'/data/data033/exo/data/nEXO/licciard/CompareBb2n/data/SensBb2nDef_FullLXe.root'
#globalArgs['ResolSS'] = ROOT.TF1('resol_ss',card.Get('ResolSS'),0,100000)
#globalArgs['ResolSS'].SetParameters(array('d',[float(i) for i in card.Get('ResolPss').split(',')]))
resol = ROOT.TF1('resol',"%f*TMath::Sqrt([0]*[0]*x + [1]*[1] + [2]*[2]*x*x)"%(0.5988),0,100000)
resol.SetParameters(0.,36.9,8.8e-3)

rootFileIn = ROOT.TFile.Open(rootFileInName)
#histoIn = rootFileIn.Get('h_SSEnergy')
histoIn = rootFileIn.Get('h_StandoffVsEnergySS')
histoIn.Rebin2D(1,10)
#histoOut = SmearHisto1D(histoIn, resol, 14)
histoOut = sh.SmearHisto2D(histoIn, resol, 10)
histoOut.Draw()

raw_input('end')

