
import ROOT
import sys

#ROOT.gInterpreter.GenerateDictionary("std::map<TString, TH1*>","map;TString.h;TH1.h")

ROOT.gSystem.Load('../lib/libnEXOSensitivity.so')

sens = ROOT.nEXOSensitivity()

sens.fCreateWsp = False
sens.GenAndFitData(1)

sys.exit(1)

sens.LoadExcelTree('../tables/Summary_01042016_Si_Cu_v2.root')
sens.ReadExcelTree()

#sens.fExcelTree.Draw("fGroup.Data()","","goff")
#print sens.fExcelTree.GetVar1().GetNdata()
#raw_input('end')

sens.fHistoFileName = "../histos/combined/QuickTestNew2_nEXO_Histograms_01042016_Si_Cu_v2.root"
sens.LoadComponentHistograms()
sens.fComponentHistos["h_LXe_bb2n_ss"].Draw("colz")
#raw_input('continue')

sens.MakeGroupHistograms()
sens.fGroupHistos["h_FullTpcCo60_ss"].Draw("colz")

#raw_input('end')

sens.fWriteWsp = True
sens.BuildWorkspace()
sens.fWsp.Print()

#sens.GenAndFitData() #Int_t nRuns = 1, Double_t yrs = 5., Double_t signalCounts = 0., Int_t fileNum = 0, Int_t initialSeed = 0, Bool_t runMinos=true, Bool_t runBkgdOnlyFit=false, Bool_t withStandoff = true);
