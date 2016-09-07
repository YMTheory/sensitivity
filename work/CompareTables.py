
import ROOT

ROOT.gSystem.Load('../lib/libnEXOSensitivity.so')

fl = ROOT.TFile.Open('../tables/Summary_01042016_Si_Cu_v2.root')
tl = fl.Get('ExcelTableValues')
fr = ROOT.TFile.Open('../tables/Summary_01042016_original_v2.root')
tr = fr.Get('ExcelTableValues')

#table = ROOT.ExcelTableValues()

nl = tl.GetEntries()
nr = tr.GetEntries()

print nl, nr

n = nl

for i in range(n):
    tl.GetEntry(i)
    tr.GetEntry(i)

    print tl.fPdf, tl.fCountsCV[0], tr.fPdf, tr.fCountsCV[0]
    print tl.fPdf, tl.fCountsCV[1], tr.fPdf, tr.fCountsCV[1]
