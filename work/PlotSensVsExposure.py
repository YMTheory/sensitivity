
import ROOT
import copy

ROOT.gSystem.Load('../lib/libnEXOSensitivity.so')

ryans = {0.5:1.388,1:2.203,2:3.429,3:4.408,4:5.246,5:5.997}
dbv68s = {0.5:1.35,1:2.167,2.5:3.957,5:6.149,10:9.4}
reds = {0.5:1.085,1:1.694,2.5:3.122,5:4.831,10:7.119}

gpt = ROOT.TGraph()
gpe = ROOT.TGraph()
gpa = ROOT.TGraph()
for t in ryans:
    gpt.SetPoint(gpt.GetN(),t,ryans[t])
    gpe.SetPoint(gpe.GetN(),t*4.66*0.9,ryans[t])
    gpa.SetPoint(gpe.GetN(),t*4.66*0.9,ryans[t]/0.9)

gpt.SetMarkerStyle(25)
gpt.SetMarkerColor(4)
gpe.SetMarkerStyle(25)
gpe.SetMarkerColor(4)
gpa.SetMarkerStyle(25)
gpa.SetMarkerColor(4)

gct = ROOT.TGraph()
gce = ROOT.TGraph()
gca = ROOT.TGraph()
for t in dbv68s:
    gct.SetPoint(gct.GetN(),t,dbv68s[t])
    gce.SetPoint(gce.GetN(),t*3.74*0.9,dbv68s[t])
    gca.SetPoint(gce.GetN(),t*3.74*0.9,dbv68s[t]/0.813)
gct.SetMarkerStyle(24)
gct.SetMarkerColor(1)
gce.SetMarkerStyle(24)
gce.SetMarkerColor(1)
gca.SetMarkerStyle(24)
gca.SetMarkerColor(1)

grt = ROOT.TGraph()
gre = ROOT.TGraph()
gra = ROOT.TGraph()
for t in reds:
    grt.SetPoint(grt.GetN(),t,reds[t])
    gre.SetPoint(gre.GetN(),t*3.74*0.9,reds[t])
    gra.SetPoint(gre.GetN(),t*3.74*0.9,reds[t]/0.82)
grt.SetMarkerStyle(26)
grt.SetMarkerColor(3)
gre.SetMarkerStyle(26)
gre.SetMarkerColor(3)
gra.SetMarkerStyle(26)
gra.SetMarkerColor(3)

canvas = ROOT.TCanvas()
leg = ROOT.TLegend(0.1,0.7,0.48,0.9)

canvas.Divide(3,1)

canvas.cd(1)
gct.GetXaxis().SetTitle('Detector Livetime (yrs)')
gct.GetYaxis().SetTitle('Sensitivity @ 90% CL (yrs)')
gct.Draw('AP')
leg.AddEntry(gct,'Current Model','P')
gpt.Draw('P')
leg.AddEntry(gpt,'Old Model','P')
grt.Draw('P')
leg.AddEntry(grt,'Current Model w/ Bkgd = Old','P')
leg.Draw()

canvas.cd(2)
gce.GetXaxis().SetTitle('Exposure (t.yrs)')
gce.GetYaxis().SetTitle('Sensitivity @ 90% CL (yrs)')
gce.Draw('AP')
gpe.Draw('P')
gre.Draw('P')

canvas.cd(3)
gca.GetXaxis().SetTitle('Exposure (t.yrs)')
gca.GetYaxis().SetTitle('Sensitivity / Signal Efficiency @ 90% CL (yrs)')
gca.Draw('AP')
gpa.Draw('P')
gra.Draw('P')

raw_input('end')
