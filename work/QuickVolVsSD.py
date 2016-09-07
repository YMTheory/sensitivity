
import ROOT

gr = ROOT.TGraph()

gr.SetPoint(0,90,3)
gr.SetPoint(1,122,2.5)
gr.SetPoint(2,159,2)
gr.SetPoint(3,202,1.5)
gr.SetPoint(4,256,1)
gr.SetPoint(5,333,0.5)

gr.SetMarkerStyle(20)
c = ROOT.TCanvas()
gr.Draw('APL')

raw_input('end')
