
import ROOT

gr = ROOT.TGraph()

gr.SetPoint(0,0,28)
#gr.SetPoint(1,1.39,8.1)
gr.SetPoint(1,1.82,7.3)
gr.SetPoint(2,3.63,6.3)
gr.SetPoint(3,15.0,4.3)
gr.SetPoint(4,79.5,2.9)


gr.SetMarkerStyle(20)
c = ROOT.TCanvas()
gr.Draw('APL')

raw_input('end')
