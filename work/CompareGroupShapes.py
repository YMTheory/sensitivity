
import ROOT
import copy

libRelPath = '../lib/libnEXOSensitivity.so'
ROOT.gSystem.Load(libRelPath)

table = ROOT.TChain('ExcelTableValues')
table.Add('../tables/Summary_v61_2016-02-24_Si_Cu_0nu_tpc.root')

all_groups = ['LXeBb0n','LXeXe137','LXeBb2n','LXeRn222','FullTpcK40','FullTpcCo60','Far','InternalU238','InternalTh232','VesselTh232','VesselU238']
discards = ['LXeBb0n','LXeXe137','LXeBb2n','FullTpcK40','FullTpcCo60']

groups = []
comps = {}

groupPdf_ss = {}
groupPdf_ms = {}
compPdf_ss = {}
compPdf_ms = {}

rb = 20

nentries = table.GetEntries()
for ientry in range(nentries):
    table.GetEntry(ientry)

    group = table.fGroup.Data()

    if group in discards:
        continue

    comp = table.fPdf.Data()
    ul_ss = table.fCountsUL[0]
    ul_ms = table.fCountsUL[1]

    pdf = ROOT.TFile.Open(table.fFileName.Data())
    
    print comp, group, pdf

    if not group in comps:
        comps[group] = []

    if not comp in comps[group]:
        comps[group].append(comp)

    compPdf_ss[comp] = copy.copy(pdf.Get('h_SSEnergy_Smear'))
    compPdf_ss[comp].Sumw2()
    norm = compPdf_ss[comp].Integral()
    if norm > 0:
        compPdf_ss[comp].Scale(ul_ss/norm)
    compPdf_ms[comp] = copy.copy(pdf.Get('h_MSEnergy_Smear'))
    compPdf_ms[comp].Sumw2()
    norm = compPdf_ms[comp].Integral()
    if norm > 0:
        compPdf_ms[comp].Scale(ul_ms/norm)

    if not group in groups:
        groups.append(group)
        groupPdf_ss[group] = copy.copy(compPdf_ss[comp])
        groupPdf_ms[group] = copy.copy(compPdf_ms[comp])
    else:
        groupPdf_ss[group].Add(compPdf_ss[comp])
        groupPdf_ms[group].Add(compPdf_ms[comp])             

    pdf.Close()

cs = {}
legs = []

cs['ss'] = ROOT.TCanvas('ss')
cs['ms'] = ROOT.TCanvas('ms')
cs['rss'] = ROOT.TCanvas('rss')

c = cs['ss']
c.Divide(5-len(discards)/2,2)
for i in range(len(groups)):
    c.cd(i+1)
    histo = groupPdf_ss[groups[i]]
    histo.Rebin(rb)
    #histo.Scale(1/histo.Integral())
    histo.SetTitle('%s SS'%groups[i])
    histo.SetLineColor(1)
    histo.Draw('hist')
    legs.append(ROOT.TLegend(0.1,0.7,0.48,0.9))
    leg = legs[-1]
    for k,comp in enumerate(comps[groups[i]]):
        h = compPdf_ss[comp]
        h.Rebin(rb)
        #h.SetMarkerColor(k+2)
        #h.Scale(histo.Integral()/h.Integral())
        #h.Divide(histo)
        h.SetLineColor(k+2)
        h.Draw('hist same')
        name = comp[:comp.rfind('_')]
        leg.AddEntry(h,name,'L')
    leg.Draw()
        

c = cs['ms']
c.Divide(5-len(discards)/2,2)
for i in range(len(groups)):
    c.cd(i+1)
    histo = groupPdf_ms[groups[i]]
    histo.Rebin(rb)
    #histo.Scale(1/histo.Integral())
    histo.SetTitle('%s MS'%groups[i])
    histo.SetLineColor(1)
    histo.Draw('hist')
    legs.append(ROOT.TLegend(0.1,0.7,0.48,0.9))
    leg = legs[-1]
    for k,comp in enumerate(comps[groups[i]]):
        h = compPdf_ms[comp]
        h.Rebin(rb)
        #h.SetMarkerColor(k+2)
        h.SetLineColor(k+2)
        h.Draw('hist same')
        name = comp[:comp.rfind('_')]
        leg.AddEntry(h,name,'L')
    leg.Draw()

c = cs['rss']
c.Divide(5-len(discards)/2,2)
for i in range(len(groups)):
    c.cd(i+1)
    histo = groupPdf_ss[groups[i]]
    legs.append(ROOT.TLegend(0.1,0.7,0.48,0.9))
    leg = legs[-1]
    for k,comp in enumerate(comps[groups[i]]):
        h = compPdf_ss[comp]
        h.SetMarkerColor(k+2)
        h.Scale(histo.Integral()/h.Integral())
        h.Divide(histo)
        if k == 0:
            h.SetTitle('%s SS'%groups[i])
            h.GetYaxis().SetRangeUser(0,5)
            h.Draw()
        else:
            h.Draw('same')
        name = comp[:comp.rfind('_')]
        leg.AddEntry(h,name,'L')
    leg.Draw()


raw_input('end')

