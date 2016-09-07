
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

    hss = pdf.Get('h_StandoffVsEnergySS_Smear')
    hss.Sumw2()
    nss = hss.Integral()
    if nss > 0:
        hss.Scale(ul_ss/nss)
    compPdf_ss[comp] = copy.copy(hss.ProjectionX('_px',1))
    hms = pdf.Get('h_StandoffVsEnergyMS_Smear')
    hms.Sumw2()
    nms = hms.Integral()
    if nms > 0:
        hms.Scale(ul_ms/nms)
    compPdf_ms[comp] = copy.copy(hms.ProjectionX('_px',1))


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

c = cs['ss']
c.Divide(5-len(discards)/2,2)
for i in range(len(groups)):
    c.cd(i+1)
    histo = groupPdf_ss[groups[i]]
    #histo.Scale(1/histo.Integral())
    histo.SetTitle('%s SS'%groups[i])
    histo.SetLineColor(1)
    histo.Draw('hist')
    legs.append(ROOT.TLegend(0.1,0.7,0.48,0.9))
    leg = legs[-1]
    for k,comp in enumerate(comps[groups[i]]):
        h = compPdf_ss[comp]
        #h.SetMarkerColor(k+2)
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
    #histo.Scale(1/histo.Integral())
    histo.SetTitle('%s MS'%groups[i])
    histo.SetLineColor(1)
    histo.Draw('hist')
    legs.append(ROOT.TLegend(0.1,0.7,0.48,0.9))
    leg = legs[-1]
    for k,comp in enumerate(comps[groups[i]]):
        h = compPdf_ms[comp]
        #h.SetMarkerColor(k+2)
        h.SetLineColor(k+2)
        h.Draw('hist same')
        name = comp[:comp.rfind('_')]
        leg.AddEntry(h,name,'L')
    leg.Draw()


raw_input('end')

