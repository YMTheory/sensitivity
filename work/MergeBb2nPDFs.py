
import ROOT
import copy
import SmearHisto as sh

# volume used for normalization (90 = inner 3t)
sdcut = 90
elo = 2200
ehi = 2300
resqcorr = 0.5748 # scaling was 0.5988 -> 0.0104 resol, but 0.5748 -> 0.0100

# resolution function
resol = ROOT.TF1('resol',"%f*TMath::Sqrt([0]*[0]*x + [1]*[1] + [2]*[2]*x*x)"%(resqcorr),0,100000)
resol.SetParameters(0.,36.9,8.8e-3)
print resol.Eval(2458)

# read in files
files = {}
histos = {}
htemp = {}

files['nls'] =  ROOT.TFile.Open('../histos/bb2n/SensBb2nDef_FullLXe.root')
htemp['nls_ss'] = files['nls'].Get('h_StandoffVsEnergySS')
htemp['nls_ms'] = files['nls'].Get('h_StandoffVsEnergyMS')
histos['nls'] = copy.copy(htemp['nls_ss'])
histos['nls'].Add(htemp['nls_ms'])

files['nhs'] = ROOT.TFile.Open('../histos/bb2n/PlotsSensEndPointStatNew_FullLXe.root')
htemp['nhs_ss'] = files['nhs'].Get('h_StandoffVsEnergySS')
htemp['nhs_ms'] = files['nhs'].Get('h_StandoffVsEnergyMS')
histos['nhs'] = copy.copy(htemp['nhs_ss'])
histos['nhs'].Add(htemp['nhs_ms'])

# use projections onto energy for normalization
projs = {}
for name in histos:
    histo = histos[name]
    print name, histo.GetNbinsX(), histo.GetNbinsY(), histo.GetXaxis().GetBinCenter(1),  histo.GetYaxis().GetBinCenter(37), histo.Integral(1,2800,1,37), histo.Integral(1,2800,37,91), histo.Integral(1,2800,91,123), histo.Integral(1,2800,90,123), histo.Integral(1,-1,1,-1), histo.Integral(1,2800,123,160), histo.Integral(1,2800,160,203)
    projs[name] = histo.ProjectionX(name+'_ene',sdcut,-1)

# combine bb2n spectra (low + high stat end point)
normlo = projs['nls'].FindBin(elo)
normhi = projs['nls'].FindBin(ehi)

il = projs['nls'].Integral(normlo,normhi)
ih = projs['nhs'].Integral(normlo,normhi)
norm = il*1./ih
print 'Using normalization', norm

# create file with combined histograms
files['ncs'] = ROOT.TFile.Open('../histos/bb2n/nEXO_Histos_bb2n.root','recreate')

for hname in ['h_SSEnergy','h_MSEnergy']:
    print 'working on', hname

    hls = files['nls'].Get(hname)
    hhs = files['nhs'].Get(hname)

    hcs = copy.copy(hhs)
    hcs.Scale(norm)
    for b in range(1,hcs.GetXaxis().GetNbins()+1):
        if hcs.GetBinCenter(b) < elo:
            hcs.SetBinContent(b,hls.GetBinContent(b))

    files['ncs'].cd()
    hcs.Write()

    shcs = sh.SmearHisto1D(hcs,resol,1)
    shcs_name = shcs.GetName()
    shcs.SetName(shcs_name.replace('GausSmear','Smear'))
    shcs.Write()

for hname in ['h_StandoffDistSS','h_StandoffDistMS']:
    print 'working on', hname

    hls = copy.copy(files['nls'].Get(hname))
    files['ncs'].cd()
    hls.Write()

for hname in ['h_StandoffVsEnergySS','h_StandoffVsEnergyMS']:
    print 'working on', hname

    hls = files['nls'].Get(hname)
    hhs = files['nhs'].Get(hname)

    hcs = copy.copy(hhs)
    hcs.Scale(norm)
    for bb in range(1,hcs.GetYaxis().GetNbins()+1):
        for b in range(1,hcs.GetXaxis().GetNbins()+1):
            if hcs.GetXaxis().GetBinCenter(b) < elo:
                hcs.SetBinContent(b,bb,hls.GetBinContent(b,bb))

    files['ncs'].cd()
    hcs.Write()

    hcs.Rebin2D(1,10)
    shcs = sh.SmearHisto2D(hcs,resol,10)
    shcs_name = shcs.GetName()
    shcs.SetName(shcs_name.replace('GausSmear','Smear'))
    shcs.Write()

# close files
files['ncs'].Close()
files['nhs'].Close()
files['nls'].Close()
